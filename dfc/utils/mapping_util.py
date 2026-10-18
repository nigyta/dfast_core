#! /usr/bin/env python
# coding: UTF8

"""Helpers for mapping nucleotide references (CARD/VFDB) onto the genome with blastn."""

import re
from dataclasses import dataclass, field

from Bio.Data import CodonTable
from Bio.Seq import Seq

BLAST_OUTFMT = "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen btop"
END_SLACK = 3  # bases of the reference that may be left unaligned at each end of a full-length hit
STOP_EXTENSION = 100  # bases searched downstream of a region for an in-frame stop codon
TAIL_STOP_CODONS = 2  # number of codons at the 3' end that may be lost without making a CDS truncated

BTOP_ITEM = re.compile(r"\d+|[A-Z-]{2}")


@dataclass
class MappingHit:
    ref_id: str
    seq_id: str
    pident: float
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int
    btop: str

    @property
    def strand(self) -> int:
        return 1 if self.sstart <= self.send else -1

    @property
    def left0(self) -> int:
        """0-based start of the region extrapolated to the full reference."""
        if self.strand == 1:
            return self.sstart - 1 - (self.qstart - 1)
        return self.send - 1 - (self.qlen - self.qend)

    @property
    def right(self) -> int:
        """Exclusive end of the region extrapolated to the full reference."""
        if self.strand == 1:
            return self.send + (self.qlen - self.qend)
        return self.sstart + (self.qstart - 1)


def read_blast_table(file_name: str) -> list[MappingHit]:
    hits = []
    with open(file_name) as f:
        for line in f:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            hits.append(MappingHit(ref_id=c[0], seq_id=c[1], pident=float(c[2]),
                                   qstart=int(c[6]), qend=int(c[7]), sstart=int(c[8]), send=int(c[9]),
                                   evalue=float(c[10]), bitscore=float(c[11]),
                                   qlen=int(c[12]), slen=int(c[13]), btop=c[14]))
    return hits


def is_full_length(hit: MappingHit, pident_cutoff: float) -> bool:
    return (hit.qstart <= 1 + END_SLACK and hit.qend >= hit.qlen - END_SLACK
            and hit.pident >= pident_cutoff and 0 <= hit.left0 < hit.right <= hit.slen)


def indel_analysis(hit: MappingHit) -> tuple[int, list[int]]:
    """Return (shift, positions) from the btop string.

    shift is the bases inserted in the genome minus the bases deleted from it, relative to the reference.
    positions are the 1-based genome coordinates where each run of gaps starts.
    """
    shift = 0
    positions = []
    offset = 0  # subject bases consumed so far
    prev_kind = None
    for item in BTOP_ITEM.findall(hit.btop):
        if item.isdigit():
            offset += int(item)
            prev_kind = None
            continue
        if item[1] == "-":  # deletion: base missing in the genome
            kind, step, delta = "del", 0, -1
        elif item[0] == "-":  # insertion: extra base in the genome
            kind, step, delta = "ins", 1, 1
        else:  # mismatch
            offset += 1
            prev_kind = None
            continue
        if kind != prev_kind:
            positions.append(hit.sstart + offset if hit.strand == 1 else hit.sstart - offset)
        shift += delta
        offset += step
        prev_kind = kind
    return shift, positions


def group_by_locus(hits: list[MappingHit]) -> list[list[MappingHit]]:
    """Group hits whose extrapolated regions overlap on the same sequence and strand.

    Each group is sorted by (bitscore, pident) descending; groups are ordered by (seq_id, left0).
    """
    groups = []
    right = 0  # largest right end in the current group
    for h in sorted(hits, key=lambda h: (h.seq_id, h.strand, h.left0)):
        last = groups[-1][0] if groups else None
        if last and (last.seq_id, last.strand) == (h.seq_id, h.strand) and h.left0 < right:
            groups[-1].append(h)
            right = max(right, h.right)
        else:
            groups.append([h])
            right = h.right
    groups.sort(key=lambda g: (g[0].seq_id, min(h.left0 for h in g)))
    for g in groups:
        g.sort(key=lambda h: (h.bitscore, h.pident), reverse=True)
    return groups


def is_valid_reference(nucl_seq: str, transl_table: int) -> tuple[bool, bool]:
    """Return (valid, has_terminal_stop) for a reference CDS.

    valid: a multiple of 3, only ACGT, starts with a start codon and has no internal stop codon.
    has_terminal_stop: the last codon is a stop codon.
    """
    table = CodonTable.unambiguous_dna_by_id[transl_table]
    seq = nucl_seq.upper()
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    has_terminal_stop = len(seq) % 3 == 0 and bool(codons) and codons[-1] in table.stop_codons
    body = codons[:-1] if has_terminal_stop else codons
    valid = (len(seq) % 3 == 0 and bool(codons) and set(seq) <= set("ACGT")
             and codons[0] in table.start_codons
             and not any(c in table.stop_codons for c in body))
    return valid, has_terminal_stop


@dataclass
class MappedRegion:
    seq_id: str
    left0: int
    right: int
    strand: int
    kind: str  # "valid", "no_start", "frameshift", "internal_stop" or "stop_lost"
    translation: str = ""  # protein without the stop; empty for frameshift, internal_stop and stop_lost
    notes: list[str] = field(default_factory=list)
    stop_codons: set[tuple[int, int, int]] = field(default_factory=set)  # (start0, end, strand) in genome coordinates
    indels: set[int] = field(default_factory=set)  # 1-based genome positions where gap runs start


def resolve_region(record_seq: Seq, hit: MappingHit, transl_table: int, ref_has_stop: bool) -> MappedRegion | None:
    """Classify the genome region covered by a full-length hit as a valid CDS or one of the defect kinds.

    Return None for a reference without a terminal stop whose region has no in-frame stop within STOP_EXTENSION.
    """
    table = CodonTable.unambiguous_dna_by_id[transl_table]
    left0, right, strand = hit.left0, hit.right, hit.strand
    region = MappedRegion(hit.seq_id, left0, right, strand, "valid")

    shift, positions = indel_analysis(hit)
    if shift % 3 != 0:
        region.kind = "frameshift"
        region.notes = [f"frameshift at around {p}" for p in positions]
        region.indels = set(positions)
        return region

    # Work in gene orientation; the sequence beyond the region (downstream) is kept for stop codon extension.
    length = right - left0
    if strand == 1:
        seq = str(record_seq[left0:right + STOP_EXTENSION]).upper()
    else:
        beyond = max(left0 - STOP_EXTENSION, 0)
        seq = str(record_seq[beyond:right].reverse_complement()).upper()

    def set_length(new_length):
        region.left0, region.right = (left0, left0 + new_length) if strand == 1 else (right - new_length, right)

    # Scan the first in-frame stop inside the region before any extension.
    codons = [seq[i:i + 3] for i in range(0, length, 3)]
    n_sense = len(codons) - 1 if ref_has_stop else len(codons)  # codons of the reference before its stop
    stops = [i for i, c in enumerate(codons) if c in table.stop_codons]
    if not stops:  # no stop in the region: look for one downstream
        stop = next((i for i in range(length, len(seq) - 2, 3)
                     if seq[i:i + 3] in table.stop_codons and i + 3 - length <= STOP_EXTENSION), None)
        if stop is None and ref_has_stop:
            region.kind = "stop_lost"
            region.notes.append("stop codon lost")
            return region
        if stop is None:
            return None
        if ref_has_stop:
            region.notes.append(f"stop codon of the reference is not conserved; extended by {stop + 3 - length} bp")
        codons = [seq[i:i + 3] for i in range(0, stop, 3)]
        set_length(stop + 3)
    elif stops[0] < n_sense:  # a stop before the reference's own stop
        if stops[0] > 0 and stops[0] >= n_sense - TAIL_STOP_CODONS:
            region.notes.append(f"{n_sense - stops[0]} codon(s) shorter than the reference")
            codons = codons[:stops[0]]
            set_length(3 * (stops[0] + 1))
        else:
            region.kind = "internal_stop"
            for i in (i for i in stops if i < n_sense):
                start0 = region.left0 + 3 * i if strand == 1 else region.right - 3 * i - 3
                region.stop_codons.add((start0, start0 + 3, strand))
                region.notes.append(f"internal stop codon at {start0 + 1}")
            return region
    else:  # the stop of the reference is conserved
        codons = codons[:n_sense]

    translation = str(Seq("".join(codons)).translate(table=transl_table))
    if codons[0] not in table.start_codons:
        region.kind = "no_start"
        region.notes.append("no start codon")
    else:
        translation = "M" + translation[1:]
    region.translation = translation
    return region


def classify_cds(cds, region: MappedRegion, tolerance: int) -> str | None:
    """Relate a predicted CDS on the same contig to a mapped region.

    Returns "same_gene" (same strand and same stop coordinate), "overlapping" (overlap longer than tolerance),
    "neighbour" (overlap of 1..tolerance bases) or None (no overlap).
    """
    start, end = int(cds.location.start), int(cds.location.end)
    if cds.location.strand == region.strand and (end == region.right if region.strand == 1 else start == region.left0):
        return "same_gene"
    overlap = min(end, region.right) - max(start, region.left0)
    if overlap <= 0:
        return None
    return "overlapping" if overlap > tolerance else "neighbour"
