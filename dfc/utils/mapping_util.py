#! /usr/bin/env python
# coding: UTF8

"""Helpers for mapping nucleotide references (CARD/VFDB) onto the genome with blastn."""

import re
from dataclasses import dataclass

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
