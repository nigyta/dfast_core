#! /usr/bin/env python
# coding: UTF8

import os
import json
from Bio.SeqFeature import (FeatureLocation, ExactPosition,
                            BeforePosition, AfterPosition)
from ..models.bio_feature import ExtendedFeature
from .base_tools import ContigAnnotationTool

# JSON short type code -> INSDC /mobile_element_type value (registered as mobile_element)
_MOBILE_ELEMENT_TYPE = {
    "is": "insertion sequence",
    "mite": "MITE",
    "tn": "transposon",            # unit transposon
    "cn": "transposon",            # composite transposon (DB hit)
    "integron": "integron",
    "retrotransposon": "retrotransposon",
}

# JSON value of PredictionEvidence.PUTATIVE (corresponds to `r["evidence"] == 2` in io.py)
PUTATIVE_EVIDENCE = 2

# JSON short type code -> human-readable full name (for report/amr_summary)
_TYPE_NAME = {
    "is": "insertion sequence",
    "mite": "MITE",
    "tn": "unit transposon",
    "cn": "composite transposon",
    "ice": "integrative conjugative element",
    "aice": "actinomycete integrative conjugative element",
    "ime": "integrative mobilizable element",
    "cime": "cis-mobilizable element",
    "mic": "mobile insertion cassette",
    "iscr": "IS common region",
    "integron": "integron",
    "retrotransposon": "retrotransposon",
    "other": "mobile genetic element",
}


def classify_mge(type_code, evidence):
    """Return (feature_key, mobile_element_type, is_putative) from a MEF short type code and evidence.

    - Putative composite transposons (type=="cn" and evidence==2) are misc_feature.
    - Types with a clear INSDC /mobile_element_type counterpart are mobile_element.
    - Everything else (ICE/IME/CIME/MIC etc., or unknown) is misc_feature.
    """
    is_putative = (type_code == "cn" and evidence == PUTATIVE_EVIDENCE)
    if is_putative:
        return ("misc_feature", None, True)
    met = _MOBILE_ELEMENT_TYPE.get(type_code)
    if met is not None:
        return ("mobile_element", met, False)
    return ("misc_feature", None, False)


# Descriptive labels for types registered as misc_feature
_MISC_LABEL = {
    "ice": "integrative conjugative element (ICE)",
    "aice": "actinomycete integrative conjugative element (AICE)",
    "ime": "integrative mobilizable element (IME)",
    "cime": "cis-mobilizable element (CIME)",
    "mic": "mobile insertion cassette",
    "iscr": "IS common region (ISCR)",
    "other": "mobile genetic element",
}


def build_qualifiers(entry):
    """Build GenBank qualifiers (dict) from a result entry of the MEF JSON."""
    name = entry["name"]
    type_code = entry["type"]
    evidence = entry.get("evidence", 1)
    feature_key, met, is_putative = classify_mge(type_code, evidence)

    identity = float(entry.get("identity", 0)) * 100
    coverage = float(entry.get("coverage", 0)) * 100
    accession = (entry.get("template") or {}).get("accession", "")

    detail = "MobileElementFinder: {name}; identity:{id:.1f}%, coverage:{cov:.1f}%".format(
        name=name, id=identity, cov=coverage)
    if accession:
        detail += "; similar to {acc} (MGEdb)".format(acc=accession)
    notes = [detail]
    if is_putative:
        notes.append("putative composite transposon (predicted by MobileElementFinder)")
    elif feature_key == "misc_feature":
        label = _MISC_LABEL.get(type_code, "mobile genetic element")
        notes.append("{label}: {name}".format(label=label, name=name))

    qualifiers = {"note": notes}
    if met is not None:
        qualifiers["mobile_element_type"] = ["{met}:{name}".format(met=met, name=name)]
    return qualifiers


def _zero_based_left(start, end, allele_seq_length):
    """Return the 0-based left coordinate for biopython from the MEF start/end.

    MobileElementFinder's coordinate system is inconsistent. Putative composite transposons
    (`cn` evidence=2) report `start` as **0-based** (internally it stores the slice index
    `flank_IS_start - 1`, used to cut out the sequence, as start). IS elements and
    DB-hit composites report `start` as **1-based**. So only composites have a start that is
    1 too small, and at the beginning of a contig start=0 yields the invalid coordinate `<0`.

    The coordinate system is decided not by type or evidence but by consistency with
    `allele_seq_length`, which MEF itself reports (the true 1-based length `end - start + 1`,
    computed independently of the buggy start):
        end - start + 1 == allele_seq_length  -> start is 1-based -> left = start - 1
        end - start     == allele_seq_length  -> start is 0-based -> left = start
    If MEF later fixes composites to be 1-based, those entries are automatically classified
    as 1-based and get the `-1`, so no double correction happens.
    If the length is missing or matches neither, the documented 1-based start is assumed.
    """
    start, end = int(start), int(end)
    if allele_seq_length is not None:
        asl = int(allele_seq_length)
        if end - start == asl:          # 0-based start (MEF composite, already -1)
            left = start
        elif end - start + 1 == asl:    # 1-based start (normal / after an upstream fix)
            left = start - 1
        else:
            left = start - 1            # inconsistent -> assume the documented 1-based start
    else:
        left = start - 1
    return left


def _location(start, end, strand, trunc_5p, trunc_3p, allele_seq_length=None):
    """Build a FeatureLocation from the MEF start/end, strand and truncation.

    biopython uses 0-based half-open coordinates. _zero_based_left() absorbs the coordinate
    system difference for left. Truncation is mapped, taking the strand into account,
    from 5'/3' to partial (Before/After) left/right genomic ends.

    MEF's trunc_5p is the alignment start position on the reference (1-based), 1 when
    complete to the 5' end. trunc_3p is the number of unaligned bases on the 3' side,
    0 when complete. So 5' truncated means trunc_5p > 1 and 3' truncated means trunc_3p > 0.
    """
    left = _zero_based_left(start, end, allele_seq_length)
    right = int(end)
    is_5p_trunc = int(trunc_5p) > 1
    is_3p_trunc = int(trunc_3p) > 0
    # strand=+1: left end=5', right end=3' / strand=-1: left end=3', right end=5'
    if strand == -1:
        left_trunc, right_trunc = is_3p_trunc, is_5p_trunc
    else:
        left_trunc, right_trunc = is_5p_trunc, is_3p_trunc
    # If it goes past the contig start (left<0), clamp to 0 and treat it as partial.
    # Guards against passing an invalid `<0` to biopython (1.83 turns it into location=None).
    if left < 0:
        left = 0
        left_trunc = True
    left_pos = BeforePosition(left) if left_trunc else ExactPosition(left)
    right_pos = AfterPosition(right) if right_trunc else ExactPosition(right)
    return FeatureLocation(left_pos, right_pos, strand=strand)


def entry_to_feature(entry, index):
    """Convert one result entry of the MEF JSON to (seq_id, ExtendedFeature)."""
    seq_id = entry["contig"].split()[0]
    feature_key, _met, is_putative = classify_mge(entry["type"], entry.get("evidence", 1))
    location = _location(entry["start"], entry["end"], entry["strand"],
                         entry.get("trunc_5p", 0), entry.get("trunc_3p", 0),
                         entry.get("allele_seq_length"))
    feature = ExtendedFeature(location=location, type=feature_key,
                              id="MGE_{0}".format(index), seq_id=seq_id)
    feature.qualifiers = build_qualifiers(entry)
    if is_putative:
        # Internal marker (not written to output). Excluded from the DDBJ ann by default
        # (they often overlap/duplicate and are unsuitable for submission), but kept in gbk/gff.
        feature.annotations["mge_putative_composite"] = True
    return seq_id, feature


def parse_mge_results(data):
    """Convert the whole MEF JSON to {seq_id: [ExtendedFeature, ...]}."""
    result = {}
    for index, entry in enumerate(data.get("result", []), start=1):
        seq_id, feature = entry_to_feature(entry, index)
        result.setdefault(seq_id, []).append(feature)
    return result


class MobileElementFinder(ContigAnnotationTool):
    """MobileElementFinder

    Tool type: contig annotation (feature)
    URL: https://bitbucket.org/mhkj/mge_finder
    REF: Johansson et al., 2021
    """
    version = None
    TYPE = "feature"
    NAME = "MobileElementFinder"
    # mefinder prints a pkg_resources DeprecationWarning to stderr.
    # base_tools.setVersion() treats any stderr output as "not found", so stderr is discarded
    # only for the version check (no effect on actual runs, which are judged by the return code).
    VERSION_CHECK_CMD = ["mefinder", "--version", "2>/dev/null"]
    VERSION_PATTERN = r"(\d+\.\d+\.\d+)"

    def __init__(self, options=None, workDir="OUT"):
        if options is None:
            options = {}
        super(MobileElementFinder, self).__init__(options, workDir)
        self.cmd_options = options.get("cmd_options", "")
        self.db_path = options.get("db_path", "")
        self.output_directory = os.path.join(self.workDir, "contig_annotation", "mge_finder")
        self.result_file = os.path.join(self.output_directory, "mge.json")
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

    def getCommand(self):
        prefix = os.path.join(self.output_directory, "mge")
        cmd = ["mefinder", "find",
               "-c", self.genomeFasta,
               "--json",
               "-t", str(self.options.get("cpu", 1)),
               "--temp-dir", self.output_directory]
        if self.db_path:
            cmd += ["--db-path", self.db_path]
        if self.cmd_options:
            cmd += self.cmd_options.split()
        cmd += [prefix]
        return cmd

    def getFeatures(self):
        """Return {seq_id: [ExtendedFeature, ...]}."""
        with open(self.result_file) as fh:
            data = json.load(fh)
        return parse_mge_results(data)

    def getResult(self):
        """ContigAnnotationTool contract: (source_notes, report).

        source_notes is empty because MEF returns located features via getFeatures().
        report holds the MGE summary for the ## lines of amr_summary.tsv (same format as PlasmidFinder).
        Types are converted to human-readable full names (_TYPE_NAME).
        """
        source_notes = {}
        report = {}
        with open(self.result_file) as fh:
            data = json.load(fh)
        for entry in data.get("result", []):
            seq_id = entry["contig"].split()[0]
            type_name = _TYPE_NAME.get(entry["type"], entry["type"])
            summary = "MobileElementFinder: {name} ({type}), identity:{id:.1f}%, coverage:{cov:.1f}%, {s}..{e}".format(
                name=entry["name"], type=type_name,
                id=float(entry["identity"]) * 100, cov=float(entry["coverage"]) * 100,
                s=entry["start"], e=entry["end"])
            report.setdefault(seq_id, []).append(summary)
        return source_notes, report
