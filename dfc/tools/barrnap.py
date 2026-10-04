#! /usr/bin/env python
# coding: UTF8

import argparse
import os
import shlex

from .base_tools import StructuralAnnotationTool
from ..models.bio_feature import ExtendedFeature

# rRNA profile HMMs, downloaded into DB_ROOT/rrna by `dfast_file_downloader.py --rrna`, one file per kingdom.
# "barrnap" is db/<kingdom>.hmm of Barrnap 0.9 (same as 0.8). "rfam" is built from the Rfam seed alignments
# of RFAM_FAMILIES with hmmbuild, as pybarrnap does.
KINGDOMS = ["bac", "arc", "euk"]  # also the order that wins a tie in compete()
RFAM_RELEASE = "15.1"
RFAM_FAMILIES = {
    "bac": [("16S_rRNA", "RF00177"), ("23S_rRNA", "RF02541"), ("5S_rRNA", "RF00001")],
    "arc": [("16S_rRNA", "RF01959"), ("23S_rRNA", "RF02540"), ("5S_rRNA", "RF00001")],
    "euk": [("18S_rRNA", "RF01960"), ("28S_rRNA", "RF02543"), ("5S_rRNA", "RF00001")],
}
# model name: (file name pattern for a kingdom, label used in the inference qualifier)
RRNA_MODELS = {
    "barrnap": ("barrnap_0.9_{}.hmm", "BarrnapHmm:0.9"),
    "rfam": ("Rfam_" + RFAM_RELEASE + "_{}.hmm", "Rfam:" + RFAM_RELEASE),
}
# Expected rRNA lengths used by Barrnap to reject short hits and tag partial ones.
# Other models (5.8S in the Barrnap arc/euk files) are ignored.
EXPECTED_LENGTH = {"5S_rRNA": 119, "16S_rRNA": 1585, "23S_rRNA": 3232, "18S_rRNA": 1869, "28S_rRNA": 2912}
SUBUNIT = {"5S_rRNA": "5S", "16S_rRNA": "SSU", "18S_rRNA": "SSU", "23S_rRNA": "LSU", "28S_rRNA": "LSU"}
EUKARYOTIC_RRNA = {"18S_rRNA", "28S_rRNA"}  # reported as possible contamination
MAX_LENGTH = 3878  # Barrnap's nhmmer --w_length: int(1.2 * 3232)


def read_nhmmer_tblout(tblout, kingdom):
    """Hits in nhmmer --tblout as dicts."""
    hits = []
    with open(tblout) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            x = line.split()
            if x[2] not in EXPECTED_LENGTH:
                continue
            ali_from, ali_to = int(x[6]), int(x[7])
            left, right, strand = (ali_from, ali_to, "+") if ali_from < ali_to else (ali_to, ali_from, "-")
            hits.append({"seq_id": x[0], "gene": x[2], "left": left, "right": right, "strand": strand,
                         "evalue": x[12], "score": float(x[13]), "kingdom": kingdom})
    return hits


def compete(hits):
    """
    Among overlapping hits of the same subunit on the same strand from different kingdoms, keep the one
    with the best bit score (as in the clan competition of Rfam). Ties go to the earlier kingdom in KINGDOMS.
    Hits of the same kingdom are left as they are, so a single kingdom gives the same result as Barrnap.
    """
    kept = []
    for h in sorted(hits, key=lambda h: (-h["score"], KINGDOMS.index(h["kingdom"]))):
        if not any(k["kingdom"] != h["kingdom"] and k["seq_id"] == h["seq_id"] and k["strand"] == h["strand"]
                   and SUBUNIT[k["gene"]] == SUBUNIT[h["gene"]] and k["left"] <= h["right"] and h["left"] <= k["right"]
                   for k in kept):
            kept.append(h)
    return kept


def apply_barrnap_rules(hits, reject=0.5, lencutoff=0.8):
    """
    Barrnap 0.8 rules: reject hits shorter than reject * expected length, and tag hits shorter than
    lencutoff * expected length as partial ("aligned only N percent ...").
    Returns (seq_id, left, right, strand, gene, product, note, evalue), sorted by seq_id and left like Barrnap.
    """
    results = []
    for h in hits:
        gene, expected, length = h["gene"], EXPECTED_LENGTH[h["gene"]], h["right"] - h["left"] + 1
        if length < int(reject * expected):
            continue
        product = gene.replace("_r", " ribosomal ", 1)  # 16S_rRNA -> 16S ribosomal RNA
        note = None
        if length < int(lencutoff * expected):
            note = "aligned only {0} percent of the {1}".format(int(100 * length / expected), product)
        results.append((h["seq_id"], h["left"], h["right"], h["strand"], gene, product, note, h["evalue"]))
    return sorted(results, key=lambda r: (r[0], r[1]))


class Barrnap(StructuralAnnotationTool):
    """
    rRNA prediction with the Barrnap method: nhmmer (HMMER) is run directly against rRNA profile HMMs.
    This is a Python reimplementation of Barrnap by Torsten Seemann (GPL-3.0), written with reference to
    the source code of Barrnap 0.8, and the default models are those of Barrnap 0.9.

    Tool type: rRNA prediction
    URL: https://github.com/tseemann/barrnap
    REF:

    """
    version = None
    TYPE = "rRNA"
    NAME = "nhmmer"
    VERSION_CHECK_CMD = ["nhmmer", "-h"]
    VERSION_PATTERN = r"# HMMER (\S+)"

    def __init__(self, options=None, workDir="OUT"):
        """
        options:
            model: "barrnap" (default) or "rfam", see RRNA_MODELS
            kingdoms: prokaryotic models to search, "bac" and/or "arc" (default: ["bac"]).
            check_eukaryotic: also search the euk models as a contamination check (default: True).
                Overlapping hits of different kingdoms compete by bit score; 18S/28S hits that win are
                reported as possible contamination.
            db_dir: directory holding the model files (DB_ROOT/rrna)
            cmd_options: Barrnap options "--reject", "--lencutoff" and "--evalue"
        """
        options = options or {}
        super(Barrnap, self).__init__(options, workDir)
        model = options.get("model", "barrnap")
        file_pattern, self.model_label = RRNA_MODELS[model]
        kingdoms = options.get("kingdoms", ["bac"])
        unknown = set(kingdoms) - {"bac", "arc"}
        if unknown or not kingdoms:
            self.logger.error("Invalid rRNA kingdoms: {0}. Choose bac and/or arc.".format(", ".join(kingdoms)))
            exit(1)
        self.kingdoms = list(kingdoms) + (["euk"] if options.get("check_eukaryotic", True) else [])
        self.hmm_files = {k: os.path.join(options.get("db_dir", ""), file_pattern.format(k)) for k in self.kingdoms}
        parser = argparse.ArgumentParser(prog="Barrnap cmd_options", add_help=False)
        parser.add_argument("--reject", type=float, default=0.5)  # Barrnap 0.8 default
        parser.add_argument("--lencutoff", type=float, default=0.8)
        parser.add_argument("--evalue", type=float, default=1e-6)
        self.params = parser.parse_args(shlex.split(options.get("cmd_options", "")))
        for hmm_file in self.hmm_files.values():
            if not os.path.exists(hmm_file):
                self.logger.error("rRNA model file not found: {0}. Run 'dfast_file_downloader.py --rrna {1}'.".format(
                    hmm_file, model))
                exit(1)

    def tblout(self, kingdom):
        return self.outputFile.replace(".txt", ".{0}.txt".format(kingdom))

    def getCommand(self, kingdom="bac"):
        return ["nhmmer", "--cpu", "1", "-E", str(self.params.evalue), "--w_length", str(MAX_LENGTH),
                "-o", os.devnull, "--tblout", self.tblout(kingdom), self.hmm_files[kingdom], self.genomeFasta]

    def run(self):
        for kingdom in self.kingdoms:
            self.executeCommand(self.getCommand(kingdom))

    def getFeatures(self):
        D = {}
        hits = compete([h for k in self.kingdoms for h in read_nhmmer_tblout(self.tblout(k), k)])
        for i, (sequence, left, right, strand, rRNA_type, product, note, _) in enumerate(
                apply_barrnap_rules(hits, self.params.reject, self.params.lencutoff), 1):
            # Coordinates are kept as reported. Contig ends and assembly gaps are handled by
            # FeatureUtil.adjust_rrna_features(), which also treats partial hits ("aligned only N percent")
            # and turns eukaryotic rRNAs into misc_feature.
            location = self.getLocation(left, right, strand)
            annotations = {"partial_flag": "00", "rRNA_type": rRNA_type}
            if note:
                annotations["partial"] = True
            if rRNA_type in EUKARYOTIC_RRNA:
                annotations["eukaryotic_rrna"] = True
            feature = ExtendedFeature(location=location, type="rRNA", id="{0}_{1}".format(self.__class__.__name__, i),
                                      seq_id=sequence, annotations=annotations)
            feature.qualifiers = {
                "product": [product],
                "inference": ["COORDINATES:profile:{0}".format(self.model_label)],
            }
            if note:
                feature.qualifiers["note"] = [note]
            D.setdefault(sequence, []).append(feature)
        return D
