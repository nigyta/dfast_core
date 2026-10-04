#! /usr/bin/env python
# coding: UTF8

import argparse
import os
import shlex

from .base_tools import StructuralAnnotationTool
from ..models.bio_feature import ExtendedFeature

# rRNA profile HMMs, downloaded into DB_ROOT/rrna by `dfast_file_downloader.py --rrna`.
# "barrnap" is db/bac.hmm of Barrnap 0.9 (same as 0.8). "rfam" is built from the Rfam seed alignments
# of RFAM_FAMILIES with hmmbuild, as pybarrnap does.
RFAM_RELEASE = "15.1"
RFAM_FAMILIES = [("16S_rRNA", "RF00177"), ("23S_rRNA", "RF02541"), ("5S_rRNA", "RF00001")]
# model name: (file name, label used in the inference qualifier)
RRNA_MODELS = {
    "barrnap": ("barrnap_0.9_bac.hmm", "BarrnapHmm:0.9"),
    "rfam": ("Rfam_{0}_bac.hmm".format(RFAM_RELEASE), "Rfam:" + RFAM_RELEASE),
}
# Expected rRNA lengths used by Barrnap to reject short hits and tag partial ones.
EXPECTED_LENGTH = {"5S_rRNA": 119, "16S_rRNA": 1585, "23S_rRNA": 3232}
MAX_LENGTH = 3878  # Barrnap's nhmmer --w_length: int(1.2 * 3232)


def parse_nhmmer_tblout(tblout, reject=0.5, lencutoff=0.8):
    """
    Read nhmmer --tblout and apply the Barrnap 0.8 rules.
    Yields (seq_id, left, right, strand, gene, product, note, evalue), sorted by seq_id and left like Barrnap.
    """
    hits = []
    with open(tblout) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            x = line.split()
            ali_from, ali_to = int(x[6]), int(x[7])
            left, right, strand = (ali_from, ali_to, "+") if ali_from < ali_to else (ali_to, ali_from, "-")
            seq_id, gene, evalue = x[0], x[2], x[12]
            product = gene.replace("_r", " ribosomal ", 1)  # 16S_rRNA -> 16S ribosomal RNA
            expected, length = EXPECTED_LENGTH[gene], right - left + 1
            if length < int(reject * expected):
                continue
            note = None
            if length < int(lencutoff * expected):
                note = "aligned only {0} percent of the {1}".format(int(100 * length / expected), product)
            hits.append((seq_id, left, right, strand, gene, product, note, evalue))
    return sorted(hits, key=lambda h: (h[0], h[1]))


class Barrnap(StructuralAnnotationTool):
    """
    rRNA prediction with the Barrnap method: nhmmer (HMMER) is run directly against rRNA profile HMMs.

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
            db_dir: directory holding the model files (DB_ROOT/rrna)
            cmd_options: Barrnap options "--reject", "--lencutoff" and "--evalue"
        """
        options = options or {}
        super(Barrnap, self).__init__(options, workDir)
        model = options.get("model", "barrnap")
        file_name, self.model_label = RRNA_MODELS[model]
        self.hmm_file = os.path.join(options.get("db_dir", ""), file_name)
        parser = argparse.ArgumentParser(prog="Barrnap cmd_options", add_help=False)
        parser.add_argument("--reject", type=float, default=0.5)  # Barrnap 0.8 default
        parser.add_argument("--lencutoff", type=float, default=0.8)
        parser.add_argument("--evalue", type=float, default=1e-6)
        self.params = parser.parse_args(shlex.split(options.get("cmd_options", "")))
        if not os.path.exists(self.hmm_file):
            self.logger.error("rRNA model file not found: {0}. Run 'dfast_file_downloader.py --rrna {1}'.".format(
                self.hmm_file, model))
            exit(1)

    def getCommand(self):
        return ["nhmmer", "--cpu", "1", "-E", str(self.params.evalue), "--w_length", str(MAX_LENGTH),
                "-o", os.devnull, "--tblout", self.outputFile, self.hmm_file, self.genomeFasta]

    def getFeatures(self):
        D = {}
        hits = parse_nhmmer_tblout(self.outputFile, self.params.reject, self.params.lencutoff)
        for i, (sequence, left, right, strand, rRNA_type, product, note, _) in enumerate(hits, 1):
            # Coordinates are kept as reported. Contig ends and assembly gaps are handled by
            # FeatureUtil.adjust_rrna_features(), which also treats partial hits ("aligned only N percent").
            location = self.getLocation(left, right, strand)
            annotations = {"partial_flag": "00", "rRNA_type": rRNA_type}
            if note:
                annotations["partial"] = True
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
