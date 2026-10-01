#! /usr/bin/env python
# coding: UTF8

from .base_tools import StructuralAnnotationTool
# from Bio.SeqFeature import SeqFeature, FeatureLocation
from ..models.bio_feature import ExtendedFeature


class Barrnap(StructuralAnnotationTool):
    """
    Barrnap

    Tool type: rRNA prediction
    URL: https://github.com/tseemann/barrnap
    REF:

    """
    version = None
    TYPE = "rRNA"
    NAME = "Barrnap"
    VERSION_CHECK_CMD = ["barrnap", "--version", "2>&1"]
    VERSION_PATTERN = r"barrnap (.+)$"
    VERSION_ERROR_MSG = "This may happen if Time::Piece cannot be found. " + \
      "If you are using CentOS/RedHat, try 'sudo yum install perl-Time-Piece'."
    SHELL = True

    def __init__(self, options=None, workDir="OUT"):
        """
        """

        super(Barrnap, self).__init__(options, workDir)
        self.cmd_options = options.get("cmd_options", "")
        # Barrnap's own default differs by version (0.8: 0.5, 0.9: 0.25), so set it explicitly.
        if "--reject" not in self.cmd_options:
            self.cmd_options = ("--reject 0.5 " + self.cmd_options).strip()

    def getCommand(self):
        """barrnap --threads 1 genome.fna > out.gff 2> out.log"""
        cmd = ["barrnap", "--threads", "1", self.cmd_options, self.genomeFasta, ">", self.outputFile, "2>", self.logFile]
        return cmd

    def getFeatures(self):
        """Barrnap generates standard GFF format."""

        def _parseResult():
            with open(self.outputFile) as f:
                for line in f:
                    if line.startswith("#"):
                        continue

                    sequence, toolName, featureType, left, right, score, strand, _, qualifiers = line.strip("\n").split("\t")
                    """ex) Name=23S_rRNA;product=23S ribosomal RNA (partial);note=aligned only 30 percent of the 23S ribosomal RNA"""
                    qualifiers = [x.split("=")[1] for x in qualifiers.split(";")]
                    rRNA_type = qualifiers[0]
                    product = qualifiers[1].replace(" (partial)", "")
                    if len(qualifiers) == 3:  # partial
                        note = qualifiers[2]
                        partial = True
                    else:
                        note = None
                        partial = False

                    yield sequence, toolName, featureType, left, right, strand, rRNA_type, product, note, partial

        D = {}
        i = 0
        for sequence, toolName, featureType, left, right, strand, rRNA_type, product, note, partial in _parseResult():
            # Coordinates are kept as reported. Contig ends and assembly gaps are handled by
            # FeatureUtil.adjust_rrna_features(), which also treats Barrnap partial hits ("aligned only N percent").
            location = self.getLocation(left, right, strand)
            i += 1

            annotations = {"partial_flag": "00", "rRNA_type": rRNA_type}
            if partial:
                annotations["partial"] = True

            feature = ExtendedFeature(location=location, type="rRNA", id="{0}_{1}".format(self.__class__.__name__, i),
                                      seq_id=sequence, annotations=annotations)
            feature.qualifiers = {
                "product": [product],
                "inference": ["COORDINATES:profile:{0}:{1}".format(self.__class__.NAME, self.__class__.version)],
            }
            if note:
                feature.qualifiers["note"] = [note]

            D.setdefault(sequence, []).append(feature)
        return D

