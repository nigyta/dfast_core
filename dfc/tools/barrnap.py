#! /usr/bin/env python
# coding: UTF8

from .base_tools import StructuralAnnotationTool
from Bio import SeqIO
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

        def _getLengthDict(fileName):
            R = list(SeqIO.parse(open(fileName), "fasta"))
            return {r.id: len(r) for r in R}

        def _checkPartial(left, right, strand, seqLength):
            left_flag, right_flag = "0", "0"
            if int(left) <= 10:
                left_flag = "1"
                left = 1
            if seqLength - int(right) <= 10:
                right_flag = "1"
                right = seqLength
            partial_flag = left_flag + right_flag
            return left, right, partial_flag

        dict_length = _getLengthDict(self.genomeFasta)

        D = {}
        i = 0
        for sequence, toolName, featureType, left, right, strand, rRNA_type, product, note, partial in _parseResult():
            left, right, partial_flag = _checkPartial(left, right, strand, dict_length[sequence])
            location = self.getLocation(left, right, strand, partial_flag)
            i += 1

            annotations = {"partial_flag": partial_flag, "rRNA_type": rRNA_type}
            if partial or partial_flag != "00":
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

