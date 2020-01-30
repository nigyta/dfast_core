#! /usr/bin/env python
# coding: UTF8

import re
from Bio import SeqIO
from .base_tools import StructuralAnnotationTool
from ..models.bio_feature import ExtendedFeature


class RNAmmer(StructuralAnnotationTool):
    """
    Barrnap

    Tool type: rRNA prediction
    URL: https://github.com/tseemann/barrnap
    REF:

    """
    version = None
    TYPE = "rRNA"
    NAME = "RNAmmer"
    VERSION_CHECK_CMD = ["rnammer", "-v"]
    VERSION_PATTERN = r"This rnammer (.+), running from"
    SHELL = True

    def __init__(self, options=None, workDir="OUT"):
        """
        """

        super(RNAmmer, self).__init__(options, workDir)
        self.model = options.get("model", "bac")
        self.cmd_options = options.get("cmd_options", "")

    def getCommand(self):
        """rnammer -S arc/bac/euk (-multi) (-m tsu,lsu,ssu) (-f) (-k) (-gff [gff file]) (-xml [xml file]) (-f [fasta file]) (-h [HMM report]) [sequence]"""
        cmd = ["rnammer", "-S", self.model, "-m tsu,lsu,ssu", self.cmd_options, "-gff", self.outputFile, self.genomeFasta]
        return cmd

    def getFeatures(self):
        """RNAmmer generates standard GFF format."""
        pat_product = re.compile(r"(\d+)s_rRNA")
        def _parseResult():
            with open(self.outputFile) as f:
                for line in f:
                    if line.startswith("#"):
                        continue

                    sequence, toolName, featureType, left, right, score, strand, _, qualifiers = line.strip("\n\t").split("\t")
                    """ex) 5s_rRNA"""

                    m = pat_product.search(qualifiers)
                    if m:
                        rRNAtype = m.group(1)
                        product = rRNAtype + "S ribosomal RNA"
                    else:
                        self.logger.warning("[RNAmmer] Unknown feature name for rRNA. '{0}'".format(qualifiers))
                        continue

                    yield sequence, toolName, featureType, left, right, strand, product

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
        for sequence, toolName, featureType, left, right, strand, product in _parseResult():
            left, right, partial_flag = _checkPartial(left, right, strand, dict_length[sequence])
            location = self.getLocation(left, right, strand, partial_flag)
            i += 1

            annotations = {"partial_flag": partial_flag}
            if partial_flag != "00":
                annotations["partial"] = True

            feature = ExtendedFeature(location=location, type="rRNA", id="{0}_{1}".format(self.__class__.__name__, i),
                                      seq_id=sequence, annotations=annotations)
            feature.qualifiers = {
                "product": [product],
                "inference": ["COORDINATES:profile:{0}:{1}".format(self.__class__.NAME, self.__class__.version)],
            }


            D.setdefault(sequence, []).append(feature)
        return D


if __name__ == '__main__':
    pass