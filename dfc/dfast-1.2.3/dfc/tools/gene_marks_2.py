#! /usr/bin/env python
# coding: UTF8

from .base_tools import StructuralAnnotationTool
from Bio import SeqIO
# from Bio.SeqFeature import SeqFeature, FeatureLocation
from ..models.bio_feature import ExtendedFeature


class GeneMarkS2(StructuralAnnotationTool):
    """
    GeneMarkS2

    Tool type: CDS prediction
    URL: 
    REF:

    """
    version = None
    TYPE = "CDS"
    NAME = "GeneMarkS2"
    VERSION_CHECK_CMD = ["/home/apfennig3/Team1-GenePrediction/bin/gms2.pl | tail -n 1", "1>&14"]
    VERSION_PATTERN = r"Version: (.+).lic"
    SHELL = True

    def __init__(self, options=None, workDir="OUT"):
        """
        """

        super(GeneMarkS2, self).__init__(options, workDir)
        self.transl_table = options.get("transl_table", 11)
        self.cmd_options = options.get("cmd_options", "")

    def getCommand(self):
        """/home/apfennig3/Team1-GenePrediction/bin/gms2.pl --seq $genome --genome-type bacteria --output $tmp_dir/output.gtf --format gff"""
        cmd = ["/home/apfennig3/Team1-GenePrediction/bin/gms2.pl", self.cmd_options, "--genome-type", "bacteria", "--gcode", str(self.transl_table), "--format", "gff", "--seq", self.genomeFasta, "--ouput", self.outputFile]
        return cmd


    def getFeatures(self):
        """GeneMarkS2 generates standard GFF format.
           

        """

        def _parseResult():
            with open(self.outputFile) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    sequence, toolName, featureType, left, right, score, strand, _, qualifiers = line.strip("\n").split("\t")
                    qualifiers = dict([x.split("=") for x in qualifiers.strip(";").split(";")])
                    yield sequence, toolName, featureType, left, right, strand, qualifiers

        def _getLengthDict(fileName):
            R = list(SeqIO.parse(open(fileName), "fasta"))
            return {r.id: len(r) for r in R}

        def _get_feature(left, right, strand, partial_flag, seq_length, i):
            left_flag, right_flag = partial_flag[0], partial_flag[1]
            left, right, codon_start = int(left), int(right), 1

            if left_flag == "1" and left <= 3:
                if strand == "+":
                    codon_start = left
                left = 1
            if right_flag == "1" and seq_length - right <= 2:
                if strand == "-":
                    codon_start = seq_length - right + 1
                right = seq_length
            location = self.getLocation(left, right, strand, partial_flag)

            annotations = {"partial_flag": partial_flag}
            if partial_flag != "00":
                annotations["partial"] = True
            feature = ExtendedFeature(location=location, type="CDS", id="{0}_{1}".format(self.__class__.__name__, i),
                                      seq_id=sequence, annotations=annotations)
            feature.qualifiers = {
                "product": ["hypothetical protein"],
                "inference": ["COORDINATES:ab initio prediction:{0}:{1}".format(self.NAME, self.version)],
                "transl_table": [str(self.transl_table)],
                "codon_start": [codon_start]
            }
            return feature

        dict_length = _getLengthDict(self.genomeFasta)

        D = {}
        i = 0
        for sequence, toolName, featureType, left, right, strand, qualifiers in _parseResult():
            partial_flag = qualifiers.get("partial", "00")
            seq_length = dict_length[sequence]

            i += 1
            feature = _get_feature(left, right, strand, partial_flag, seq_length, i)

            if qualifiers.get("rbs_motif", "None") != "None":
                feature.annotations["rbs"] = qualifiers["rbs_motif"]
            D.setdefault(sequence, []).append(feature)
        return D


