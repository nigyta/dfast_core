#! /usr/bin/env python
# coding: UTF8

from .base_tools import StructuralAnnotationTool
from Bio import SeqIO
# from Bio.SeqFeature import SeqFeature, FeatureLocation
from ..models.bio_feature import ExtendedFeature


class Prodigal(StructuralAnnotationTool):
    """
    Prodigal

    Tool type: CDS prediction
    URL: 
    REF:

    """
    version = None
    TYPE = "CDS"
    NAME = "Prodigal"
    VERSION_CHECK_CMD = ["prodigal", "-v", "2>&1"]
    VERSION_PATTERN = r"Prodigal V(.+):"
    SHELL = True

    def __init__(self, options=None, workDir="OUT"):
        """
        """

        super(Prodigal, self).__init__(options, workDir)
        self.transl_table = options.get("transl_table", 11)
        self.cmd_options = options.get("cmd_options", "")

    def getCommand(self):
        """prodigal -f gff genome.fna > out.gff 2> out.log"""
        cmd = ["prodigal", self.cmd_options, "-g", str(self.transl_table), "-f", "gff", "-i", self.genomeFasta, ">", self.outputFile, "2>", self.logFile]
        return cmd


    def getFeatures(self):
        """Prodigal generates standard GFF format.
           sequence001	Prodigal_v2.6.3	CDS	166312	166752	57.7	+	0	ID=1_163;partial=00;start_type=ATG;rbs_motif=AGGAGG;rbs_spacer=5-10bp;gc_cont=0.463;conf=100.00;score=57.67;cscore=39.08;sscore=18.59;rscore=15.73;uscore=-0.22;tscore=3.72;

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


