#! /usr/bin/env python
# coding: UTF8

import subprocess

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
    # Single mode trains on the input itself; Prodigal itself recommends at least 20 kb for training.
    MIN_SINGLE_LENGTH = 20000

    def __init__(self, options=None, workDir="OUT"):
        """
        """

        super(Prodigal, self).__init__(options, workDir)
        self.transl_table = options.get("transl_table", 11)
        self.cmd_options = options.get("cmd_options", "")

        with open(self.genomeFasta) as f:
            self.total_length = sum(len(r) for r in SeqIO.parse(f, "fasta"))
        if self.total_length < self.MIN_SINGLE_LENGTH and str(self.transl_table) != "11":
            self.logger.error("The total input length ({0:,} bp) is below {1:,} bp, where Prodigal must use its "
                              "pre-trained metagenomic models, which support only genetic code 11 (transl_table is {2}). "
                              "Provide a longer input or use genetic code 11. Aborted...".format(
                                  self.total_length, self.MIN_SINGLE_LENGTH, self.transl_table))
            exit(1)

    def getCommand(self):
        """prodigal -f gff genome.fna > out.gff 2> out.log
        Inputs shorter than MIN_SINGLE_LENGTH use meta mode (-p meta) unless -p is given in cmd_options."""
        cmd = ["prodigal", self.cmd_options]
        if self.total_length < self.MIN_SINGLE_LENGTH and "-p" not in self.cmd_options.split():
            cmd += ["-p", "meta"]
            self.logger.info("Input is shorter than {0:,} bp ({1:,} bp); Prodigal will run with '-p meta'.".format(
                self.MIN_SINGLE_LENGTH, self.total_length))
        cmd += ["-g", str(self.transl_table), "-f", "gff", "-i", self.genomeFasta, ">", self.outputFile, "2>", self.logFile]
        return cmd

    def run(self):
        """Prodigal's stderr goes to the log file, so the generic stderr check cannot see a failure; check the exit code."""
        cmd = " ".join(self.getCommand())
        self.logger.debug('Running command "{0}" ({1})'.format(cmd, self.__class__.__name__))
        if subprocess.run(cmd, shell=True).returncode != 0:
            self.logger.error("Process aborted due to an error in {0}.".format(self.__class__.__name__))
            with open(self.logFile) as f:
                self.logger.error(f.read())
            exit(1)

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
            if "partial" not in qualifiers:
                self.logger.error("The Prodigal output has no 'partial' attribute, which is needed to detect genes "
                                  "cut at sequence ends. This Prodigal build is broken (e.g. the Ubuntu apt package); "
                                  "install Prodigal from Bioconda. Aborted...")
                exit(1)
            partial_flag = qualifiers["partial"]
            seq_length = dict_length[sequence]

            i += 1
            feature = _get_feature(left, right, strand, partial_flag, seq_length, i)

            if qualifiers.get("rbs_motif", "None") != "None":
                feature.annotations["rbs"] = qualifiers["rbs_motif"]
            D.setdefault(sequence, []).append(feature)
        return D


