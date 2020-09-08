#! /usr/bin/env python
# coding: UTF8

from .base_tools import StructuralAnnotationTool
from ..models.bio_feature import ExtendedFeature
from Bio import SeqIO


class tRNAscan(StructuralAnnotationTool):
    """
    tRNAscan class extended from StructuralAnnotationTool class

    Tool type: tRNA
    URL: 
    REF: 

    """
    version = None
    TYPE = "tRNA"
    NAME = "tRNAscan-SE"
    VERSION_CHECK_CMD = ["tRNAscan-SE", "-h", "2>&1"]
    VERSION_PATTERN = r"tRNAscan-SE (.+) \("
    SHELL = True
    def __init__(self, options=None, workDir="OUT"):
        """
        tRNAscan options:
          -B  --bact            : search for bacterial tRNAs (use bacterial tRNA model)
          -A  --arch            : search for archaeal tRNAs (use archaeal tRNA model)
          -O  --organ           : search for organellar (mitochondrial/chloroplast) tRNAs
          -G  --general         : use general tRNA model (cytoplasmic tRNAs from all 3 domains included)

          -D  --nopseudo        : disable pseudogene checking
          -o  --output <file>   : save final results in <file>
          -b  --brief           : brief output format (no column headers)
          -Q   --forceow  : do not prompt user before overwriting pre-existing result files  (for batch processing)
        :param options:
        :param workDir:
        """
        super(tRNAscan, self).__init__(options, workDir)
        self.model = options.get("model", "--bact")
        self.cmd_options = options.get("cmd_options", "")

    def getCommand(self):
        """
        tRNAscan-SE -B -b -o tRNAscanSE.tsv genome.fna 2> tRNAscanSE.log
        """
        # TODO: check pseudogene output
        if self.__class__.version[0] == "2":
            default_option = "--brief --forceow --thread 1"
        else:
            default_option = "--brief --forceow"
        cmd = ["tRNAscan-SE", self.model, self.cmd_options, default_option, "--output", self.outputFile,
               self.genomeFasta, "2>&1", self.logFile]
        return cmd

    def getFeatures(self):
        def get_length_dict():
            dict_length = {r.id: len(r) for r in SeqIO.parse(self.genomeFasta, "fasta")}
            return dict_length

        def _parseResult():
            with open(self.outputFile) as f:
                for line in f:
                    cols = line.strip("\n\t ").split("\t")
                    sequence, _, start, end, aa, anticodon, _, _, score = cols[0:9]
                    pseudo = cols[9] if len(cols) == 10 else ""

                    # Remove pseudo
                    if aa == "Undet" or aa == "Pseudo" or pseudo == "Pseudo":
                        self.logger.info("'Pseudo or Undet' tRNA found. Skipping... {}:{}..{}({})".format(
                            sequence, start, end, aa))
                        continue

                    # clean irregular
                    if aa == "Ile2":
                        aa = "Ile"
                    elif aa == "SeC(p)":
                        aa = "SeC"
                    elif aa == "(Ile)":
                        aa = "Ile"

                    if int(start) <= int(end):
                        left, right, strand = start, end, "+"
                    else:
                        left, right, strand = end, start, "-"
                    sequence = sequence.strip()  # Remove trailing blank
                    yield sequence, left, right, strand, aa, anticodon, score

        dict_length = get_length_dict()
        D = {}
        i = 0
        for sequence, left, right, strand, aa, anticodon, score in _parseResult():
            seq_length = dict_length[sequence]
            if int(left) < 1:
                left = 1
                location = self.getLocation(left, right, strand, partial_flag="10")
            elif int(right) > seq_length:
                right = seq_length
                location = self.getLocation(left, right, strand, partial_flag="01")
            else:
                location = self.getLocation(left, right, strand)


            i += 1

            type_ = "tRNA"
            annotations = {"anticodon": anticodon}
            feature = ExtendedFeature(location=location, type=type_, id="{0}_{1}".format(self.__class__.__name__, i),
                                      seq_id=sequence, annotations=annotations)
            feature.qualifiers = {
                "product": ["tRNA-{0}".format(aa)],
                "inference": ["COORDINATES:ab initio prediction:{0}:{1}".format(self.__class__.NAME, self.__class__.version)],
            }
            D.setdefault(sequence, []).append(feature)
        return D


if __name__ == '__main__':
    pass
