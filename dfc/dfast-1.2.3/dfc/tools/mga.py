#! /usr/bin/env python
# coding: UTF8

from .base_tools import StructuralAnnotationTool
# from Bio.SeqFeature import SeqFeature, FeatureLocation
from ..models.bio_feature import ExtendedFeature


class MGA(StructuralAnnotationTool):
    """
    MetaGeneAnnotator

    Tool type: CDS prediction
    URL: http://metagene.nig.ac.jp/metagene/download_mga.html
    REF: https://academic.oup.com/dnaresearch/article/15/6/387/512877/MetaGeneAnnotator-Detecting-Species-Specific

    """
    version = None
    TYPE = "CDS"
    NAME = "MetaGeneAnnotator"
    VERSION_CHECK_CMD = None
    VERSION_PATTERN = None

    SHELL = True  # Run through shell

    def __init__(self, options=None, workDir="OUT"):
        super(MGA, self).__init__(options, workDir)
        self.transl_table = options.get("transl_table", 11)  # MGA always uses transl_table = 11
        self.cmd_options = options.get("cmd_options", "-s")  # -s for single species, -m for multiple species

    def setVersion(self):
        """
        MGA does not have an option to check the version.
        The bundled version of MGA is released 2008/08/19.
        """
        version = "2008/08/19"
        self.logger.info("Checking {0}... ".format(self.__class__.NAME))
        out, err = self.executeCommand(["mga 2>&1"], shell=True)
        if out.decode("utf-8").startswith("usage: mga"):
            self.__class__.version = version
            self.logger.info("{self.__class__.NAME} initialized. (Version {self.__class__.version})".format(self=self))
        else:
            self.logger.error("{self.__class__.NAME} not found in PATH. Aborted...".format(self=self))
            exit()

    def getCommand(self):
        """mga -s testDir/INPUT/genome.fna"""
        cmd = ["mga", self.cmd_options, self.genomeFasta, ">", self.outputFile]
        return cmd

    def getFeatures(self):

        def _parseResult():
            with open(self.outputFile) as f:
                counter = 0
                for line in f:
                    if counter == 0 and line.startswith("#"):
                        counter += 1
                        seq_name = line[2:].strip()  # sequence1
                    elif counter == 1 and line.startswith("#"):
                        counter += 1
                    elif counter == 2 and line.startswith("#"):
                        counter = 0
                    else:
                        yield [seq_name] + line.strip("\n").split("\t")

        D = {}
        i = 0
        for sequence, _, left, right, strand, frame, complete_flag, gene_score, \
            model, rbs_left, rbs_right, rbs_score in _parseResult():

            """MGA completeflag= {11:both ends exist, 10:stop codon missing, 01:start codon missing, 00:both missing}
               To make partial flag, number should be inversed (e.g. 0 to 1 and 1 to 0),
               and the order should be reversed if strand=-1
            """
            partial_flag = complete_flag.replace("1", "9").replace("0", "1").replace("9", "0")
            if strand == "-":
                partial_flag = partial_flag[::-1]
            location = self.getLocation(left, right, strand, partial_flag)

            annotations = {"partial_flag": partial_flag}
            if partial_flag != "00":
                annotations["partial"] = True

            i += 1
            feature = ExtendedFeature(location=location, type="CDS", id="{0}_{1}".format(self.__class__.__name__, i),
                                      seq_id=sequence, annotations=annotations)
            feature.qualifiers = {
                "product": ["hypothetical protein"],
                # "translation": [""],
                "inference": ["COORDINATES:ab initio prediction:{0}".format(self.NAME)],
                "transl_table": [str(self.transl_table)],
                "codon_start": [int(frame) + 1]
            }

            if partial_flag != "00":
                feature.annotations["partial"] = True
            if rbs_score != "." and rbs_score != "-":
                rbs_location = self.getLocation(rbs_left, rbs_right, strand, partial_flag="00")
                rbs = ExtendedFeature(location=rbs_location, type="regulatory",
                                      id="{0}-RBS_{1}".format(self.__class__.__name__, i).format(i), seq_id=sequence)
                rbs.qualifiers = {"regulatory_class": ["ribosome_binding_site"],
                                  "inference": ["COORDINATES:ab initio prediction:{0}".format(self.NAME)],
                                  }
                feature.annotations["rbs"] = rbs
            D.setdefault(sequence, []).append(feature)
        return D


if __name__ == '__main__':
    from logging import getLogger, INFO, DEBUG, StreamHandler

    logger = getLogger("app")

    logger.setLevel(DEBUG)

    handler = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)

    tool = MGA(workDir="/Users/tanizawa/projects/newpipe/dfast_core/OUT")
    tool.run()
    print(tool.getFeatures())
