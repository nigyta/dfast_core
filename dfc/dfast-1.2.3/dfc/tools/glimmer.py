#! /usr/bin/env python
# coding: UTF8

from .base_tools import StructuralAnnotationTool
from Bio import SeqIO
# from Bio.SeqFeature import SeqFeature, FeatureLocation
from ..models.bio_feature import ExtendedFeature


class Glimmer(StructuralAnnotationTool):
    """
    Barrnap

    Tool type: CDS prediction
    URL: 
    REF:

    """
    version = None
    TYPE = "CDS"
    NAME = "Glimmer"
    VERSION_CHECK_CMD = [""]
    VERSION_PATTERN = r""
    SHELL = True

    def __init__(self, options=None, workDir="OUT"):
        """
        Not implemented.
        """
        raise NotImplementedError
        super(Glimmer, self).__init__(options, workDir)
        # self.cmd_options = options.get("cmd_options", "-l")

    def setVersion(self):
        """
        Glimmer does not have an option to check the version.
        As of May 1 2017, The latest version is 3.02.
        """
        version = "3.02"
        self.logger.info("Checking {0}... ".format(self.__class__.NAME))
        out, err = self.executeCommand(["mga 2>&1"], shell=True)
        if out.decode("utf-8").startswith("usage: mga"):
            self.__class__.version = version
            self.logger.info("{self.__class__.NAME} initialized. (Version {self.__class__.version})".format(self=self))
        else:
            self.logger.error("{self.__class__.NAME} not found in PATH. Aborted...".format(self=self))
            exit()

    def getCommand(self):
        """
        """

        cmd_str = """mkdir -p {self.workDir}/structural/glimmer ; cd {self.workDir}/structural/glimmer ;
        long-orfs -n -t 1.15 ../../input/genome.fna result.longorfs 2> glimmer.log; 
        extract -t ../../input/genome.fna result.longorfs > result.train 2>> glimmer.log; 
        build-icm -r result.icm < result.train 2>> glimmer.log; 
        glimmer3 -o50 -g110 -t30 ../../input/genome.fna result.icm result 2>> glimmer.log;
        cp result.predict ../glimmer.txt 2>> glimmer.log""".format(self=self)

        cmd = [cmd_str]
        return cmd

    def getFeatures(self):
        pass