#! /usr/bin/env python
# coding: UTF8

import subprocess
import re
import os.path
from logging import getLogger, StreamHandler, FileHandler, DEBUG, WARN, INFO, Formatter
from Bio.SeqFeature import FeatureLocation, ExactPosition, BeforePosition, AfterPosition


class Tool(object):
    """
    Abstract Base Tool class.
    """
    version = None
    NAME = "AbstractBaseTool"  # Should be overridden in child classes.
    VERSION_CHECK_CMD = ["echo",
                         "version 1.0.0"]  # Command to get the tool version. Should be overridden in child classes.
    VERSION_PATTERN = r"version (.+)"  # Regex pattern that matches the tool version. Should be overridden in child classes.
    VERSION_ERROR_MSG = None

    def __init__(self, options=None):
        if not options:
            options = {}
        self.options = options
        self.logger = getLogger(__name__)
        if self.__class__.version is None:
            self.setVersion()
        if self.options:
            self.logger.info("Setting {0} options. {1}".format(self.__class__.NAME, str(self.options) ))


    def setVersion(self):
        """
        Check the tool version.
        This method uses VERSION_CHECK_CMD and VERSION_PATTERN.
        """

        self.logger.info("Checking {0} version... ".format(self.__class__.NAME))
        out, err = self.executeCommand(self.__class__.VERSION_CHECK_CMD, shell=True)
        if err:
            self.logger.error("{self.NAME} not found in PATH. Aborted...".format(self=self))
            exit(1)
        else:
            m = re.search(self.__class__.VERSION_PATTERN, out.decode('utf-8'))
            if m:
                version = m.group(1)
                self.__class__.version = version
                self.logger.info("{self.__class__.NAME} initialized. (Version {self.version})".format(self=self))
            else:
                self.logger.error("{self.__class__.NAME} version could not be determined. Aborted...".format(self=self))
                if self.__class__.VERSION_ERROR_MSG:
                    self.logger.error(self.__class__.VERSION_ERROR_MSG)
                exit(1)

    def getVersion(self):
        return "{tool.__class__.NAME} {tool.version}".format(tool=self)

    def executeCommand(self, cmd, shell=False):
        """
        Any output to standard error will be handled as an error.
        Some tools write logs to standard error, which must be redirected to standard output.
        :param cmd:
        :param shell:
        :return:
        """
        self.logger.debug('Running command "{0}" ({1})'.format(" ".join(cmd), self.__class__.__name__))
        if shell:
            cmd = " ".join(cmd)
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, shell=shell)
        out, err = p.communicate()
        # self.logger.debug(out)
        if err:
            self.logger.error("Process aborted due to an error in {self.__class__.__name__}.".format(self=self))
            self.logger.error(err.decode("utf8"))
            exit(1)
        # self.logger.debug(err)
        return out, err


class StructuralAnnotationTool(Tool):
    NAME = "StructuralAnnotationTool"  # Should be overridden in child classes.
    TYPE = ""  # Must be one of the features. e.g. CDS, tRNA, rRNA
    VERSION_CHECK_CMD = ["echo",
                         "version 1.0.0"]  # Command to get the tool version. Should be overridden in child classes.
    VERSION_PATTERN = r"version (.+)"  # Regex pattern that matches the tool version. Should be overridden in child classes.
    SHELL = False

    def __init__(self, options=None, workDir="OUT"):
        super(StructuralAnnotationTool, self).__init__(options)
        self.workDir = workDir
        self.genomeFasta = os.path.join(self.workDir, "input", "genome.fna")

        self.outputFile = os.path.join(self.workDir, "structural", "{0}.txt".format(self.__class__.__name__))
        self.logFile = os.path.join(self.workDir, "structural", "{0}.log".format(self.__class__.__name__))

    def getFeatures(self):
        """Override this method"""
        raise NotImplementedError

    def getCommand(self):
        """Override this method"""
        raise NotImplementedError

    def run(self):
        cmd = self.getCommand()
        self.executeCommand(cmd, shell=self.__class__.SHELL)

    def getLocation(self, left, right, strand, partial_flag="00"):
        """partialFlag = {00:both ends existing, 10:left-end missing, 01:right-end missing, 00:both-ends missing}"""
        strand = 1 if (strand == "+" or strand == "1" or strand == 1) else -1
        leftPosition = BeforePosition(int(left) - 1) if partial_flag[0] == "1" else ExactPosition(int(left) - 1)
        rightPosition = AfterPosition(int(right)) if partial_flag[1] == "1" else ExactPosition(int(right))
        return FeatureLocation(leftPosition, rightPosition, strand=strand)


class Aligner(Tool):
    NAME = "AbstractAlignerClass"  # Should be overridden in child classes.
    VERSION_CHECK_CMD = ["echo",
                         "version 1.0.0"]  # Command to get the tool version. Should be overridden in child classes.
    VERSION_PATTERN = r"version (.+)"  # Regex pattern that matches the tool version. Should be overridden in child classes.

    def format_db_command(self, db_fasta_file, db_name):
        raise NotImplementedError

    def format_db(self, db_fasta_file, db_name):
        cmd = self.format_db_command(db_fasta_file, db_name)
        self.executeCommand(cmd, shell=False)


class JavaWrapper(Tool):
    NAME = "Java"  # Should be overridden in child classes.
    VERSION_PATTERN = r'version "(.+)"'
    
    SHELL = False

    def __init__(self, options=None, java_options=""):

        self.java_options = java_options

        super(JavaWrapper, self).__init__(options)

    def setVersion(self):
        """
        Check the tool version.
        This method uses VERSION_CHECK_CMD and VERSION_PATTERN.
        """
        version_check_command_ = ["java", self.java_options, "-version", "2>&1 /dev/null"]

        self.logger.info("Checking {0} version... ".format(self.__class__.NAME))
        out, err = self.executeCommand(version_check_command_, shell=True)
        m = re.search(self.__class__.VERSION_PATTERN, out.decode('utf-8'))
        if m:
            version = m.group(1)
            self.__class__.version = version
            self.logger.info("Java initialized. (Version {self.version})".format(self=self))
        else:
            self.logger.error("Java version could not be determined. Aborted...\n" +
                              "Please check if Java is installed in your system.\n" +
                              "Depending on your system, Java requires additional options to run, such as '-Xmx512m'. " +
                              "You can set them to java_options in the config file or the environmental variable _JAVA_OPTIONS.")
            exit()

    def getCommand(self):
        return ["java", self.java_options]

if __name__ == '__main__':
    logger = getLogger(__name__)

    logger.setLevel(DEBUG)

    handler = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)

    tool = JavaWrapper()
    print(tool.version)

