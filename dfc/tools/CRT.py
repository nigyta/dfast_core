#! /usr/bin/env python
# coding: UTF8

import os
from logging import getLogger, DEBUG, StreamHandler
from collections import namedtuple, Counter
from Bio.SeqFeature import FeatureLocation, BeforePosition, AfterPosition
from Bio import SeqIO
from .base_tools import JavaWrapper, StructuralAnnotationTool
from ..models.bio_feature import ExtendedFeature

CrtSeq = namedtuple("CrtSeq", ["id", "length", "start", "end"])

class CRT(StructuralAnnotationTool):
    """
    CRT (CRISPR recognition tool)

    Tool type: CRISPR prediction
    URL: http://www.room220.com/crt/
    REF: Bland C, Ramsey TL, Sabree F, Lowe M, Brown K, Kyrpides NC, Hugenholtz P
         CRISPR Recognition Tool (CRT): a tool for automatic detection of clustered regularly interspaced palindromic repeats.
         BMC Bioinformatics. 2007 Jun 18;8(1):209

        usage:  java -cp CRT1.2-CLI.jar crt query.fna result.txt >& result.log 

    """

    version = None
    TYPE = "CRISPR (later will be replaced by repeat_region)"
    NAME = "CRT"
    VERSION_CHECK_CMD = ""
    VERSION_PATTERN = ""
    SHELL = True

    INTERVAL = 200
    # CRT does not accept a draft genome consisting of multi-fasta.
    # Sequences are concatenated with interval of 200 Ns before runnning CRT.
    # 20230613: Interval was changed from 100 to 200 as CRT finds spacer sequences containing 100Ns 

    def __init__(self, options=None, workDir="OUT"):
        """
        """
        if options is None:
            options = {}
        jar_file = options.get("jar_file", "")
        self.jar_file = jar_file

        super(CRT, self).__init__(options, workDir)
        self.cmd_options = options.get("cmd_options", "")

        java_options = options.get("java_options", "")
        self.java = JavaWrapper(java_options=java_options)
        self.cmd_options = options.get("cmd_options", "")
        self.concatenated_query, self.seq_info = self.prepare_concatenated_query()

    def setVersion(self):
        """
        CRT does not have a command to get the version number.
        """
        self.logger.info("Checking {0} version... ".format(self.__class__.NAME))
        if not os.path.exists(self.jar_file):
            self.logger.error("Jar file for CRT not found. Aborted. [jar_file={0}]".format(self.jar_file))
            exit(1)

        self.__class__.version = "1.2"
        self.logger.info("{self.__class__.NAME} initialized. (Version {self.__class__.version})".format(self=self))

    def getCommand(self):
        cmd = self.java.getCommand() + [
            "-cp", self.jar_file, "crt", self.cmd_options, self.concatenated_query, self.outputFile, ">", self.logFile, "2>&1"]
        return cmd

    def prepare_concatenated_query(self):

        start = 1
        seq_info = []
        seqs = []
        concatenated_query = os.path.join(self.workDir, "structural", "crt_query.fna")
        for i, record in enumerate(SeqIO.parse(open(self.genomeFasta), "fasta")):
            end = start + len(record) - 1
            crt_seq = CrtSeq(record.id, len(record), start, end)
            seq_info.append(crt_seq)
            start = end + self.__class__.INTERVAL + 1
            seqs.append(str(record.seq))
        with open(concatenated_query, "w") as f:
            seq_interval = "n" * self.__class__.INTERVAL
            seq = ">query\n" + seq_interval.join(seqs) + "\n"
            f.write(seq)
        return concatenated_query, seq_info

    def getFeatures(self):

        def _get_feature(line):
            # modified at 2017.8.15. Suppress error when a crispr is detected at the edge of a contig.
            cols = line.strip().split()
            start, end = int(cols[3]), int(cols[5])

            extracted = [x for x in self.seq_info if x.start - self.__class__.INTERVAL <= start and end <= x.end + self.__class__.INTERVAL]
            # extracted = [x for x in self.seq_info if x.start <= start and end <= x.end]
            if len(extracted) != 1:
                print(extracted)
                # print(self.seq_info)
                print(start, end)
            assert len(extracted) == 1
            crt_seq = extracted[0]
            seq_id = crt_seq.id
            start = start - crt_seq.start
            end = end - crt_seq.start + 1

            start_position = BeforePosition(0) if start <= 0 else start
            end_position = AfterPosition(crt_seq.length) if end >= crt_seq.length else end
            location = FeatureLocation(start_position, end_position, strand=1)
            return ExtendedFeature(location=location, type="CRISPR", seq_id=seq_id)

        def _get_repeat_sequence(fh):
            L = []
            next(fh)
            next(fh)
            line = next(fh)
            while not line.startswith("-"):
                L.append(line.split()[1])
                line = next(fh)
            most_common, count = Counter(L).most_common()[0]
            return most_common

        i = 0
        D = {}
        with open(self.outputFile) as fh:
            for line in fh:
                if line.startswith("ORGANISM"):
                    # seq_id = line.strip().split()[1]
                    pass
                if line.startswith("CRISPR "):
                    i += 1
                    feature = _get_feature(line)
                    feature.id = "{0}_{1}".format(self.__class__.__name__, i)
                    repeat = _get_repeat_sequence(fh)
                    feature.qualifiers["inference"] = \
                        ["COORDINATES:alignment:{0}:{1}".format(self.__class__.NAME, self.__class__.version)]
                    feature.qualifiers["rpt_family"] = ["CRISPR"]
                    feature.qualifiers["rpt_type"] = ["direct"]
                    feature.qualifiers["rpt_unit_seq"] = [repeat.lower()]
                    D.setdefault(feature.seq_id, []).append(feature)
        return D

if __name__ == '__main__':
    logger = getLogger(__name__)

    logger.setLevel(DEBUG)

    handler = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)
    logger.debug("test")
    tool = CRT()
    print(tool.version)