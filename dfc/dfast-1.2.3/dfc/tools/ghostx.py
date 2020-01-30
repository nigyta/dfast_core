#! /usr/bin/env python
# coding: UTF8

from .base_tools import Aligner


class Ghostx(Aligner):
    """
    Ghostz by Suzuki, S. et al.

    Tool type: Protein-protein alignment
    URL: http://www.bi.cs.titech.ac.jp/ghostz/
    REF: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421


    """
    version = None
    NAME = "GHOSTX"
    VERSION_CHECK_CMD = ["ghostx", "-h"]
    VERSION_PATTERN = r"GHOSTX \- homology search tool\. version (.+)"

    def format_db_command(self, db_fasta_file, db_name):
        return ["ghostx", "db", "-i", db_fasta_file, "-o", db_name]

    def get_command(self, query_file, db_name, result_file):
        return ["ghostx", "aln", "-i", query_file, "-d", db_name, "-o", result_file, "-b 1 -v 1"]

    def get_self_alignment_command(self, query_file, db_name, result_file):
        return ["ghostx", "aln", "-i", query_file, "-d", db_name, "-o", result_file, "-b 100 -v 1"]

if __name__ == '__main__':
    from logging import getLogger, INFO, DEBUG, StreamHandler

    logger = getLogger("app")

    logger.setLevel(DEBUG)

    handler = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)

    tool = Ghostz()
    tool2 = Ghostz()

