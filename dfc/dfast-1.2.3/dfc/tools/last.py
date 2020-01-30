#! /usr/bin/env python
# coding: UTF8


from .base_tools import Tool


class Lastdb(Tool):
    """
    http://last.cbrc.jp

    """
    version = None
    NAME = "lastdb"
    VERSION_CHECK_CMD = ["lastdb", "-V"]
    VERSION_PATTERN = r"^lastdb (.+)$"

    def __init__(self, options=None):
        super(Lastdb, self).__init__(options=options)
        if options is None:
            options = {}
        # self.evalue_cutoff = options.get("evalue_cutoff", 1e-5)

    def get_command(self, query_file, db_name):
        # ./lastdb -p ../fd_ref /Users/tanizawa/projects/newpipe/dfast_core/OUT/FrameShiftDetection/reference.fasta
        return ["lastdb", "-p", db_name, query_file]

class Lastal(Tool):
    version = None
    NAME = "lastal"
    VERSION_CHECK_CMD = ["lastal", "-V"]
    VERSION_PATTERN = r"^lastal (.+)$"

    def __init__(self, options=None):
        if options is None:
            options = {}
        super(Lastal, self).__init__(options=options)
        self.transl_table = options.get("transl_table", 11)
        # self.genetic_code_file = options.get("genetic_code_file", "")
        # self.evalue_cutoff = options.get("evalue_cutoff", 1e-5)

    def get_command(self, query_file, db_name, result_file):
        # ./lastal -f MAF -F15 fd_ref fd_query0.fasta > result.out
        if self.transl_table == 4 or self.transl_table == 25:
            return ["lastal4", "-f MAF", "-F 15", db_name, query_file, ">", result_file]
        else:
            return ["lastal", "-f MAF", "-F 15", db_name, query_file, ">", result_file]

if __name__ == '__main__':
    from logging import getLogger, INFO, DEBUG, StreamHandler

    logger = getLogger(__name__)

    logger.setLevel(DEBUG)

    handler = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)

    Lastdb()