#! /usr/bin/env python
# coding: UTF8


from .base_tools import Tool


class RPSblast(Tool):
    version = None
    NAME = "RPSblast"
    VERSION_CHECK_CMD = ["rpsblast", "-version"]
    VERSION_PATTERN = r"rpsblast: (.+)\+"

    def __init__(self, options=None):
        super(RPSblast, self).__init__(options=options)
        self.evalue_cutoff = options.get("evalue_cutoff", 1e-6)

    def get_command(self, query_file, db_name, result_file):
        return ["rpsblast", "-query", query_file, "-db", db_name, "-out", result_file, "-outfmt 5",
                "-evalue", str(self.evalue_cutoff)]
        # rpsblast -query query0.fasta -db Cog -out rpsblast.out -outfmt 5 -evalue 1e-5


if __name__ == '__main__':
    from logging import getLogger, INFO, DEBUG, StreamHandler

    logger = getLogger(__name__)

    logger.setLevel(DEBUG)

    handler = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)

    RPSblast()
