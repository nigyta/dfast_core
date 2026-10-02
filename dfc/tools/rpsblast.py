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
        # ASN.1 archive (-outfmt 11) is the input format recommended for rpsbproc; XML (-outfmt 5) is deprecated.
        return ["rpsblast", "-query", query_file, "-db", db_name, "-out", result_file, "-outfmt 11",
                "-evalue", str(self.evalue_cutoff)]


if __name__ == '__main__':
    from logging import getLogger, INFO, DEBUG, StreamHandler

    logger = getLogger(__name__)

    logger.setLevel(DEBUG)

    handler = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)

    RPSblast()
