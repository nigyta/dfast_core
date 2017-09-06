#! /usr/bin/env python
# coding: UTF8


from .base_tools import Tool


class Hmmer_hmmscan(Tool):
    version = None
    NAME = "hmmscan"
    VERSION_CHECK_CMD = ["hmmscan", "-h"]
    VERSION_PATTERN = r"# HMMER (.+) \("

    def __init__(self, options=None):
        super(Hmmer_hmmscan, self).__init__(options=options)
        if options is None:
            options = {}
        self.evalue_cutoff = options.get("evalue_cutoff", 1e-5)

    def get_command(self, query_file, db_name, result_file):
        return ["hmmscan", "--noali", "--notextw", "--cpu", "1", "-E", str(self.evalue_cutoff), "--tblout", result_file,
                db_name, query_file]


class Hmmer_hmmpress(Tool):
    version = None
    NAME = "hmmpress"
    VERSION_CHECK_CMD = ["hmmpress", "-h"]
    VERSION_PATTERN = r"# HMMER (.+) \("

    def __init__(self, options=None):
        super(Hmmer_hmmpress, self).__init__(options=options)
        # if options is None:
        #     options = {}

    def get_command(self, input_file):
        return ["hmmpress", "-f", input_file]

    def run(self, input_file):
        cmd = self.get_command(input_file)
        self.executeCommand(cmd, shell=False)


class Hmmer_hmmsearch(Tool):
    version = None
    NAME = "hmmsearch"
    VERSION_CHECK_CMD = ["hmmsearch", "-h"]
    VERSION_PATTERN = r"# HMMER (.+) \("

    def __init__(self, options=None):
        super(Hmmer_hmmsearch, self).__init__(options=options)
        if options is None:
            options = {}
        self.evalue_cutoff = options.get("evalue_cutoff", 1e-5)

    def get_command(self, sequence_file, hmm_profile, result_file):
        return ["hmmsearch", "--noali", "--notextw", "--cpu", "1", "-E", str(self.evalue_cutoff), "--tblout", result_file,
                "-o", result_file + ".raw", hmm_profile, sequence_file]


if __name__ == '__main__':
    from logging import getLogger, INFO, DEBUG, StreamHandler

    logger = getLogger(__name__)

    logger.setLevel(DEBUG)

    handler = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)

    Hmmer_hmmscan()
