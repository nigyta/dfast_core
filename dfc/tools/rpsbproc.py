#! /usr/bin/env python
# coding: UTF8

from .base_tools import Tool


class Rpsbproc(Tool):
    version = None
    NAME = "rpsbproc"
    # rpsbproc 0.5 or later is required: it reads the ASN.1 archive output of rpsblast (-outfmt 11).
    # Older versions (e.g. v0.11) do not print "rpsbproc: <version>", so they fail the version check.
    VERSION_CHECK_CMD = ["rpsbproc", "-version", "2>&1"]
    VERSION_PATTERN = r"(?m)^rpsbproc: (\d+\.\d+\S*)"
    VERSION_ERROR_MSG = ("rpsbproc 0.5 or later is required. Install it from Bioconda (conda install -c bioconda rpsbproc) "
                         "or https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/rpsbproc/current/")

    def __init__(self, options=None):
        super(Rpsbproc, self).__init__(options=options)
        self.rpsbproc_data = self.options.get("rpsbproc_data")

    def get_command(self, rpsblast_result, result_file):
        # rpsblast_result is an ASN.1 archive (-outfmt 11). With XML input and -q, rpsbproc 0.5 writes an empty file.
        # By default log will be written std-err, so output must be redirected to std-out.
        return ["rpsbproc", "-i", rpsblast_result, "-o", result_file, "-d", self.rpsbproc_data, "-q", "2>&1"]

if __name__ == '__main__':
    pass
