#! /usr/bin/env python
# coding: UTF8

from .base_tools import Tool


class Rpsbproc(Tool):
    version = None
    NAME = "rpsbproc"
    VERSION_CHECK_CMD = ["rpsbproc", "-h"]
    VERSION_PATTERN = r"Post-RPSBLAST Processing Utility v(.+)\n"

    def __init__(self, options=None):
        super(Rpsbproc, self).__init__(options=options)
        self.rpsbproc_data = self.options.get("rpsbproc_data")

    def get_command(self, rpsblast_result, result_file):
        return ["rpsbproc", "-i", rpsblast_result, "-o", result_file, "-d", self.rpsbproc_data, "-q", "2>&1"]
        # return ["rpsbproc", "-i", rpsblast_result, "-o", result_file, "-q", "2>&1"]
        # By default log will be written std-err, so output must be redirected to std-out.
        # rpsbproc - i rpsblast.out -o rpsbproc.out - d /Users/ytanizaw/project/dfast/newpipe/tools/bin/rpsbproc_data - q

if __name__ == '__main__':
    pass
