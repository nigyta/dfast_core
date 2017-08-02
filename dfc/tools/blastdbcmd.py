#! /usr/bin/env python
# coding: UTF8

from .base_tools import Tool


class Blastdbcmd(Tool):
    """
    Blastp from NCBI-BLAST+ package.

    Tool type: 
    URL: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
    REF: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421

    """
    version = None
    NAME = "Blastdbcmd"
    VERSION_CHECK_CMD = ["blastdbcmd", "-version"]
    VERSION_PATTERN = r"blastdbcmd: (.+)\+"

    def get_command(self, entry_file, db_name, result_file, error_log_file):
        # blastdbcmd -db uniprot_sprot.fasta -entry_batch list.txt -out a.fasta 2> e.log
        return ["blastdbcmd", "-db", db_name, "-entry_batch", entry_file, "-out", result_file, "2>", error_log_file]



