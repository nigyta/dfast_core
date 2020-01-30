#! /usr/bin/env python
# coding: UTF8

from .base_tools import Aligner


class Blastp(Aligner):
    """
    Blastp from NCBI-BLAST+ package.

    Tool type: Protein-protein alignment
    URL: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
    REF: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421

    """
    version = None
    NAME = "Blastp"
    VERSION_CHECK_CMD = ["blastp", "-version"]
    VERSION_PATTERN = r"blastp: (.+)\+"


    def format_db_command(self, db_fasta_file, db_name):
        return ["makeblastdb", "-dbtype", "prot", "-in", db_fasta_file, "-out", db_name, "-hash_index", "-parse_seqids"]

    def get_command(self, query_file, db_name, result_file):

        return ["blastp", "-query", query_file, "-db", db_name, "-out", result_file,
                "-outfmt 6",
                "-num_alignments 1 -max_hsps 1"]

        # modified 2017.4.28
        # "-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'",

    def get_self_alignment_command(self, query_file, db_name, result_file):
        return ["blastp", "-query", query_file, "-db", db_name, "-out", result_file,
                "-outfmt 6",
                "-max_hsps 1"]

    def get_extended_command(self, query_file, db_name, result_file):
        outfmt = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle'"
        return ["blastp", "-query", query_file, "-db", db_name, "-out", result_file,
                "-outfmt", outfmt,
                "-num_alignments 1 -max_hsps 1"]

if __name__ == '__main__':
    from logging import getLogger, INFO, DEBUG, StreamHandler

    logger = getLogger(__name__)

    logger.setLevel(DEBUG)

    handler = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)

    tool = Blastp()
    tool2 = Blastp()

