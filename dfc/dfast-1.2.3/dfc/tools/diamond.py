#! /usr/bin/env python
# coding: UTF8

from .base_tools import Aligner


class Diamond(Aligner):
    """
    Diamond

    Tool type: Protein-protein alignment
    URL: 
    REF: 

    diamond makedb --in DFAST-default.faa --db DFAST-default
    # Will output the Diamond's database file "DFAST-default.dmnd"
    
    diamond blastp --db DFAST-default.dmnd --query protein.faa -o dmnd.out.tsv --threads 1 --max-hsps 1 --max-target-seqs 1 --evalue 1e-5 
    """
    version = None
    NAME = "Diamond"
    VERSION_CHECK_CMD = ["diamond", "version"]
    VERSION_PATTERN = r"diamond version (.+)"

    def format_db_command(self, db_fasta_file, db_name):
        return ["diamond", "makedb", "--in", db_fasta_file, "--db", db_name]

    def get_command(self, query_file, db_name, result_file):
        return ["diamond", "blastp", "--query", query_file, "--db", db_name, "--out", result_file, "--threads 1 --max-hsps 1 --max-target-seqs 1"]

    def get_self_alignment_command(self, query_file, db_name, result_file):
        return ["diamond", "blastp", "--query", query_file, "--db", db_name, "--out", result_file, "--threads 1 --max-hsps 1 --max-target-seqs 100"]

if __name__ == '__main__':
    from logging import getLogger, INFO, DEBUG, StreamHandler

    logger = getLogger("app")

    logger.setLevel(DEBUG)

    handler = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)

    tool = Diamond()

