#! /usr/bin/env python
# coding: UTF8

import os
import sys
import re
from logging import getLogger, DEBUG, INFO, StreamHandler
from Bio import SeqIO
# from .dbfile_util import prepare_ghostx_database


logger = getLogger(__name__)


pat_ncbi = re.compile(r"^(.+) \[(.+)\]$")

def ncbi_parser(s_id, title):
    if title.startswith(s_id):
        title = title.replace(s_id, "").strip()
    if "MULTISPECIES: " in title:
        title = title.replace("MULTISPECIES: ", "")
    if "|" in s_id:
        s_id = s_id.split("|")[1]
    m = pat_ncbi.match(title)
    if m:
        product, organism = m.groups()
    else:
        product, organism = title, ""
    gene = ""
    source_db = "RefSeq" if s_id[2] == "_" else "INSD"
    ec_number = ""
    return s_id, product, organism, gene, source_db, ec_number


def uniprot_parser(s_id, title):
    if title.startswith(s_id):
        title = title.replace(s_id, "").strip()
    if "|" in s_id:
        s_id = s_id.split("|")[1]
    source_db = "UniprotKB"

    splitted_title = title.split("=")
    if len(splitted_title) == 1:
        return title, "", "", source_db
    key = "product"
    D = {}
    for i in range(len(splitted_title) - 1):

        words = splitted_title[i].split()
        value = " ".join(words[:-1])
        D[key] = value
        key = words[-1]
    else:
        D[key] = splitted_title[-1]

    product = D["product"]
    organism = D.get("OS") or D.get("Tax", "")
    gene = D.get("GN", "")
    ec_number = ""
    return s_id, product, organism, gene, source_db, ec_number


def plain_fasta_parser(s_id, title):
    if title.startswith(s_id):
        title = title.replace(s_id, "").strip()
    if "|" in s_id:
        s_id = s_id.split("|")[1]
    product = title
    gene = ""
    source_db = ""
    organism = ""
    ec_number = ""
    return s_id, product, organism, gene, source_db, ec_number


    
def prokka_fasta_parser(s_id, title):
    if title.startswith(s_id):
        title = title.replace(s_id, "").strip()
    if "~~~" in title:
        ec_number, gene, product = title.split("~~~")
    else:
        ec_number, gene, product = "", "", title
    source_db = ""
    organism = ""
    return s_id, product, organism, gene, source_db, ec_number


def auto_fasta_parser(s_id, title):
    if "~~~" in title:
        return prokka_fasta_parser(s_id, title)
    elif " OS=" in title or " Tax=" in title or " n=" in title or " GN=" in title:
        return uniprot_parser(s_id, title)
    elif title.endswith("]") and " [" in title:
        return ncbi_parser(s_id, title)
    else:
        return plain_fasta_parser(s_id, title)


fasta_parsers = {
    "ncbi": ncbi_parser,
    "uniprot": uniprot_parser,
    "plain": plain_fasta_parser,
    "prokka": prokka_fasta_parser,
    "auto": auto_fasta_parser
}


def get_source_db(prot_id):
    if "." in prot_id:
        if prot_id[2] == "_":
            return "RefSeq"
        else:
            return "INSD"
    elif len(prot_id) == 6 or len(prot_id) == 10:
        return "UniProtKB"
    else:
        return ""






class RefUtil(object):

    def __init__(self):
        pass

    def read_protein_from_gbk(self, file_name):
        pass
        # R = list(SeqIO.parse(open(gbkfile), "genbank"))

if __name__ == '__main__':
    pass
