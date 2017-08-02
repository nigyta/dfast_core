#! /usr/bin/env python
# coding: UTF8

import os
import re
from logging import getLogger
from Bio import SeqIO

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
    organism = D.get("OS", "")
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
    elif " OS=" in title or " n=" in title or " GN=" in title:
        return uniprot_parser(s_id, title)
    elif title.endswith("]") and " [" in title:
        return ncbi_parser(s_id, title)
    else:
        return plain_fasta_parser(s_id, title)


fasta_parsers = {
    "ncbi": ncbi_parser,
    "uniprot": uniprot_parser,
    "plain": plain_fasta_parser,
    "prokka": auto_fasta_parser,
    "auto": auto_fasta_parser
}


def check_db_file(db_name, aligner):
    ref_file = db_name + ".ref"

    if not os.path.exists(ref_file):
        logger.error("Reference file ({}) does not exist. Aborting...".format(ref_file))
        exit(1)
    if aligner == "blastp":
        file_extensions = [".pin", ".phr", ".psq"]
    elif aligner in ["ghostx", "ghostz"]:
        file_extensions = [".inf"]
    else:
        file_extensions = []
    for file_ext in file_extensions:
        file_name = db_name + file_ext
        if not os.path.exists(file_name):
            logger.error("Database file ({}) does not exist. Aborting...".format(file_name))
            exit(1)


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


def read_db_attributes(ref_file_name):
    '''Read the first line of a reference file'''
    pat_attribute = re.compile(r"\[(.+?)=(.*?)\]")
    with open(ref_file_name) as f:
        line = f.readline().strip()
    attributes = pat_attribute.findall(line)
    return dict(attributes)


class RefUtil(object):

    def __init__(self):
        pass

    def read_protein_from_gbk(file_name):
        R = list(SeqIO.parse(open(gbkfile), "genbank"))

if __name__ == '__main__':
    pass
