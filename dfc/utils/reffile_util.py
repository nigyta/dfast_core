#! /usr/bin/env python
# coding: UTF8

import os
import sys
import re
from logging import getLogger, DEBUG, INFO, StreamHandler

from ..tools.ghostx import Ghostx
from ..tools.ghostz import Ghostz
from ..tools.diamond import Diamond
from ..tools.blastp import Blastp
from ..tools.hmmer import Hmmer_hmmpress
from ..models.protein import Protein

from Bio import SeqIO


logger = getLogger(__name__)


def check_db_file(db_name, aligner):
    ref_file = db_name + ".ref"

    if not os.path.exists(ref_file):
        logger.error(
            "Reference file ({}) does not exist. Aborting...".format(ref_file))
        exit(1)
    if aligner == "blastp":
        file_ext = ".pin"
        file_name = db_name + file_ext
        if not os.path.exists(file_name):
            logger.warning("BLAST index files do not exist.")
            prepare_blast_database(ref_file)
    elif aligner == "ghostx":
        file_ext = ".inf"
        file_name = db_name + file_ext
        if not os.path.exists(file_name):
            logger.warning("GHOSTX index files do not exist.")
            prepare_ghostx_database(ref_file)
    elif aligner == "ghostz":
        pass
    elif aligner == "diamond":
        file_ext = ".dmnd"
        file_name = db_name + file_ext
        if not os.path.exists(file_name):
            logger.warning("Diamond index file does not exist.")
            prepare_database_dmnd(ref_file)
    else:
        logger.error(
            "Database file ({}) does not exist. Aborting...".format(file_name))
        exit(1)


def check_hmm_db_file(db_file):
    if not os.path.exists(db_file):
        logger.error(
            "HMM Profile ({}) does not exist. Aborting...".format(db_file))
        exit(1)
    pressed_hmm_file = db_file + ".h3i"
    if not os.path.exists(pressed_hmm_file):
        logger.warning("HMM compressed files do not exist.")
        run_hmmpress(db_file)


def read_db_attributes(ref_file_name):
    '''Read the first line of a reference file'''
    pat_attribute = re.compile(r"\[(.+?)=(.*?)\]")
    with open(ref_file_name) as f:
        line = f.readline().strip()
    attributes = pat_attribute.findall(line)
    return dict(attributes)


def fasta2dfast(input_file, output_file=None, no_header=False, source_db=None, attributes=None):
    if output_file is None:
        fh = sys.stdout
    else:
        dir_name = os.path.dirname(output_file)
        if dir_name and not os.path.exists(dir_name):
            os.makedirs(dir_name)
        fh = open(output_file, "w")
    D = Protein.read_from_fasta(input_file, parser_type="auto")
    if not no_header:
        if attributes:
            fh.write("# {}\n".format(attributes))
        else:
            fh.write("#\n")
        fh.write(
            "#id\tdescription\tgene\tEC_number\tflag\torganism\tsource_DB\tsequence\n")
    for protein in D.values():
        if source_db is None:
            infer_source_db = False
        else:
            if source_db == "auto":
                infer_source_db = True
            else:
                infer_source_db = False
                protein.source_db = source_db
        fh.write(protein.to_tsv(infer_source_db=infer_source_db))
    fh.close()


def dfast2fasta(input_file, output_file=None, with_description=False):

    D = Protein.read_from_dfast_reference(input_file)
    Buffer = ""
    for key, protein in D.items():
        if key.startswith("_"):
            continue  # _ means special attributes.
        if protein.id == "" or protein.sequence == "":
            continue  # skip empty line
        Buffer += protein.to_fasta(with_description=with_description)
    if output_file is None:
        logger.info("Converting DFAST reference '{0}' to FASTA".format(
            input_file))
        fh = sys.stdout
    else:
        logger.info("Converting DFAST reference '{0}' to FASTA '{1}'".format(
            input_file, output_file))

        dir_name = os.path.dirname(output_file)
        if dir_name and not os.path.exists(dir_name):
            os.makedirs(dir_name)
        fh = open(output_file, "w")
    fh.write(Buffer)
    fh.close()
    return output_file


def genbank2dfast(input_file, output_file=None, no_header=False, attributes=None):
    if output_file is None:
        fh = sys.stdout
    else:
        dir_name = os.path.dirname(output_file)
        if dir_name and not os.path.exists(dir_name):
            os.makedirs(dir_name)
        fh = open(output_file, "w")
    D = Protein.read_from_genbank(input_file)
    if not no_header:
        if attributes:
            fh.write("# {}\n".format(attributes))
        else:
            fh.write("#\n")
        fh.write(
            "#id\tdescription\tgene\tEC_number\tflag\torganism\tsource_DB\tsequence\n")
    for protein in D.values():
        infer_source_db = False
        fh.write(protein.to_tsv(infer_source_db=infer_source_db))
    fh.close()


def create_blast_idx(file_name, db_name=None):
    if db_name is None:
        db_name = file_name
    blastp = Blastp()
    logger.info("Creating index files for BLAST. {0}".format(db_name))
    logger.setLevel(DEBUG)

    blastp.format_db(file_name, db_name)
    logger.setLevel(INFO)
    logger.info("Done")


def create_ghostx_idx(file_name, db_name=None):
    if db_name is None:
        db_name = file_name
    ghostx = Ghostx()
    logger.info("Creating index files for GHOSTX. {0}".format(db_name))
    logger.setLevel(DEBUG)
    ghostx.format_db(file_name, db_name)
    logger.setLevel(INFO)
    logger.info("Done")


def prepare_database(file_name):
    base_name, _ext = (os.path.splitext(file_name))
    output_file = base_name + ".faa"
    fasta_file = dfast2fasta(file_name, output_file)
    create_ghostx_idx(fasta_file, base_name)
    create_blast_idx(fasta_file, base_name)


def prepare_ghostx_database(file_name):
    base_name, _ext = (os.path.splitext(file_name))
    output_file = base_name + ".faa"
    fasta_file = dfast2fasta(file_name, output_file)
    create_ghostx_idx(fasta_file, base_name)


def prepare_blast_database(file_name):
    base_name, _ext = (os.path.splitext(file_name))
    output_file = base_name + ".faa"
    fasta_file = dfast2fasta(file_name, output_file)
    create_blast_idx(fasta_file, base_name)


def run_hmmpress(file_name):
    hmmpress = Hmmer_hmmpress()
    logger.info(
        "Preparing HMM database using hmmpress. [{0}]".format(file_name))
    logger.setLevel(DEBUG)
    hmmpress.run(file_name)
    logger.setLevel(INFO)
    logger.info("Done")


def prepare_database_dmnd(file_name):
    base_name, _ext = (os.path.splitext(file_name))
    output_file = base_name + ".faa"
    # logger.info("Converting DFAST reference '{0}' to FASTA '{1}'".format(
        # file_name, output_file))
    fasta_file = dfast2fasta(file_name, output_file)
    diamond = Diamond()
    logger.info(
        "Preparing a database for Diamond. (DB_FILE_NAME={0}.dmnd)".format(base_name))
    logger.setLevel(DEBUG)
    diamond.format_db(fasta_file, base_name)
    logger.setLevel(INFO)
    logger.info("Done")
