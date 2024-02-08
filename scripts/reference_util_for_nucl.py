#! /usr/bin/env python
import sys
import os
import json
# import platform
from logging import getLogger, DEBUG, INFO, StreamHandler

app_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
logger = getLogger("")
logger.setLevel(INFO)

# srcPath = os.path.join(os.path.dirname(__file__), "..", "dfc")
sys.path.append(app_root)
from dfc.utils.path_util import set_binaries_path
# from dfc.tools.ghostx import Ghostx
# from dfc.tools.ghostz import Ghostz
# from dfc.tools.diamond import Diamond
from dfc.tools.blastn import Blastn
from dfc.models.nucref import NucRef
from dfc.models.card import CARD
from dfc.models.vfdb import VFDB
# from dfc.tools.hmmer import Hmmer_hmmpress
# from dfc.models.protein import Protein
# from dfc.utils.reffile_util import fasta2dfast, dfast2fasta, genbank2dfast, prepare_database, prepare_database_dmnd, run_hmmpress
set_binaries_path(app_root)




def plasmiddb2ref(ref_file_name, db_name):
    D = NucRef.read_ref_file(ref_file_name)
    nucl_fasta_file_name = db_name + ".nucl.fasta"
    prot_fasta_file_name = db_name + ".prot.fasta"
    ref_tsv_file_name = db_name + ".nucl.ref"

    NucRef.write_nucl_fasta_file(D, nucl_fasta_file_name)
    NucRef.write_prot_fasta_file(D, prot_fasta_file_name)
    NucRef.write_tsv_file(D, ref_tsv_file_name)
    make_nucl_blast_db(nucl_fasta_file_name, db_name)


def make_nucl_blast_db(fasta_file_name, db_name):
    logger.setLevel(DEBUG)
    logger.info("Creating index files for BLASTN. {0}".format(db_name))
    blastn = Blastn()
    blastn.format_db(fasta_file_name, db_name)

def card2ref(card_json_file, db_name):
    D = CARD.parse_card_json_file(card_json_file)
    nucl_fasta_file_name = db_name + ".nucl.fasta"
    prot_fasta_file_name = db_name + ".prot.fasta"
    ref_tsv_file_name = db_name + ".nucl.ref"

    CARD.write_nucl_fasta_file(D, nucl_fasta_file_name)
    CARD.write_prot_fasta_file(D, prot_fasta_file_name)
    CARD.write_tsv_file(D, ref_tsv_file_name)
    make_nucl_blast_db(nucl_fasta_file_name, db_name)


def vfdb2ref(vfdb_nucl_fasta_file, vfdb_prot_fasta_file, db_name):
    D = VFDB.parse_vfdb_fasta_files(vfdb_nucl_fasta_file, vfdb_prot_fasta_file)
    nucl_fasta_file_name = db_name + ".nucl.fasta"
    prot_fasta_file_name = db_name + ".prot.fasta"
    ref_tsv_file_name = db_name + ".nucl.ref"

    VFDB.write_nucl_fasta_file(D, nucl_fasta_file_name)
    VFDB.write_prot_fasta_file(D, prot_fasta_file_name)
    VFDB.write_tsv_file(D, ref_tsv_file_name)
    make_nucl_blast_db(nucl_fasta_file_name, db_name)


if __name__ == '__main__':

    # ref_file_name = "/work/db/plasmid_referece.nucref"
    # db_name = "/work/db/PLASMID_DB"
    # plasmiddb2ref(ref_file_name, db_name)

    # card_json_file = "/work/db/card.json"
    # db_name = "/work/db/CARD"
    # card2ref(card_json_file, db_name)

    vfdb_nucl_fasta_file = "/work/db/VFDB_setB_nt.fas"
    vfdb_prot_fasta_file = "/work/db/VFDB_setB_pro.fas"
    db_name = "db/VFDB"
    vfdb2ref(vfdb_nucl_fasta_file, vfdb_prot_fasta_file, db_name)
