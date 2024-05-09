#! /usr/bin/env python
import sys
import os
import json
import shutil
import subprocess
import argparse

from logging import getLogger, DEBUG, INFO, StreamHandler

CARD_VERSION = os.environ.get("CARD_VERSION", "3.2.9")
VFDB_VERSION = os.environ.get("VFDB_VERSION", "2024-05-03")

app_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
logger = getLogger("")
logger.setLevel(INFO)
logger.addHandler(StreamHandler())


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



app_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
# description = """\
# Reference Utility for Nucleotide databases\n\

#     --protein, --cdd, --hmm: For DFAST reference libraries. 
#         Files will be downloaded to DB root directory by default.
#         DB root can be specified with "--dbroot" option.

#     --assembly, --assembly_fasta: For Reference genomes
#         Reference genome file will be downloaded from NCBI Assembly Database either in GenBank or Fasta format.
#         Files will be written to the current directory or the directory specified with "--out" option.

# """

parser = argparse.ArgumentParser(description="Reference Utility for Nucleotide databases",
                                 usage=None, epilog=None,
                                 # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 add_help=True #, # allow_abbrev=False
                                 # argument_default=argparse.SUPPRESS
                                 # formatter_class=argparse.RawTextHelpFormatter
                                 )

# group_basic = parser.add_argument_group("Basic options")
parser.add_argument("--card", action="store_true", help="Prepare reference data for CARD")
parser.add_argument("--vfdb", action="store_true", help="Prepare reference data for VFDB")
# parser.add_argument("--plasmid", action="store_true", help="Prepare reference data for PlasmidDB")
parser.add_argument("-d", "--dbroot", help="DB root directory (default: APP_ROOT/db.\nFor --protein, --cdd, --hmm, --plasmidfinder. Not allowed with argument --out)", type=str, metavar="PATH")

args = parser.parse_args()

if not any([args.card, args.vfdb]):
    parser.print_help()
    sys.stderr.write(f"\n'--card' or '--vfdb' must be specified.\n\n")
    exit()

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

def write_version_file(db_name, version):
    # db_name should be the same as the -db option of BLAST
    version_file = db_name + ".version"
    logger.info(f"Writing the version file to {version_file} ({version})")
    with open(version_file, "w") as f:
        f.write(version)


def card2ref(card_json_file, db_name):
    D = CARD.parse_card_json_file(card_json_file)
    nucl_fasta_file_name = db_name + ".nucl.fasta"
    prot_fasta_file_name = db_name + ".prot.fasta"
    ref_tsv_file_name = db_name + ".nucl.ref"

    CARD.write_nucl_fasta_file(D, nucl_fasta_file_name)
    CARD.write_prot_fasta_file(D, prot_fasta_file_name)
    CARD.write_tsv_file(D, ref_tsv_file_name)
    make_nucl_blast_db(nucl_fasta_file_name, db_name)
    write_version_file(db_name, CARD_VERSION)

def vfdb2ref(vfdb_nucl_fasta_file, vfdb_prot_fasta_file, db_name):
    D = VFDB.parse_vfdb_fasta_files(vfdb_nucl_fasta_file, vfdb_prot_fasta_file)
    nucl_fasta_file_name = db_name + ".nucl.fasta"
    prot_fasta_file_name = db_name + ".prot.fasta"
    ref_tsv_file_name = db_name + ".nucl.ref"

    VFDB.write_nucl_fasta_file(D, nucl_fasta_file_name)
    VFDB.write_prot_fasta_file(D, prot_fasta_file_name)
    VFDB.write_tsv_file(D, ref_tsv_file_name)
    make_nucl_blast_db(nucl_fasta_file_name, db_name)
    write_version_file(db_name, VFDB_VERSION)

def get_db_root(args):
    """
    Priority of DB_ROOT: command argument (--dbroot), env variable DFAST_DB_ROOT, default
    """
    db_root_env = os.getenv('DFAST_DB_ROOT')

    if args.dbroot:
        db_root = args.dbroot
    elif db_root_env:
        db_root = db_root_env
    else:
        db_root = os.path.join(app_root, "db")
    return db_root

if args.card:
    # Retrieve and convert CARD reference data
    db_root = get_db_root(args)
    tmp_dir = os.path.join(db_root, "card_temp")
    out_dir = os.path.join(db_root, "card")
    logger.info(f"Preparing reference data for CARD in {out_dir} [ver. {CARD_VERSION}]")
    logger.info(f"Created temporary working directory: {tmp_dir}")

    if os.path.exists(tmp_dir):
        logger.info(f'Will delete existing directory: "{tmp_dir}"')
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)
    cmd = f"cd {tmp_dir} && curl -LO https://card.mcmaster.ca/download/0/broadstreet-v{CARD_VERSION}.tar.bz2 && "
    cmd += f"tar xf broadstreet-v{CARD_VERSION}.tar.bz2 && "
    cmd += f"rm broadstreet-v{CARD_VERSION}.tar.bz2"
    logger.info(f'Running command: "{cmd}"')
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()

    card_json_file = os.path.join(tmp_dir, "card.json")
    db_name = os.path.join(out_dir, "CARD")
    if os.path.exists(out_dir):
        logger.info(f'Will delete existing directory: "{out_dir}"')
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    card2ref(card_json_file, db_name)
    shutil.rmtree(tmp_dir)


if args.vfdb:
    # Retrieve and convert VFDB reference data
    db_root = get_db_root(args)
    tmp_dir = os.path.join(db_root, "vfdb_temp")
    out_dir = os.path.join(db_root, "vfdb")
    logger.info(f"Preparing reference data for VFDB in {out_dir} [ver. {VFDB_VERSION}]")
    logger.info(f"Created temporary working directory: {tmp_dir}")

    if os.path.exists(tmp_dir):
        logger.info(f'Will delete existing directory: "{tmp_dir}"')
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)

    # Download nt fasta
    cmd = f"cd {tmp_dir} && curl -LO http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz && gunzip VFDB_setB_nt.fas.gz"
    logger.info(f'Running command: "{cmd}"')
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()

    # Download prot fasta
    cmd = f"cd {tmp_dir} && curl -LO http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz && gunzip VFDB_setB_pro.fas.gz"
    logger.info(f'Running command: "{cmd}"')
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, shell=True)

    out, err = p.communicate()
    vfdb_nucl_fasta_file = os.path.join(tmp_dir, "VFDB_setB_nt.fas")
    vfdb_prot_fasta_file = os.path.join(tmp_dir, "VFDB_setB_pro.fas")
    db_name = os.path.join(out_dir, "VFDB")
    if os.path.exists(out_dir):
        logger.info(f'Will delete existing directory: "{out_dir}"')
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    vfdb2ref(vfdb_nucl_fasta_file, vfdb_prot_fasta_file, db_name)

if __name__ == '__main__':
    pass
    # ref_file_name = "/work/db/plasmid_referece.nucref"
    # db_name = "/work/db/PLASMID_DB"
    # plasmiddb2ref(ref_file_name, db_name)

    # card_json_file = "/work/db/card.json"
    # db_name = "/work/db/CARD"
    # card2ref(card_json_file, db_name)

    # vfdb_nucl_fasta_file = "/work/db/VFDB_setB_nt.fas"
    # vfdb_prot_fasta_file = "/work/db/VFDB_setB_pro.fas"
    # db_name = "db/VFDB"
    # vfdb2ref(vfdb_nucl_fasta_file, vfdb_prot_fasta_file, db_name)
