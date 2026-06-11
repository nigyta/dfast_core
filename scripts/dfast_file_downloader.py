#! /usr/bin/env python
# coding: UTF8

import argparse
import sys
import os
import gzip
from logging import getLogger, StreamHandler, INFO
from ftplib import FTP
import tarfile
import shutil
import subprocess
import tempfile

version = sys.version_info.major
if version == 3:
    from urllib import request
else:
    from six.moves.urllib import request

app_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
logger = getLogger(__name__)
logger.setLevel(INFO)
logger.addHandler(StreamHandler())

sys.path.append(app_root)
from dfc.utils.reffile_util import prepare_database, run_hmmpress
from dfc.utils.path_util import set_binaries_path

set_binaries_path(app_root)

try:
    # This part was added to avoid error under a specific environment.
    import ssl
    ssl._create_default_https_context = ssl._create_unverified_context
except Exception as e:
    logger.warning("SSL configuration failed. Will continue processing without SSL configuration.")

app_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))

ncbi_ftp_server = "ftp.ncbi.nlm.nih.gov"
cdd_directory = "/pub/mmdb/cdd/little_endian/"
host_dfast = "https://dfast.ddbj.nig.ac.jp"

db_urls = {
    "dfast": host_dfast + "/dfc/distribution/DFAST-default.ref.gz",
    "bifido": host_dfast + "/dfc/distribution/DFAST-BIFIDO.ref.gz",
    "cyanobase": host_dfast + "/dfc/distribution/cyanobase_aa.ref.gz",
    "ecoli": host_dfast + "/dfc/distribution/DFAST-ECOLI.ref.gz",
    "lab": host_dfast + "/dfc/distribution/DFAST-LAB.ref.gz",
    "pylori": host_dfast + "/dfc/distribution/DFAST-Hpylori.ref.gz",
}

hmm_urls = {
    "Pfam": "ftp://ftp.ebi.ac.uk//pub/databases/Pfam/releases/Pfam37.0/Pfam-A.hmm.gz",
    "dbCAN": "http://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V12.txt",
    "TIGR": "https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/TIGRFAMs_15.0_HMM.LIB.gz"
}

cdd_url = "https://ftp.ncbi.nlm.nih.gov//pub/mmdb/cdd/little_endian/DBNAME_LE.tar.gz"

description = """\
DFAST file downloader\n\

    --protein, --cdd, --hmm: For DFAST reference libraries. 
        Files will be downloaded to DB root directory by default.
        DB root can be specified with "--dbroot" option.

    --assembly, --assembly_fasta: For Reference genomes
        Reference genome file will be downloaded from NCBI Assembly Database either in GenBank or Fasta format.
        Files will be written to the current directory or the directory specified with "--out" option.

"""
parser = argparse.ArgumentParser(description=description,
                                 usage=None, epilog=None,
                                 # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 add_help=True, # allow_abbrev=False
                                 # argument_default=argparse.SUPPRESS
                                 formatter_class=argparse.RawTextHelpFormatter
                                 )

# group_basic = parser.add_argument_group("Basic options")
parser.add_argument("--protein", nargs='+', choices=list(db_urls.keys()),
                         help="DFAST reference databases. [{}]".format("|".join(list(db_urls.keys()))), metavar="STR")
parser.add_argument("--cdd", nargs='+', choices=["Cdd", "Cdd_NCBI", "Cog", "Kog", "Pfam", "Prk", "Smart", "Tigr"],
                         help="Preformatted RPS-BLAST database. [Cdd|Cdd_NCBI|Cog|Kog|Pfam|Prk|Smart|Tigr]", metavar="STR")
parser.add_argument("--hmm", nargs='+', choices=["Pfam", "TIGR", "dbCAN"],
                         help="Preformatted RPS-BLAST database. [Pfam|TIGR|dbCAN]", metavar="STR")
parser.add_argument("--assembly", nargs='*', metavar="ACCESSION",
                         help="Accession(s) for NCBI Assembly DB. eg. GCF_000091005.1 GCA_000008865.1")
parser.add_argument("--assembly_fasta", nargs='*', metavar="ACCESSION",
                         help="Accession(s) for NCBI Assembly DB. eg. GCF_000091005.1 GCA_000008865.1")
parser.add_argument("--plasmidfinder", action="store_true",
                         help="Reference data for PlasmidFinder")
parser.add_argument("--mefinder", action="store_true",
                         help="Install MobileElementFinder and build its MGEdb BLAST index")
parser.add_argument("--no_indexing", action="store_true",
                         help="Do not perform database indexing")
group_out = parser.add_mutually_exclusive_group()
group_out.add_argument("-o", "--out", help="Output directory (default: current directory.\nFor --assembly, --assembly_fasta. Not allowed with argument --dbroot)", type=str, metavar="PATH")
group_out.add_argument("-d", "--dbroot", help="DB root directory (default: APP_ROOT/db.\nFor --protein, --cdd, --hmm, --plasmidfinder, --mefinder. Not allowed with argument --out)", type=str, metavar="PATH")

args = parser.parse_args()


if all(x is None for x in [args.protein, args.cdd, args.hmm, args.assembly, args.assembly_fasta, args.plasmidfinder]) and not args.mefinder:
    parser.print_help()
    exit()

def retrieve_dfast_reference(db_name, out_dir="."):
    target_url = db_urls[db_name]
    target_file = os.path.basename(target_url)
    output_file = os.path.join(out_dir, target_file)
    request.urlretrieve(target_url, output_file)
    logger.info("\tTarget URL: {}".format(target_url))
    return output_file


def retrieve_hmm(db_name, out_dir="."):
    target_url = hmm_urls[db_name]
    target_file = os.path.basename(target_url)
    output_file = os.path.join(out_dir, target_file)
    request.urlretrieve(target_url, output_file)
    logger.info("\tDownloading {}".format(target_url))
    return output_file

# deprecated
def retrieve_cdd_ftp(db_name, out_dir="."):
    ftp = FTP(host=ncbi_ftp_server)
    logger.info("\tLogging in to the FTP server. {}".format(ncbi_ftp_server + cdd_directory))
    ftp.login()
    # ftp.cwd(directory)

    target_file = "{}_LE.tar.gz".format(db_name)
    logger.info("\tDownloading {}".format(ncbi_ftp_server + cdd_directory + target_file))
    output_file = os.path.join(out_dir, target_file)
    with open(output_file, "wb") as f:
        ftp.retrbinary("RETR " + cdd_directory + target_file, f.write)
    ftp.quit()
    return output_file

def retrieve_cdd(db_name, out_dir="."):
    if db_name != "Cog":
        logger.warning("Currently, only 'cog' is supported for CDD database.")
        exit(1)

    # target_url = cdd_url.replace("DBNAME", db_name)
    target_url = "https://ddbj.nig.ac.jp/public/software/dfast/cog.tar.gz"
    logger.warning("Cog reference data will be downloaded from https://ddbj.nig.ac.jp/public/software/dfast/cog.tar.gz")
    target_file = os.path.basename(target_url)
    output_file = os.path.join(out_dir, target_file)
    request.urlretrieve(target_url, output_file)
    logger.info("\tDownloading {}".format(target_url))
    return output_file

def retrieve_assembly(accession, out_dir="."):
    def _get_ftp_directory(accession):
        path1, path2, path3, path4 = accession[0:3], accession[4:7], accession[7:10], accession[10:13]
        return "/".join(["/genomes", "all", path1, path2, path3, path4])

    directory = _get_ftp_directory(accession)
    ftp = FTP(host=ncbi_ftp_server)
    logger.info("\tLogging in to the FTP server. {}".format(ncbi_ftp_server + directory))
    ftp.login()
    ftp.cwd(directory)
    L = ftp.nlst()
    if len(L) == 0:
        logger.warning("\tFile not found. Skip retrieving file for {}".format(accession))
        return None
    asm_name = sorted([x for x in L if x.startswith(accession)])[-1]
    target_file = "/".join([directory, asm_name, asm_name + "_genomic.gbff.gz"])
    logger.info("\tDownloading {}".format(ncbi_ftp_server + directory + "/" + asm_name + "/" + asm_name + "_genomic.gbff.gz"))
    output_file = os.path.join(out_dir, "_".join(asm_name.split("_")[0:2]) + ".gbk.gz")
    with open(output_file, "wb") as f:
        ftp.retrbinary("RETR " + target_file, f.write)
    ftp.quit()
    return output_file

def retrieve_assembly_fasta(accession, out_dir="."):
    def _get_ftp_directory(accession):
        path1, path2, path3, path4 = accession[0:3], accession[4:7], accession[7:10], accession[10:13]
        return "/".join(["/genomes", "all", path1, path2, path3, path4])

    directory = _get_ftp_directory(accession)
    ftp = FTP(host=ncbi_ftp_server)
    logger.info("\tLogging in to the FTP server. {}".format(ncbi_ftp_server + directory))
    ftp.login()
    ftp.cwd(directory)
    L = ftp.nlst()
    if len(L) == 0:
        logger.warning("\tFile not found. Skip retrieving file for {}".format(accession))
        return None
    asm_name = sorted([x for x in L if x.startswith(accession)])[-1]
    target_file = "/".join([directory, asm_name, asm_name + "_genomic.fna.gz"])
    logger.info("\tDownloading {}".format(ncbi_ftp_server + directory + "/" + asm_name + "/" + asm_name + "_genomic.fna.gz"))
    output_file = os.path.join(out_dir, "_".join(asm_name.split("_")[0:2]) + ".fna.gz")
    with open(output_file, "wb") as f:
        ftp.retrbinary("RETR " + target_file, f.write)
    ftp.quit()
    return output_file

plasmidfinder_extra_databases = [
    # {
    #     # figshare dataset: https://figshare.com/articles/dataset/repP_database_v1_0/31883911
    #     "url": "https://ndownloader.figshare.com/files/63271042",
    #     "filename": "repP_database_v2.fsa",
    #     "db_prefix": "repP_database_v2",
    #     "name": "repP",
    #     "description": "repP plasmid replicons",
    #     "fix_headers": False,
    # },
    {
        # figshare dataset: https://figshare.com/articles/dataset/repP_database_fsa/26778175
        "url": "https://ndownloader.figshare.com/files/63270934",
        "filename": "PseudomonasRepDB_v1.fsa",
        "db_prefix": "PseudomonasRepDB_v1",
        "name": "PseudomonasRepDB",
        "description": "PseudomonasRepDB database v1.0",
        "fix_headers": False,
    },
    {
        # figshare dataset: https://figshare.com/articles/dataset/WHRepDB_database_v1_0/31883887
        "url": "https://ndownloader.figshare.com/files/63271000",
        "filename": "WHRepDB_v1.fsa",
        "db_prefix": "WHRepDB_v1",
        "name": "WHRepDB",
        "description": "WHRepDB database v1.0",
        "fix_headers": False,
    },
    {
        # figshare dataset: https://figshare.com/articles/dataset/AcinetobacterPlasmidTyping_database_v3_0/30426412
        # Github: https://github.com/MehradHamidian/AcinetobacterPlasmidTyping
        "url": "https://ndownloader.figshare.com/files/58980199",
        "filename": "AcinetobacterPlasmidTyping_v3.fsa",
        "db_prefix": "AcinetobacterPlasmidTyping_v3",
        "name": "APT DB",
        "description": "AcinetobacterPlasmidTyping database v3",
        # Headers contain spaces (e.g. ">R3-T1_1_NC010605.1 pACICU1_c38"), which
        # blastn truncates at the space, causing ref_acc to be empty in PlasmidFinder v3.
        # Replace spaces with underscores in FASTA headers to fix parsing.
        "fix_headers": True,
    },
]


def retrieve_plasmidfinder_reference(out_dir="."):
    plasmid_db_dir = os.path.join(out_dir, "plasmidfinder_db")
    if os.path.exists(plasmid_db_dir):
        logger.info(f'Will delete existing directory: "{plasmid_db_dir}"')
        shutil.rmtree(plasmid_db_dir)

    # Clone the database repository
    clone_cmd = f"cd {out_dir} && git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git"
    logger.info(f'Cloning PlasmidFinder database: "{clone_cmd}"')
    p = subprocess.Popen(clone_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    if p.returncode != 0 and err:
        logger.error("Failed to clone PlasmidFinder database!")
        logger.error(err.decode("utf8"))
        exit(1)

    # Download extra databases and update config
    _download_plasmidfinder_extra_databases(plasmid_db_dir)
    logger.info("PlasmidFinder database setup completed.")


def retrieve_mefinder_reference(out_dir=".", no_indexing=False):
    mefinder_db_dir = os.path.join(out_dir, "mefinder_db")

    # 1. Install MobileElementFinder if not already available
    check = subprocess.run(["mefinder", "--version"],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if check.returncode == 0:
        logger.info("MobileElementFinder is already installed.")
    else:
        logger.info("Installing MobileElementFinder (clone, relax biopython pin, pip install)...")
        tmpdir = tempfile.mkdtemp()
        clone_dir = os.path.join(tmpdir, "mge_finder")
        if subprocess.run(["git", "clone",
                           "https://bitbucket.org/mhkj/mge_finder.git",
                           clone_dir]).returncode != 0:
            logger.error("Failed to clone MobileElementFinder repository!")
            shutil.rmtree(tmpdir, ignore_errors=True)
            exit(1)
        # Relax the biopython upper-bound pin so MEF works with DFAST's biopython
        setup_cfg = os.path.join(clone_dir, "setup.cfg")
        with open(setup_cfg) as f:
            cfg = f.read()
        new_cfg = cfg.replace("biopython<=1.80", "biopython")
        if new_cfg == cfg:
            logger.warning("biopython pin 'biopython<=1.80' not found in setup.cfg; skipping patch.")
        with open(setup_cfg, "w") as f:
            f.write(new_cfg)
        # setuptools<81 keeps pkg_resources available, which MobileElementFinder imports.
        if subprocess.run([sys.executable, "-m", "pip", "install",
                           clone_dir, "setuptools<81"]).returncode != 0:
            logger.error("Failed to pip install MobileElementFinder!")
            shutil.rmtree(tmpdir, ignore_errors=True)
            exit(1)
        shutil.rmtree(tmpdir, ignore_errors=True)

    if no_indexing:
        logger.info("Skipping MGEdb BLAST index build (--no_indexing).")
        return

    # 2. Build BLAST index from the bundled MGEdb into DB_ROOT/mefinder_db
    if not os.path.exists(mefinder_db_dir):
        os.makedirs(mefinder_db_dir)
    logger.info(f"Building MGEdb BLAST index into: {mefinder_db_dir}")
    if subprocess.run(["mefinder", "index", "--db-path", mefinder_db_dir]).returncode != 0:
        logger.error("Failed to build MobileElementFinder BLAST index!")
        exit(1)

    # 3. Verify index files exist
    required = [os.path.join(mefinder_db_dir, "mge_records." + ext)
                for ext in ("nin", "nhr", "nsq")]
    missing = [f for f in required if not os.path.exists(f)]
    if missing:
        logger.error(f"BLAST index incomplete, missing: {missing}")
        exit(1)
    logger.info(f"MobileElementFinder setup completed. Index at: {mefinder_db_dir}")


def _download_plasmidfinder_extra_databases(plasmid_db_dir):
    # PlasmidFinder v3 uses blastn with -subject <db>.fsa directly,
    # so no pre-built index (makeblastdb/KMA) is needed.
    config_file = os.path.join(plasmid_db_dir, "config")
    for db in plasmidfinder_extra_databases:
        output_file = os.path.join(plasmid_db_dir, db["filename"])
        logger.info(f'\tDownloading extra database: {db["filename"]} from {db["url"]}')
        try:
            request.urlretrieve(db["url"], output_file)
        except Exception as e:
            logger.warning(f'\tFailed to download {db["filename"]}: {e}. Skipping.')
            continue
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            logger.info(f'\tDownloaded {db["filename"]} ({os.path.getsize(output_file)} bytes)')
            if db.get("fix_headers"):
                _fix_fasta_headers(output_file)
                logger.info(f'\tFixed FASTA headers in {db["filename"]} (spaces replaced with underscores)')
            entry = f'{db["db_prefix"]}\t{db["name"]}\t{db["description"]}\n'
            with open(config_file, "a") as f:
                f.write(entry)
        else:
            logger.warning(f'\tDownloaded file is empty: {db["filename"]}. Skipping.')


def _fix_fasta_headers(fasta_file):
    """Reformat FASTA headers for PlasmidFinder v3 compatibility.

    AcinetobacterPlasmidTyping_v3 headers use the format:
      gene_variant_accession extra_info
    but PlasmidFinder expects underscore-delimited gene_variant_note_accession.
    blastn truncates subject headers at spaces, so the accession (position 2) must
    be moved to the last field. We use double-underscore (empty note) to preserve it:
      >R3-T1_1_NC010605.1 pACICU1_c38  →  >R3-T1_1__NC010605.1
    """
    with open(fasta_file, "r") as f:
        lines = f.readlines()
    with open(fasta_file, "w") as f:
        for line in lines:
            if line.startswith(">") and " " in line:
                header = line[1:].strip()
                before_space = header.split(" ")[0]  # e.g. "R3-T1_1_NC010605.1"
                parts = before_space.split("_", 2)   # ["R3-T1", "1", "NC010605.1"]
                if len(parts) == 3:
                    gene, variant, accession = parts
                    line = f">{gene}_{variant}__{accession}\n"
                else:
                    line = line.replace(" ", "_")
            f.write(line)


def gunzip_file(input_file, output_file, cleanup=True):
    with open(output_file, "w") as fw:
        with gzip.open(input_file, 'r') as fr:
            fw.write(fr.read().decode("UTF-8"))
    if cleanup:
        os.remove(input_file)

def extract_tar_file(input_file, output_dir, cleanup=True):
    tar = tarfile.open(input_file)
    tar.extractall(path=output_dir)
    tar.close()
    if cleanup:
        os.remove(input_file)

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


if args.protein:
    db_root = get_db_root(args)
    out_dir = os.path.join(db_root, "protein")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    logger.info("Downloading DFAST reference. Files will be written into '{}'".format(out_dir))
    for db_name in args.protein:
        # db_name = "A"
        retrieved_file = retrieve_dfast_reference(db_name, out_dir)
        if retrieved_file:
            output_file = retrieved_file.replace(".ref.gz", ".ref")
            gunzip_file(retrieved_file, output_file)
            logger.info("\tDownloaded to {}".format(os.path.abspath(output_file)))
            if not args.no_indexing:
                prepare_database(output_file)

if args.cdd:
    db_root = get_db_root(args)
    out_dir = os.path.join(db_root, "cdd")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    logger.info("Trying to fetch RPS-BLAST databases from NCBI-CDD. Files will be written into '{}'".format(out_dir))
    for db_name in args.cdd:
        logger.info("Downloading a database for {}...".format(db_name))
        retrieved_file = retrieve_cdd(db_name, out_dir)
        if retrieved_file:
            logger.info("\tDownloaded to {}".format(os.path.abspath(retrieved_file)))
            extract_tar_file(retrieved_file, out_dir)
            logger.info("\tUnarchived {}".format(retrieved_file))

if args.hmm:
    # print(args.hmm)
    db_root = get_db_root(args)
    out_dir = os.path.join(db_root, "hmm")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    logger.info("Ttyring to fetch profile-HMM databases for {0}. Files will be written into '{1}'".format(",".join(args.hmm), out_dir))
    for db_name in args.hmm:
        logger.info("Downloading HMM database for {}...".format(db_name))
        retrieved_file = retrieve_hmm(db_name, out_dir)
        if retrieved_file:
            logger.info("\tDownloaded to {}".format(os.path.abspath(retrieved_file)))
            if retrieved_file.endswith(".gz"):
                output_file = retrieved_file.replace(".gz", "")
                gunzip_file(retrieved_file, output_file)
                logger.info("\tUnarchived to {}".format(output_file))
            else:
                output_file = retrieved_file
            if not args.no_indexing:
                run_hmmpress(output_file)

if args.assembly:
    # print(args.assembly)
    out_dir = args.out or "."
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    logger.info("Trying to fetch flat file(s) from Assembly DB. Files will be written into '{}'".format(out_dir))
    for accession in args.assembly:
        logger.info("Downloading a flat file for {}...".format(accession))
        retrieved_file = retrieve_assembly(accession, out_dir)
        if retrieved_file:
            output_genbank_file = retrieved_file.replace(".gbk.gz", ".gbk")
            gunzip_file(retrieved_file, output_genbank_file)
            logger.info("\tDownloaded to {}".format(os.path.abspath(output_genbank_file)))

if args.assembly_fasta:
    # print(args.assembly_fasta)
    out_dir = args.out or "."
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    logger.info("Trying to fetch FASTA file(s) from Assembly DB. Files will be written into '{}'".format(out_dir))
    for accession in args.assembly_fasta:
        logger.info("Downloading a FASTA file for {}...".format(accession))
        retrieved_file = retrieve_assembly_fasta(accession, out_dir)
        if retrieved_file:
            output_genbank_file = retrieved_file.replace(".fna.gz", ".fna")
            gunzip_file(retrieved_file, output_genbank_file)
            logger.info("\tDownloaded to {}".format(os.path.abspath(output_genbank_file)))

if args.plasmidfinder:
    db_root = get_db_root(args)
    out_dir = db_root  # Reference data will be downloaded to DB_ROOT/plasmidfinder_db
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    logger.info("Trying to retrieve PlasmidFinder reference data. Files will be written into '{}'".format(out_dir))
    retrieve_plasmidfinder_reference(out_dir)

if args.mefinder:
    db_root = get_db_root(args)
    out_dir = db_root  # BLAST index will be built into DB_ROOT/mefinder_db
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    logger.info("Setting up MobileElementFinder. Index will be written into '{}'".format(out_dir))
    retrieve_mefinder_reference(out_dir, no_indexing=args.no_indexing)


# test GCF_000091005.1 GCA_000008865.1
    # accession = sys.argv[1]
    # retrieved_file = retrieve_assembly(accession)
    # output_genbank_file = retrieved_file.replace(".gbk.gz", ".gbk")
    # gunzip_file(retrieved_file, output_genbank_file)
