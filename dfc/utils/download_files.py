#! /usr/bin/env python
# coding: UTF8

import sys
import os
import gzip
from logging import getLogger, StreamHandler, INFO
from ftplib import FTP
import tarfile

version = sys.version_info.major
if version == 3:
    from urllib import request
else:
    from six.moves.urllib import request

logger = getLogger(__name__)

app_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", ".."))
db_root = os.path.join(app_root, "db")
ncbi_ftp_server = "ftp.ncbi.nlm.nih.gov"
cdd_directory = "/pub/mmdb/cdd/little_endian/"
host_dfast = "https://dfast.nig.ac.jp"

db_urls = {
    "dfast": host_dfast + "/dfc/distribution/DFAST-default.ref.gz",
    "bifido": host_dfast + "/dfc/distribution/DFAST-BIFIDO.ref.gz",
    "cyanobase": host_dfast + "/dfc/distribution/cyanobase_aa.ref.gz",
    "ecoli": host_dfast + "/dfc/distribution/DFAST-ECOLI.ref.gz",
    "lab": host_dfast + "/dfc/distribution/DFAST-LAB.ref.gz",
}

hmm_urls = {
    "Pfam": "ftp://ftp.ebi.ac.uk//pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz",
    "dbCAN": "http://csbl.bmb.uga.edu/dbCAN/download/dbCAN-fam-HMMs.txt",
    "TIGR": "ftp://ftp.jcvi.org//pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.LIB.gz"
}

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


def retrieve_cdd(db_name, out_dir="."):
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



if __name__ == '__main__':
    retrieve_dfast_reference("dfast")
