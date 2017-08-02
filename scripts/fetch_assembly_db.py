#! /usr/bin/env python
# coding: UTF8

import os
import gzip
from ftplib import FTP
from six.moves import urllib

hostname = "ftp.ncbi.nlm.nih.gov"


def get_ftp_directory(accession):
    path1, path2, path3, path4 = accession[0:3], accession[4:7], accession[7:10], accession[10:13]
    return "/".join(["/genomes", "all", path1, path2, path3, path4])


def retrieve_assembly(accession, out_dir="."):
    directory = get_ftp_directory(accession)
    ftp = FTP(host=hostname)
    ftp.login()
    ftp.cwd(directory)
    L = ftp.nlst()
    print(L)
    asm_name = sorted([x for x in L if x.startswith(accession)])[-1]
    target_file = "/".join([directory, asm_name, asm_name + "_genomic.gbff.gz"])
    output_file = os.path.join(out_dir, "_".join(asm_name.split("_")[0:2]) + ".gbk.gz")
    with open(output_file, "wb") as f:
        ftp.retrbinary("RETR " + target_file, f.write)
    ftp.quit()
    return output_file


def gunzip_file(input_file, output_file):
    with open(output_file, "w") as fw:
        with gzip.open(input_file, 'r') as fr:
            fw.write(fr.read().decode("UTF-8"))




if __name__ == '__main__':
    import sys
    accession = sys.argv[1]
    retrieved_file = retrieve_assembly(accession)
    output_genbank_file = retrieved_file.replace(".gbk.gz", ".gbk")
    gunzip_file(retrieved_file, output_genbank_file)
