#!/usr/bin/env python3


# Generate a gene summary table (features.tsv) from genome.gbk file

import os
import sys
import argparse
import csv
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Generate a gene summary table (features.tsv) from genome.gbk file.',
                                     usage=None, epilog=None, add_help=True)
parser.add_argument("-d", "--dir", help="Output directory name of DFAST results.", metavar="PATH", required=True)

def parseFeature(feature):
    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
    strand = "-" if feature.location.strand == -1 else "+"
    location = f"{feature.location.start}..{feature.location.end}({strand})"
    product = feature.qualifiers.get("product", [""])[0]
    if feature.type == "repeat_region":  # for CRISPR
        product = feature.qualifiers.get("rpt_family", [""])[0]
    gene = feature.qualifiers.get("gene", [""])[0]
    translation = feature.qualifiers.get("translation", [""])[0]
    return locus_tag, location, product, gene, translation

def createTSV(genbankFileName, tsvFileName):
    csvWriter = csv.writer(
        open(tsvFileName, "w"), lineterminator="\n", delimiter='\t')
    dat = []
    header = ("number", "locus", "sequence", "location", "feature",
               "product", "gene", "nucleotide", "translation")
    dat.append(header)
    num = 0
    for gbRecord in SeqIO.parse(genbankFileName, "genbank"):
        for feature in gbRecord.features:
            if feature.type not in ["source", "gene", "assembly_gap", "repeat_region"]:
                num += 1
                locus_tag, location, product, gene, translation = parseFeature(feature)

                nucSeq = feature.location.extract(gbRecord.seq)
                nucleotide = str(nucSeq)
                dat.append((str(num), locus_tag, gbRecord.name, location,
                               feature.type, product, gene, nucleotide, translation))
    csvWriter.writerows(dat)

if __name__ == '__main__':
    args = parser.parse_args()
    dfast_output_dir = args.dir
    gbk_file = os.path.join(dfast_output_dir, "genome.gbk")
    out_tsv_file = os.path.join(dfast_output_dir, "features.tsv")
    createTSV(gbk_file, out_tsv_file)
