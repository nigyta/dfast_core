#!/usr/bin/env python3


# convert dfast gff for roary-compliant gff
# usage: python dfast2roary.py -i genome.gff > genome.mod.gff
# or, to convert multiple gff files (with an extension of .gff or .gff3) to roary style gff
# python dfast2roary.py -I input_directory -O output_directory
#
# Caveat: Do not need to run this script when DFAST is envoked using --use_locustag_as_gene_id option

import os
import sys
import argparse
import urllib.parse
import glob


def decode_attrs(file_name, out_file=None):
    if out_file:
        fh = open(out_file, "w")
    else:
        fh = sys.stdout

    in_f = open(file_name)
    for line in in_f:
        if line.startswith("##FASTA"):
            fh.write(line)
            break
        elif line.startswith("#"):
            fh.write(line)
        else:
            cols = line.strip("\n").split("\t")
            attrs = cols [8]
            attrs = urllib.parse.unquote(attrs)
            cols[8] = attrs
            ret = "\t".join(cols) + "\n"
            fh.write(ret)
            
    for line in in_f:
        fh.write(line)

def encode_attrs(file_name, out_file=None):
    if out_file:
        fh = open(out_file, "w")
    else:
        fh = sys.stdout

    in_f = open(file_name)
    for line in in_f:
        if line.startswith("##FASTA"):
            fh.write(line)
            break
        elif line.startswith("#"):
            fh.write(line)
        else:
            cols = line.strip("\n").split("\t")
            attrs = cols [8]
            attrs = attrs.split(";")
            mod_attrs = []
            for attr in attrs:
                if "=" not in attr:
                    sys.stderr.write("[Warning] No '=' in the attritbute '{}'\n".format(attr))
                    mod_attrs.append(attr)
                else:
                    if attr.count("=") > 1:
                        sys.stderr.write("[Warning] More than one '=' are found in the attritbute '{}'. This may cause an error.\n".format(attr))

                    key, value = attr.split("=", 1)
                    value = urllib.parse.quote(value)
                    if key == "ID":
                        continue  # Meant to replace ID with locus_tag
                    mod_attrs.append(key + "=" + value)
                    if key == "locus_tag":
                        mod_attrs.append("ID=" + value)  # Meant to replace ID with locus_tag


            cols[8] = ";".join(mod_attrs)
            ret = "\t".join(cols) + "\n"
            fh.write(ret)
            
    for line in in_f:
        fh.write(line)

def exec_directory(input_directory, output_directory):
    target_path = os.path.join(input_directory, "*.gff*")
    file_list = glob.glob(target_path)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    for input_file in file_list:
        base, ext = os.path.splitext(input_file)
        output_file = os.path.join(output_directory, base + ".mod" + ext)
        sys.stderr.write("Converting {} to {}\n".format(input_file, output_file))
        encode_attrs(input_file, output_file)

parser = argparse.ArgumentParser(description='Convert DFAST.gff to roary.gff by encoding GFF attributes and modifying ID.',
                                     usage=None, epilog=None, add_help=True)
group = parser.add_mutually_exclusive_group()
group.add_argument("-i", "--input", help="Input GFF file. Output will be written to STDOUT.", metavar="PATH")
group.add_argument("-I", "--input_directory", help="Input directory name", metavar="PATH")
parser.add_argument("-O", "--output_directory", help="Output directory name. Required if -I is specified.", metavar="PATH")


if __name__ == '__main__':
    args = parser.parse_args()
    if args.input:
        if args.output_directory:
            sys.stderr.write("-i and -O cannot be used together.\n")
        else:
            encode_attrs(args.input, out_file=None)
    elif args.input_directory:
        if not args.output_directory:
            sys.stderr.write("-O must be specified when -I is used.\n")
        else:
            exec_directory(args.input_directory, args.output_directory)
    else:
        parser.print_help()
