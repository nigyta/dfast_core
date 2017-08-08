#! /usr/bin/env python
import sys
import os
# import platform
from logging import getLogger, DEBUG, INFO, StreamHandler

app_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
logger = getLogger("")
logger.setLevel(INFO)

# srcPath = os.path.join(os.path.dirname(__file__), "..", "dfc")
sys.path.append(app_root)
from dfc.utils.path_util import set_binaries_path
from dfc.tools.ghostx import Ghostx
from dfc.tools.ghostz import Ghostz
from dfc.tools.blastp import Blastp
from dfc.tools.hmmer import Hmmer_hmmpress
from dfc.models.protein import Protein

set_binaries_path(app_root)



dbPath = os.path.join(os.path.dirname(__file__), "..", "db")

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
        fh.write("#id\tdescription\tgene\tEC_number\tflag\torganism\tsource_DB\tsequence\n")
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
        fh = sys.stdout
    else:
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
        fh.write("#id\tdescription\tgene\tEC_number\tflag\torganism\tsource_DB\tsequence\n")
    for protein in D.values():
        infer_source_db = False
        fh.write(protein.to_tsv(infer_source_db=infer_source_db))
    fh.close()

def make_blast_db(file_name, db_name=None):
    if db_name is None:
        db_name = file_name
    blastp = Blastp()
    logger.info("Preparing BLAST database. {0}".format(db_name))
    logger.setLevel(DEBUG)

    blastp.format_db(file_name, db_name)
    logger.setLevel(INFO)
    logger.info("Done")


def make_ghost_db(file_name, db_name=None):
    if db_name is None:
        db_name = file_name
    ghostx = Ghostx()
    logger.info("Preparing GHOSTX/Z database. {0}".format(db_name))
    logger.setLevel(DEBUG)
    ghostx.format_db(file_name, db_name)
    logger.setLevel(INFO)
    logger.info("Done")


def prepare_database(file_name):
    base_name, ext = (os.path.splitext(file_name))
    output_file = base_name + ".faa"
    logger.info("Converting DFAST reference '{0}' to FASTA '{1}'".format(file_name, output_file))
    fasta_file = dfast2fasta(file_name, output_file)
    make_ghost_db(fasta_file, base_name)
    make_blast_db(fasta_file, base_name)


def run_hmmpress(file_name):
    hmmpress = Hmmer_hmmpress()
    logger.info("Preparing HMM database using hmmpress. [{0}]".format(file_name))
    logger.setLevel(DEBUG)
    hmmpress.run(file_name)
    logger.setLevel(INFO)
    logger.info("Done")

if __name__ == '__main__':

    def f2d(args):
        if args.input is None:
            parser_toDfast.print_help()
            exit()
        fasta2dfast(args.input, args.output, args.no_header, args.source_db, args.attributes)

    def g2d(args):
        if args.input is None:
            parser_gbk2dfast.print_help()
            exit()
        genbank2dfast(args.input, args.output, args.no_header, args.attributes)

    def d2f(args):
        if args.input is None:
            parser_toFasta.print_help()
            exit()
        dfast2fasta(args.input, args.output, args.with_description)

    def formatdb(args):
        if args.input is None:
            parser_format_db.print_help()
            exit()
        prepare_database(args.input)

    def formathmm(args):
        if args.input is None:
            parser_format_hmm.print_help()
            exit()
        run_hmmpress(args.input)

    import argparse

    parser = argparse.ArgumentParser(description='DFAST reference utility script',
                                     usage=None, epilog=None, add_help=True)
    subparsers = parser.add_subparsers(help='')

    common_parser = argparse.ArgumentParser(add_help=False)

    common_parser.add_argument("-i", "--input", help="Input file", metavar="PATH")
    common_parser.add_argument("-o", "--output", help="Output file. If not specified, output is written to stdout.", metavar="PATH")

    parser_toDfast = subparsers.add_parser('fasta2dfast', help='Convert a FASTA file to a DFAST reference file.', parents=[common_parser])
    parser_toDfast.add_argument("--no_header", help="Do not add a header to output", action="store_true")
    parser_toDfast.add_argument("--source_db", help="Source database name (eg UniProtKB, INSD, RefSeq). If set to 'auto', try to infer from sequence ID.", metavar="STR")
    parser_toDfast.add_argument("--attributes", help="DB attributes. eg '[dbname=DFAST-default] [version=1.0] [contributor=yt]'", metavar="STR")
    parser_toDfast.set_defaults(func=f2d)

    parser_gbk2dfast = subparsers.add_parser('gbk2dfast', help='Convert a FASTA file to a DFAST reference file.', parents=[common_parser])
    parser_gbk2dfast.add_argument("--no_header", help="Do not add a header to output", action="store_true")
    parser_gbk2dfast.add_argument("--attributes", help="DB attributes. eg '[dbname=DFAST-default] [version=1.0] [contributor=yt]'", metavar="STR")
    parser_gbk2dfast.set_defaults(func=g2d)

    parser_toFasta = subparsers.add_parser('dfast2fasta', help='Convert a DFAST reference file to a FASTA file.', parents=[common_parser])
    parser_toFasta.add_argument("--with_description", help="Include description in FASTA headers.", action="store_true")
    parser_toFasta.set_defaults(func=d2f)

    parser_format_db = subparsers.add_parser('formatdb', help='Create BLAST and GHOSTX/Z databases from a DFAST reference file.')
    parser_format_db.add_argument("-i", "--input", help="Input reference file", metavar="PATH", required=True)
    parser_format_db.set_defaults(func=formatdb)

    parser_format_hmm = subparsers.add_parser('formathmm', help='Create profile-HMM database using hmmpress.')
    parser_format_hmm.add_argument("-i", "--input", help="Input HMM library", metavar="PATH", required=True)
    parser_format_hmm.set_defaults(func=formathmm)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        exit()

    # handler.setLevel(INFO)

    # file_name = os.path.join(dbPath, "DFAST_RefSeq.ref")
    # file_name = sys.argv[1]
    # dfast2fasta(file_name)
    # fasta2dfast(file_name)


    handler = StreamHandler()
    logger.addHandler(handler)

    args.func(args)
    # prepare_database(file_name)
    # run_hmmpress(file_name)
