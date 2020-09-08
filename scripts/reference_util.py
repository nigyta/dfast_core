#! /usr/bin/env python
import sys
import os
# import platform
from logging import getLogger, DEBUG, INFO, StreamHandler

app_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
logger = getLogger("")
logger.setLevel(INFO)

# srcPath = os.path.join(os.path.dirname(__file__), "..", "dfc")
sys.path.append(app_root)
from dfc.utils.path_util import set_binaries_path
from dfc.tools.ghostx import Ghostx
from dfc.tools.ghostz import Ghostz
from dfc.tools.diamond import Diamond
from dfc.tools.blastp import Blastp
from dfc.tools.hmmer import Hmmer_hmmpress
from dfc.models.protein import Protein
from dfc.utils.reffile_util import fasta2dfast, dfast2fasta, genbank2dfast, prepare_database, prepare_database_dmnd, run_hmmpress
set_binaries_path(app_root)



dbPath = os.path.join(os.path.dirname(__file__), "..", "db")


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

    def format_dmnddb(args):
        if args.input is None:
            parser_format_dmnddb.print_help()
            exit()
        prepare_database_dmnd(args.input)

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

    parser_format_db = subparsers.add_parser('formatdb', help='Create BLAST and GHOSTX databases from a DFAST reference file.')
    parser_format_db.add_argument("-i", "--input", help="Input reference file", metavar="PATH", required=True)
    parser_format_db.set_defaults(func=formatdb)

    parser_format_hmm = subparsers.add_parser('formathmm', help='Create profile-HMM database using hmmpress.')
    parser_format_hmm.add_argument("-i", "--input", help="Input HMM library", metavar="PATH", required=True)
    parser_format_hmm.set_defaults(func=formathmm)

    parser_format_dmnddb = subparsers.add_parser('formatdb-dmnd', help='[Beta] Create a database file for Diamond from a DFAST reference file.')
    parser_format_dmnddb.add_argument("-i", "--input", help="Input reference file", metavar="PATH", required=True)
    parser_format_dmnddb.set_defaults(func=format_dmnddb)

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
