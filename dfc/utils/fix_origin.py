#! /usr/bin/env python
# coding: UTF8

import os
import sys
import argparse
import tempfile
import shutil
from Bio import SeqIO
from logging import getLogger, DEBUG, INFO, StreamHandler, Formatter, FileHandler, CRITICAL, WARNING

app_root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
if app_root == "":
    app_root = "."

if __name__ == '__main__':
    sys.path.append(app_root)
    from dfc.utils.path_util import set_binaries_path
    set_binaries_path(app_root)

from dfc.genome import Genome
from dfc.functionalAnnotation import FunctionalAnnotation
from dfc.structuralAnnotation import StructuralAnnotation


class Config:
    GENOME_FASTA = None  # Will be overriden by -g or --genome option
    WORK_DIR = "RGU"
    CPU = 1
    FORCE_OVERWRITE = True
    DEBUG = True

    GENOME_CONFIG = {
        "complete": False,
        "use_original_name": True,
        "sort_by_length": False,
        "minimum_length": 0
    }

    GENOME_SOURCE_INFORMATION = {}

    STRUCTURAL_ANNOTATION = [
        {
            "tool_name": "MGA",
            "enabled": True,
            "options": {"cmd_options": "-s"},  # -s for single species, -m for multiple species
        },
    ]

    FUNCTIONAL_ANNOTATION = [
        {
            "component_name": "DnaAfinder",
            "enabled": True,
            "options": {
                "skipAnnotatedFeatures": False,
                "evalue_cutoff": 1e-20,
                "hmm_profile": os.path.join(app_root, "db/hmm_search//TIGR00362.hmm"),
                "offset": 0
            },
        },
    ]



def parse_arg():
    parser = argparse.ArgumentParser(description='A utility script to rotate the chromosome so that the dnaA gene comes first.\n',
                                     usage="fix_origin.py -g your_genome.fna [-o output_file] [--offset INT]", epilog=None, add_help=True)
    parser.add_argument("-g", "--genome", help="Genomic FASTA file", metavar="PATH", required=True)
    parser.add_argument("-o", "--out", help="Output file (default:stdout)", type=str, metavar="PATH")
    parser.add_argument("--offset", help="Offset from the start codon of the dnaA gene (default=0)",
                        default=0, type=int, metavar="INT")
    parser.add_argument("--debug", help="Debug mode", action="store_true")

    if len(sys.argv) == 1:
        parser.print_help()
        exit()

    return parser.parse_args()


def main():
    args = parse_arg()
    logger = getLogger(__name__)

    if args.debug:
        log_level = DEBUG
    else:
        log_level = WARNING

    handler = StreamHandler(sys.stderr)
    handler.setLevel(log_level)
    handler.setFormatter(Formatter("%(asctime)s %(message)s", datefmt='%Y/%m/%d %H:%M:%S'))
    logger.addHandler(handler)
    logger.setLevel(log_level)

    fix_origin(args.genome, output_file=args.out, offset=args.offset, logger=logger)


def fix_origin(input_file, output_file=None, offset=0, logger=None):
    if logger is None:
        logger = getLogger(__name__)
    if __name__ != '__main__':
        loglevel_bkup = logger.level
        logger.setLevel(WARNING)

    logger.warning("[WARNING] Trying to locate a dnaA gene to fix the sequence origin of the chromosome. DO NOT APPLY THIS TO A DRAFT GENOME.")
    config = Config()
    config.GENOME_FASTA = input_file
    config.FUNCTIONAL_ANNOTATION[0]["options"]["offset"] = offset
    tmp_dir = tempfile.mkdtemp()
    config.WORK_DIR = tmp_dir
    logger.info("A temporary working directory was created in {}".format(config.WORK_DIR))
    genome = Genome(config)

    sa = StructuralAnnotation(genome, config)
    fa = FunctionalAnnotation(genome, config)
    sa.execute()
    fa.execute()

    if output_file:
        f = open(output_file, "w")
        logger.info("The origin-fixed genome fasta file is written to {}.".format(output_file))
    else:
        f = sys.stdout
    SeqIO.write(list(genome.seq_records.values()), f, "fasta")

    shutil.rmtree(config.WORK_DIR)
    logger.info("The temporary working {} was cleaned up.".format(config.WORK_DIR))
    if __name__ != '__main__':
        logger.setLevel(loglevel_bkup)


if __name__ == '__main__':
  main()
