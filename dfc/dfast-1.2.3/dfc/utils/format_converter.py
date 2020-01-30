#! /usr/bin/env python
# coding: UTF8

import os
from logging import getLogger
from copy import deepcopy
from Bio import SeqIO

logger = getLogger(__name__)


def write_results(genome, config):
    fc = FormatConverter(genome, config)
    fc.write()


class FormatConverter(object):

    def __init__(self, genome, config):
        self.genome = genome
        self.work_dir = genome.workDir
        self.verbosity = config.OUTPUT_RESULT.get("verbosity", 3)
        self.logger = getLogger(__name__)

    def write(self):
        genome = deepcopy(self.genome)
        self.logger.info("Setting output verbosity level to {}.".format(self.verbosity))
        for record in genome.seq_records.values():
            for feature in record.features:
                feature.assign_hit(verbosity=self.verbosity)

        genome_fasta = os.path.join(self.work_dir, "genome.fna")
        write_genome_fasta(genome, genome_fasta)

        cds_fasta = os.path.join(self.work_dir, "cds.fna")
        write_cds_fasta(genome, cds_fasta)

        aa_fasta = os.path.join(self.work_dir, "protein.faa")
        write_aa_fasta(genome, aa_fasta)

        rna_fasta = os.path.join(self.work_dir, "rna.fna")
        write_rna_fasta(genome, rna_fasta)

        gff_file = os.path.join(self.work_dir, "genome.gff")
        write_gff(genome, gff_file)

        gbk_file = os.path.join(self.work_dir, "genome.gbk")
        write_genbank(genome, gbk_file)

        embl_file = os.path.join(self.work_dir, "genome.embl")
        write_embl(genome, embl_file)


def write_genbank(genome, file_name):
    logger.info("Writing a GenBank format file to {}".format(file_name))
    with open(file_name, "w") as f:
        SeqIO.write(list(genome.seq_records.values()), f, "genbank")


def write_embl(genome, file_name):
    logger.info("Writing an EMBL format file to {}".format(file_name))
    with open(file_name, "w") as f:
        SeqIO.write(list(genome.seq_records.values()), f, "embl")


def write_gff(genome, output_file):

    def _get_phase(feature):
        return str(int(feature.qualifiers.get("codon_start", [1])[0]) - 1) if feature.type == "CDS" else "."

    def _get_source(feature):
        for val in feature.qualifiers.get("inference", []):
            if val.startswith("COORDINATE"):
                return ":".join(val.split(":")[2:])
        return "DFAST"

    def _quote(x):
        return str(x).replace("%", "%25").replace("=", "%3d").replace(";", "%3B")

    def _create_attribute(feature):
        ignored_key = ["codon_start", "transl_table"]
        attributes = ";".join([key + "=" + ",".join(map(_quote, value))
                               for key, value in feature.qualifiers.items() if key not in ignored_key])
        return "ID=" + feature.id + ";" + attributes

    gff_buffer = "##gff-version 3\n"
    for record in genome.seq_records.values():
        for feature in record.features:
            if feature.type not in ["gene", "source"]:
                seq_id = record.name  # col0
                inference_source = _get_source(feature)  # col1
                start_pos = str(feature.location.start + 1)  # col3
                end_pos = str(int(feature.location.end))  # col4
                score = "."  # col5
                strand = "-" if feature.location.strand == -1 else "+"  # col6
                phase = _get_phase(feature)
                attributes = _create_attribute(feature)
                line = "\t".join([seq_id, inference_source, feature.type, start_pos, end_pos, score, strand, phase,
                                  attributes]) + "\n"
                gff_buffer += line

    logger.info("Writing a GFF file to {}".format(output_file))
    gff_buffer += "##FASTA\n"
    for record in genome.seq_records.values():
        gff_buffer += ">{}\n{}\n".format(record.name, str(record.seq))
    with open(output_file, "w") as f:
        f.write(gff_buffer)


def write_aa_fasta(genome, output_file):
    logger.info("Writing a protein fasta file to {}".format(output_file))
    fasta_buffer = ""
    for record in genome.seq_records.values():
        for feature in record.features:
            if feature.type == "CDS":
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                product = feature.qualifiers.get("product", [None])[0]
                translation = feature.qualifiers.get("translation", [None])[0]
                header = ">" + feature.id
                if locus_tag:
                    header += "|" + locus_tag
                if product:
                    header += " " + product
                if translation:
                    fasta_buffer += "{}\n{}\n".format(header, translation)
    with open(output_file, "w") as f:
        f.write(fasta_buffer)


def write_cds_fasta(genome, output_file):
    logger.info("Writing a CDS fasta file to {}".format(output_file))
    fasta_buffer = ""
    for record in genome.seq_records.values():
        for feature in record.features:
            if feature.type == "CDS":
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                product = feature.qualifiers.get("product", [None])[0]
                nucleotide = feature.location.extract(record.seq)
                header = ">" + feature.id
                if locus_tag:
                    header += "|" + locus_tag
                if product:
                    header += " " + product
                fasta_buffer += "{}\n{}\n".format(header, nucleotide)
    with open(output_file, "w") as f:
        f.write(fasta_buffer)


def write_rna_fasta(genome, output_file):
    logger.info("Writing an RNA fasta file to {}".format(output_file))
    fasta_buffer = ""
    for record in genome.seq_records.values():
        for feature in record.features:
            if feature.type.endswith("RNA"):
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                product = feature.qualifiers.get("product", [None])[0]
                nucleotide = feature.location.extract(record.seq)
                header = ">" + feature.id
                if locus_tag:
                    header += "|" + locus_tag
                if product:
                    header += " " + product
                fasta_buffer += "{}\n{}\n".format(header, nucleotide)
    with open(output_file, "w") as f:
        f.write(fasta_buffer)


def write_genome_fasta(genome, output_file):
    logger.info("Writing a genome fasta file to {}".format(output_file))
    fasta_buffer = ""
    for record in genome.seq_records.values():
        seq = str(record.seq)
        fasta_buffer += ">{}\n{}\n".format(record.name, seq)
    with open(output_file, "w") as f:
        f.write(fasta_buffer)

if __name__ == '__main__':
    pass
