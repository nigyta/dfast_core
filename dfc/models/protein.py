#! /usr/bin/env python
# coding: UTF8

from logging import getLogger
from Bio import SeqIO
from ..utils.ref_util import fasta_parsers, get_source_db

class Protein(object):
    def __init__(self, id_, description, gene, ec_number, flag, organism, source_db, sequence):
        self.id = id_
        self.description = description
        self.gene = gene
        self.ec_number = ec_number
        self.flag = flag
        self.organism = organism
        self.source_db = source_db
        self.sequence = sequence

    def __str__(self):
        return "<Protein:{self.id} {self.description}>".format(self=self)

    def to_fasta(self, with_description=False):
        if with_description:
            return ">{self.id} {self.description}\n{self.sequence}\n".format(self=self)
        else:
            return ">{self.id}\n{self.sequence}\n".format(self=self)

    def to_tsv(self, infer_source_db=False):
        if infer_source_db:
            self.source_db = get_source_db(self.id)
        return "\t".join([self.id, self.description, self.gene, self.ec_number,
                          self.flag, self.organism, self.source_db, self.sequence]) + "\n"

    @staticmethod
    def read_from_dfast_reference(ref_file_name):

        def _parse_header(header):
            assert header.startswith("#")
            attributes = header.strip("\n# ").split("\t")
            attributes = [_.split(":") for _ in attributes if ":" in _]
            return {"_" + key: value for key, value in attributes}

        D = {}
        with open(ref_file_name) as f:
            # first line: db description, second: field names, third-: data
            data = f.readlines()
            D.update(_parse_header(data[0]))
            assert data[1].startswith("#")
            data = data[2:]
            for row in data:
                row = row.strip("\n").split("\t")
                protein = Protein(*row)
                D[protein.id] = protein
        return D

    @staticmethod
    def read_from_genbank(gbk_file_name):
        logger = getLogger(__name__)
        D = {}
        suppress_warning = False
        for record in SeqIO.parse(open(gbk_file_name), "genbank"):
            organism = record.annotations.get("organism", "")
            source_db = "RefSeq" if "RefSeq" in record.annotations["keywords"] else "INSD"

            for feature in record.features:
                if feature.type == "CDS":
                    protein_id = feature.qualifiers.get("protein_id", [""])[0]
                    description = feature.qualifiers.get("product", [""])[0]
                    #         if "partial" in description:
                    #             print(description)
                    sequence = feature.qualifiers.get("translation", [""])[0]
                    gene = feature.qualifiers.get("gene", [""])[0]
                    ec_number = feature.qualifiers.get("EC_number", [""])[0]
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    if not protein_id:
                        if locus_tag and sequence:
                            if not suppress_warning:
                                logger.warning("Some of the CDS features have no protein_id.\n" +
                                               "Locus_tag will be used as an identifier. " + ""
                                               "This may not be compliant to the INSDC.")
                                suppress_warning = True
                            protein = Protein(locus_tag, description, gene, ec_number, "", organism, "", sequence)
                            D[protein.id] = protein
                    else:
                        protein = Protein(protein_id, description, gene, ec_number, "", organism, source_db, sequence)
                        D[protein.id] = protein
        # logger.info("{0} protein sequences were loaded from the genbank file.".format(len(D)))
        return D

    @staticmethod
    def read_from_fasta(fasta_file_name, parser_type="auto"):
        # logger = getLogger(__name__)
        parser = fasta_parsers[parser_type]
        D = {}
        for record in SeqIO.parse(open(fasta_file_name), "fasta"):
            s_id, product, organism, gene, source_db, ec_number = parser(record.id, record.description)
            sequence = str(record.seq)
            flag = ""
            protein = Protein(s_id, product, gene, ec_number, flag, organism, source_db, sequence)
            D[record.id] = protein
        return D

    @staticmethod
    def read_reference(file_name):
        """Read protein references from file. File format will be automatically inferred."""
        with open(file_name) as f:
            first_char = f.read(1)
        if first_char == ">":  # fasta format
            return Protein.read_from_fasta(file_name, parser_type="auto")
        elif first_char == "#":  # DFAST reference format
            return Protein.read_from_dfast_reference(file_name)
        elif first_char == "L":  # GenBank format (starts with "LOCUS")
            return Protein.read_from_genbank(file_name)
        else:
            logger = getLogger(__name__)
            logger.error("Unkown reference format. {0}".format(file_name))
            raise AssertionError

    @staticmethod
    def write_as_fasta(D, file_name, with_description=False):
        """
        Output protein fasta file from dictionary.

        :param D: A dictionary containing Protein objects, which is created by other methods (e.g. read_from_genbank)
        :param file_name: Output file name
        :return:
        """
        buffer = ""
        for protein in D.values():
            assert isinstance(protein, Protein)
            buffer += protein.to_fasta(with_description=with_description)
        with open(file_name, "w") as f:
            f.write(buffer)


if __name__ == '__main__':
    pass
