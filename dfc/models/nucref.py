#! /usr/bin/env python
# coding: UTF8

from logging import getLogger
from Bio import SeqIO
from ..utils.ref_util import fasta_parsers, get_source_db


class NucRefBase():

    def set_info_to_feature(self, feature):
        """
        to be implemented in the subclass
        set product, gene, etc. to feature obj
        """
        pass

    @staticmethod
    def write_tsv_file(dict_ref, out_tsv_file):
        with open(out_tsv_file, "w") as f:
            for ref in dict_ref.values():
                f.write(ref.to_tabular() + "\n")

    @classmethod
    def load_from_tsv_file(cls, tsv_file):
        D = {}
        for line in open(tsv_file):
            id_, obj = cls.parse_line(line)
            D[id_] = obj
        return D

    @staticmethod
    def write_nucl_fasta_file(dict_ref, out_fasta_file):
        with open(out_fasta_file, "w") as f:
            for ref in dict_ref.values():
                f.write(ref.to_nucl_fasta())

    @staticmethod
    def write_prot_fasta_file(dict_ref, out_fasta_file):
        with open(out_fasta_file, "w") as f:
            for ref in dict_ref.values():
                f.write(ref.to_prot_fasta())

class NucRef(NucRefBase):
    def __init__(self, protein_id, gene, product, gene_synonym, prot_seq, nucl_seq, note, accession, plasmid_name):
        self.protein_id = protein_id
        self.gene = gene
        self.product = product
        self.gene_synonym = gene_synonym
        self.prot_seq = prot_seq
        self.nucl_seq = nucl_seq
        self.note = note
        self.accession = accession
        self.plasmid_name = plasmid_name

    def __str__(self):
        return "<NucRef:{self.protein_id} {self.product}>".format(self=self)

    def set_info_to_feature(self, feature):
        feature.qualifiers["product"] = [self.product]
        if self.gene:
            feature.qualifiers["gene"] = [self.gene]
        if self.gene_synonym:
            feature.qualifiers["gene_synonym"] = [self.gene_synonym]
        note = f"similar to {self.gene} in Plasmid:{self.plasmid_name}"
        feature.qualifiers.setdefault("note", []).append(note)

    def to_tabular(self):         
        return "\t".join([self.protein_id, self.gene, self.product, self.gene_synonym, self.prot_seq, self.nucl_seq, self.note, self.accession, self.plasmid_name])

    def to_nucl_fasta(self):
        return f">{self.protein_id}\n{self.nucl_seq}\n"

    def to_prot_fasta(self):
        if self.prot_seq:
            return f">{self.protein_id}\n{self.prot_seq}\n"
        else:
            return ""

    @staticmethod
    def parse_line(line):

        cols = line.strip("\n").split("\t")
        nucref = NucRef(*cols)
        return nucref.protein_id, nucref

    @staticmethod
    def read_ref_file(ref_file_name):
        # a method for parsing the spreadsheet used for manual curation
        # for dev purpose only
        D = {}
        with open(ref_file_name) as f:
            _ = next(f)
            for line in f:
                cols = line.strip("\n").split("\t")
                accession, plasmid_name, protein_id, gene, product, gene_synonym, prot_seq, nucl_seq, note = cols
                prot_seq = prot_seq.strip("-*").upper()
                nucl_seq = nucl_seq.upper()

                nucref = NucRef(protein_id, gene, product, gene_synonym, prot_seq, nucl_seq, note, accession, plasmid_name)
                if prot_seq:
                    D[nucref.protein_id] = nucref
                else:
                    print(f"Skipping {nucref}")
        return D

    def to_dict(self):
        gene = self.gene
        if self.gene_synonym:
            gene += f" (synonym: {self.gene_synonym})"
        note = f"similar to {self.gene} in Plasmid:{self.plasmid_name}"
        if self.note:
            note += f", Note: {self.note}"
        return {"accession": f"PLADMID_DB:{self.protein_id}",
                "gene": gene,
                "product": self.product,
                "note": note
        }

class PLASMID_DB(NucRef):

    def __str__(self):
        return "<PLASMID_DB:{self.protein_id} {self.product}>".format(self=self)

    def to_dict(self):
        print("OVERRIDE")
        gene = self.gene
        if self.gene_synonym:
            gene += f" (synonym: {self.gene_synonym})"
        note = f"similar to {self.gene} in Plasmid:{self.plasmid_name}"
        if self.note:
            note += f", Note: {self.note}"
        return {"accession": f"PLADMID_DB:{self.protein_id}",
                "gene": gene,
                "product": self.product,
                "note": note
        }
