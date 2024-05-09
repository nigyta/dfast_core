#! /usr/bin/env python
# coding: UTF8

# from logging import getLogger
# from Bio import SeqIO
# from ..utils.ref_util import fasta_parsers, get_source_db

from .nucref import NucRef


class PLASMID_DB(NucRef):

    def __str__(self):
        return "<PLASMID_DB:{self.protein_id} {self.product}>".format(self=self)

    def to_dict(self):
        gene = self.gene
        if self.gene_synonym:
            gene += f" (synonym: {self.gene_synonym})"
        note = f"similar to {self.gene} in Plasmid:{self.plasmid_name}"
        if self.note:
            note += f", Note: {self.note}"
        return {"accession": f"PLASMID_DB:{self.protein_id}",
                "gene": gene,
                "product": self.product,
                "note": note
        }
