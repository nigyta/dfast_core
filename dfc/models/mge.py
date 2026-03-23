#! /usr/bin/env python
# coding: UTF8

from logging import getLogger
from Bio import SeqIO
from .nucref import NucRefBase


class MGE(NucRefBase):
    def __init__(self, mge_id, name, allele, accession, nucl_seq, note=""):
        self.mge_id = mge_id
        self.name = name
        self.allele = allele
        self.accession = accession
        self.nucl_seq = nucl_seq
        self.note = note

    def __str__(self):
        return f"<MGE:{self.mge_id} {self.name}>"

    def info(self):
        return f"{self.name} (allele {self.allele}, {self.accession})"

    def set_info_to_feature(self, feature):
        feature.qualifiers["product"] = [f"mobile element {self.name}"]
        note = f"similar to MGEdb:{self.info()}"
        feature.qualifiers.setdefault("note", []).append(note)

    def to_tabular(self):
        return "\t".join([self.mge_id, self.name, self.allele, self.accession, self.nucl_seq, self.note])

    def to_nucl_fasta(self):
        return f">{self.mge_id}\n{self.nucl_seq}\n"

    def to_prot_fasta(self):
        return ""

    @staticmethod
    def parse_line(line):
        cols = line.strip("\n").split("\t")
        mge_id, name, allele, accession, nucl_seq = cols[0], cols[1], cols[2], cols[3], cols[4]
        note = cols[5] if len(cols) > 5 else ""
        mge = MGE(mge_id, name, allele, accession, nucl_seq, note)
        return mge.mge_id, mge

    @staticmethod
    def parse_mge_fasta(fasta_file):
        """Parse MGEdb FASTA file. Headers are formatted as: >name|allele|accession"""
        D = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.description
            parts = header.split("|")
            if len(parts) >= 3:
                name = parts[0]
                allele = parts[1]
                accession = parts[2]
            else:
                name = header
                allele = ""
                accession = ""
            mge_id = header
            nucl_seq = str(record.seq).upper()
            mge = MGE(mge_id, name, allele, accession, nucl_seq)
            D[mge.mge_id] = mge
        return D

    def to_dict(self):
        note = f"similar to MGEdb:{self.name} ({self.accession})"
        if self.note:
            note += f", {self.note}"
        return {
            "accession": f"MGE:{self.mge_id}",
            "gene": self.name,
            "product": f"mobile element {self.name}",
            "note": note,
        }
