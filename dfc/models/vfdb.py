#! /usr/bin/env python
# coding: UTF8

from logging import getLogger
from Bio import SeqIO
from ..utils.ref_util import fasta_parsers, get_source_db
import re
from .nucref import NucRefBase


class VFDB(NucRefBase):
    def __init__(self, vfg, protein_id, product, gene_symbol, vf_id, vf_name, prot_seq, nucl_seq, organism, note):

        self.vfg = vfg
        self.protein_id = protein_id
        self.product = product
        self.gene_symbol = gene_symbol
        self.vf_id = vf_id
        self.vf_name = vf_name
        self.prot_seq = prot_seq
        self.nucl_seq = nucl_seq
        self.organism = organism
        self.note = note


    def __str__(self):
        return "<VFDB:{self.vfg} ({self.vf_name}) {self.product}".format(self=self)

    def info(self):
        return f"{self.vf_name} in VF_ID:{self.vf_id} ({self.organism})"

    def set_info_to_feature(self, feature):
        feature.qualifiers["product"] = [self.product]
        if self.gene_symbol:
            feature.qualifiers["gene"] = [self.gene_symbol]
        # if self.gene:
        #     feature.qualifiers["gene"] = [self.gene]
        # if self.gene_synonym:
        #     feature.qualifiers["gene_synonym"] = [self.gene_synonym]
        note = f"similar to {self.info()}"
        feature.qualifiers.setdefault("note", []).append(note)

    def to_tabular(self):         
        return "\t".join([self.vfg, self.protein_id, self.product, self.gene_symbol, self.vf_id, self.vf_name, self.prot_seq, self.nucl_seq, self.organism, self.note])

    def to_nucl_fasta(self):
        return f">{self.vfg}\n{self.nucl_seq}\n"

    def to_prot_fasta(self):
        if self.prot_seq:
            return f">{self.vfg}\n{self.prot_seq}\n"
        else:
            return ""

    @staticmethod
    def parse_line(line):
        vfg, protein_id, product, gene_symbol, vf_id, vf_name, prot_seq, nucl_seq, organism, note = line.strip("\n").split("\t")
        vfdb = VFDB(vfg, protein_id, product, gene_symbol, vf_id, vf_name, prot_seq, nucl_seq, organism, note)
        return vfdb.vfg, vfdb
    # @staticmethod
    # def read_ref_file(ref_file_name):
    #     D = {}
    #     with open(ref_file_name) as f:
    #         _ = next(f)
    #         for line in f:
    #             cols = line.strip("\n").split("\t")
    #             aro_accession, ARO_name, gene_family, prot_accession, prot_seq, nucl_accession, nucl_seq, organism, note = cols
    #             aa = aa.strip("-*").upper()
    #             nucl = nucl.upper()

    #             cds = CDS(protein_id, gene, product, gene_synonym, aa, nucl, note, accession, plasmid_name)
    #             if aa:
    #                 D[cds.protein_id] = cds
    #             else:
    #                 print(f"Skipping CDS for {cds}")
    #     return D

    @staticmethod
    def parse_vfdb_fasta_files(vfdb_nucl_fasta_file, vfdb_prot_fasta_file):
        pat_vfdb = re.compile(r"^(?P<vfg>VFG\d+)(\(gb\|(?P<prot_id>[\w_\.]+)\))? \((?P<symbol>.+)\) (?P<product>.+) \[(?P<vf_name>.+) \((?P<vf_id>VF\d+)\) - .+ \(VFC\d+\)\] \[(?P<organism>.+)\]$")
        ret = {}
        # for r_nucl, r_prot in zip(SeqIO.parse(vfdb_nucl_fasta_file, "fasta"), SeqIO.parse(vfdb_prot_fasta_file, "fasta")):

        R_prot = list(SeqIO.parse(vfdb_prot_fasta_file, "fasta"))
        dict_prot = {}
        for r in R_prot:
            r.id = r.id.split()[0].split("(")[0]
            dict_prot[r.id] = str(r.seq.upper())

        for r_nucl in SeqIO.parse(vfdb_nucl_fasta_file, "fasta"):
            descr_nucl = r_nucl.description
            # descr_prot = r_prot.description
            # assert r_nucl.id == r_prot.id

            m = pat_vfdb.match(descr_nucl)
            if m:
                vfg = m.group("vfg")
                prot_id = m.group("prot_id") or ""
                symbol = m.group("symbol")
                product = m.group("product")
                vf_name = m.group("vf_name")
                vf_id = m.group("vf_id")
                organism = m.group("organism")
                note = ""
                nucl_seq = str(r_nucl.seq.upper())
                prot_seq = dict_prot.get(vfg, "")  #str(r_prot.seq.upper())
                vfdb = VFDB(vfg, prot_id, product, symbol, vf_id, vf_name, prot_seq, nucl_seq, organism, note)
                ret[vfdb.vfg] = vfdb

            else:
                print("Regex does not match: " + descr_nucl)
                continue
        return ret

    def to_dict(self):
        accession = f"VFDB:{self.vfg}:{self.protein_id}"
        note = f"similar to {self.gene_symbol} in {self.organism}"
        if self.note:
            note += f", {self.note}"
    
        return {"accession": accession,
                "gene": self.gene_symbol,
                "product": self.product,
                "note": note
        }