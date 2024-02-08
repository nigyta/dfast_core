#! /usr/bin/env python
# coding: UTF8

from logging import getLogger
from Bio import SeqIO
from ..utils.ref_util import fasta_parsers, get_source_db
from .nucref import NucRefBase


class CARD(NucRefBase):
    def __init__(self, aro_accession, ARO_name, gene_family, prot_accession, prot_seq, nucl_accession, nucl_seq, organism, note):

        self.aro_accession = aro_accession
        self.ARO_name = ARO_name
        self.gene_family = gene_family
        self.prot_accession = prot_accession
        self.prot_seq = prot_seq
        self.nucl_accession = nucl_accession
        self.nucl_seq = nucl_seq
        self.organism = organism
        self.note = note


    def __str__(self):
        return "<CARD:{self.aro_accession} {self.gene_family}>".format(self=self)

    def set_info_to_feature(self, feature):
        feature.qualifiers["product"] = [self.gene_family]
        # if self.gene:
        #     feature.qualifiers["gene"] = [self.gene]
        # if self.gene_synonym:
        #     feature.qualifiers["gene_synonym"] = [self.gene_synonym]
        note = f"similar to {self.ARO_name} in CARD:{self.aro_accession}"
        feature.qualifiers.setdefault("note", []).append(note)


    def to_tabular(self):         
        return "\t".join([self.aro_accession, self.ARO_name, self.gene_family, self.prot_accession, self.prot_seq, self.nucl_accession, self.nucl_seq, self.organism, self.note])

    def to_nucl_fasta(self):
        return f">{self.aro_accession}\n{self.nucl_seq}\n"

    def to_prot_fasta(self):
        if self.prot_seq:
            return f">{self.aro_accession}\n{self.prot_seq}\n"
        else:
            return ""


    @staticmethod
    def parse_line(line):

        aro_accession, ARO_name, gene_family, prot_accession, prot_seq, nucl_accession, nucl_seq, organism, note = line.strip("\n").split("\t")
        card = CARD(aro_accession, ARO_name, gene_family, prot_accession, prot_seq, nucl_accession, nucl_seq, organism, note)
        return card.aro_accession, card

    @staticmethod
    def read_ref_file(ref_file_name):
        D = {}
        with open(ref_file_name) as f:
            _ = next(f)
            for line in f:
                cols = line.strip("\n").split("\t")
                aro_accession, ARO_name, gene_family, prot_accession, prot_seq, nucl_accession, nucl_seq, organism, note = cols
                aa = aa.strip("-*").upper()
                nucl = nucl.upper()

                cds = CDS(protein_id, gene, product, gene_synonym, aa, nucl, note, accession, plasmid_name)
                if aa:
                    D[cds.protein_id] = cds
                else:
                    print(f"Skipping CDS for {cds}")
        return D

    @staticmethod
    def parse_card_json_file(card_json_file):
        def _get_gene_family_name(card_model):
            dict_category = card_model["ARO_category"]
            for aro_category in dict_category.values():
                if aro_category["category_aro_class_name"] == "AMR Gene Family":
                    return aro_category["category_aro_name"]
            else:
                return ""

        def _get_sequences(card_model):
            sequences = card_model["model_sequences"]["sequence"].values()
            assert len(sequences) == 1
            sequence = list(sequences)[0]
            protein_sequence = sequence["protein_sequence"]
            dna_sequence = sequence["dna_sequence"]
            NCBI_taxonomy = sequence["NCBI_taxonomy"]
            prot_accession, prot_seq = protein_sequence["accession"], protein_sequence["sequence"]
            nucl_accession, nucl_seq = dna_sequence["accession"], dna_sequence["sequence"]
            organism = NCBI_taxonomy["NCBI_taxonomy_name"]

            return prot_accession, prot_seq, nucl_accession, nucl_seq, organism

        def _parse_card_model(card_model):
            aro_accession = "ARO:" + card_model["ARO_accession"]
            ARO_name = card_model["ARO_name"]
            gene_family = _get_gene_family_name(card_model)
            prot_accession, prot_seq, nucl_accession, nucl_seq, organism = _get_sequences(card_model)
            note = ""  # currently not used.
            return CARD(aro_accession, ARO_name, gene_family, prot_accession, prot_seq, nucl_accession, nucl_seq, organism, note)


        D = json.load(open(card_json_file))
        ret = {}
        for key, card_model in D.items():
            if key.startswith("_"):
                print("skipping...", key)
                continue
            # print(card_model)
            if card_model["model_type"] == "protein homolog model":
                card = _parse_card_model(card_model)
                ret[card.aro_accession] = card
        return ret

    def to_dict(self):
        accession = f"{self.aro_accession}:{self.prot_accession}"
        note = f"similar to CARD:{self.ARO_name} in {self.organism} ({self.nucl_accession})"
        if self.note:
            note += f", {self.note}"
    
        return {"accession": accession,
                "gene": self.ARO_name,
                "product": self.gene_family,
                "note": note
        }