#! /usr/bin/env python
# coding: UTF8

import re


class Hit(object):

    def assign(self, feature, verbosity):
        raise NotImplementedError

    def assign_as_note(self, verbosity):
        raise NotImplementedError

    def get_inference(self):
        raise NotImplementedError


class ProteinHit(Hit):
    def __init__(self, id_, description, gene, ec_number, source_db, organism, db_name, e_value, score, identity, q_cov,
                 s_cov, flag, notes=None):
        self.id = id_
        self.description = description
        self.gene = gene
        self.ec_number = ec_number
        self.source_db = source_db
        self.organism = organism
        self.db_name = db_name
        self.score = score
        self.e_value = float(e_value)
        self.identity = float(identity)
        self.q_cov = q_cov
        self.s_cov = s_cov
        self.flag = flag
        if not notes:
            notes = []
        self.notes = notes

    def __str__(self):
        gene_id = self.id
        if self.db_name:
            gene_id = self.db_name + ":" + gene_id

        if self.organism:
            ret = "{gene_id} {self.description} ({self.organism}) [pid:{self.identity:.1f}%, q_cov:{self.q_cov:.1f}%, s_cov:{self.s_cov:.1f}%, Eval:{self.e_value:.1e}]".format(gene_id=gene_id, self=self)
        else:
            ret = "{gene_id} {self.description} [pid:{self.identity:.1f}%, q_cov:{self.q_cov:.1f}%, s_cov:{self.s_cov:.1f}%, Eval:{self.e_value:.1e}]".format(gene_id=gene_id, self=self)
        if self.flag:
            ret = ret[:-1] + ", " + self.flag + "]"
        return ret

    def assign(self, feature, verbosity=2):
        feature.qualifiers["product"] = [self.description]
        # print("DEBUG setting product", self.description)
        if self.gene:
            feature.qualifiers["gene"] = [self.gene]
        if self.ec_number:
            feature.qualifiers["EC_number"] = [x.strip() for x in self.ec_number.split(",")]
        feature.qualifiers.setdefault("inference", []).append(self.get_inference())
        # todo modify note qualifier into appropriate one
        feature.qualifiers.setdefault("note", []).extend(self.notes)
        if verbosity >= 2:
            feature.qualifiers.setdefault("note", []).append(str(self))

    def assign_as_note(self, feature, verbosity=2):
        if not self.flag:
            self.flag = "alternative hit"
        if verbosity >= 2:
            feature.qualifiers.setdefault("note", []).append(str(self))

    def get_inference(self):
        if self.source_db:
            return "similar to AA sequence:{0}:{1}".format(self.source_db, self.id)
            # return "DESCRIPTION:similar to AA sequence:{0}:{1}".format(self.source_db, self.id)
        else:
            return "similar to AA sequence:{0}".format(self.id)
            # return "DESCRIPTION:similar to AA sequence:{0}".format(self.id)


class HmmHit(Hit):

    def __init__(self, accession, name, description, evalue, score, bias, db_name):
        self.accession, self.name, self.description, self.evalue, self.score, self.bias, self.db_name = \
            accession, name, description, float(evalue), float(score), float(bias), db_name

    def __repr__(self):
        return "{hmm.db_name}:{hmm.accession}; {hmm.description} [Name:{hmm.name}, Eval:{hmm.evalue:.1e}, score:{hmm.score:.1f}, bias:{hmm.bias:.1f}]".format(
            hmm=self)

    def assign(self, feature, verbosity=2):
        self.assign_as_note(feature, verbosity=verbosity)

    def assign_as_note(self, feature, verbosity=2):
        if verbosity >= 2:
            feature.qualifiers.setdefault("note", []).append(str(self))

class CddHit(Hit):
    def __init__(self, result_type, hit_type, pssm_id, from_, to_, evalue, score, accession, short_name, incomplete, description):
        self.result_type, self.hit_type, self.pssm_id, self.from_, self.to_, self.evalue, self.score, self.accession, self.short_name, self.incomplete, self.description = \
            result_type, hit_type, pssm_id, from_, to_, float(evalue), float(score), accession, short_name, incomplete, description

        # modify result for COG hit. (Replacing description of the function to a category code).
        if self.accession.startswith("COG"):
            definition, category = get_cog_definition_and_category(self.description)
            self.result_type = "COG"
            self.description = definition
            self.category = category

    def __repr__(self):
        hit_info = "Category:" + self.category if self.result_type == "COG" else self.hit_type
        if self.accession == self.short_name:
            ret = "{cdd.result_type}:{cdd.accession} {cdd.description} [{hit_info}, Aligned:{cdd.from_}-{cdd.to_}, Eval:{cdd.evalue:.1e}, score:{cdd.score:.1f}".format(
                cdd=self, hit_info=hit_info)
        else:
            ret = "{cdd.result_type}:{cdd.accession}:{cdd.short_name} {cdd.description} [{hit_info}, Aligned:{cdd.from_}-{cdd.to_}, Eval:{cdd.evalue:.1e}, score:{cdd.score:.1f}".format(
                cdd=self, hit_info=hit_info)

        if "C" in self.incomplete:
            ret += ", C-term missing"
        elif "N" in self.incomplete:
            ret += ", N-term missing"
        ret += "]"
        return ret

    def assign(self, feature, verbosity=2):
        self.assign_as_note(feature, verbosity=verbosity)

    def assign_as_note(self, feature, verbosity=2):
        # if verbosity >= 2:
        #     feature.qualifiers.setdefault("note", []).append(str(self))

        if verbosity >= 2 and self.hit_type == "Multidom":
            feature.qualifiers.setdefault("note", []).append(str(self))

        elif verbosity >= 3:

            feature.qualifiers.setdefault("note", []).append(str(self))


cog_pat = re.compile(r"(.+) \[(.+)\]\.")
def get_cog_definition_and_category(description):
    m = cog_pat.match(description)
    if m:
        definition, category = m.group(1), m.group(2)
        for code, cat in cog_categories.items():
            category = category.replace(cat, code)
        category = category.replace(", ", "")
        return definition, category
    else:
        return description, "-"

cog_categories = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}


class PseudoGene(Hit):


    def __init__(self, stop_codon, insertion, deletion):
        """
        :param stop_codons: set of (stopcodon_start, stopcodon_end, strand)
        :param insertion: set of insertion location
        :param deletion: set of deletion location
        """
        self.stop_codon = stop_codon
        self.insertion = insertion
        self.deletion = deletion

    def __repr__(self):
        return "{hmm.db_name}:{hmm.accession}; {hmm.description} [Name:{hmm.name}, Eval:{hmm.evalue:.1e}, score:{hmm.score:.1f}, bias:{hmm.bias:.1f}]".format(
            hmm=self)

    def assign(self, feature, verbosity=2):
        self.assign_as_note(feature, verbosity=verbosity)

    def assign_as_note(self, feature, verbosity=2):
        if verbosity == 1:
            if len(self.stop_codon) > 0 or len(self.insertion) > 0 or len(self.deletion) > 0:
                note = "possible pseudo"
                if len(self.stop_codon) > 0:
                    note += ", internal stop codon"
                if len(self.insertion) > 0 or len(self.deletion) > 0:
                    note += ", frameshifted"
                feature.qualifiers.setdefault("note", []).append(note)

        elif verbosity >= 2:
            if len(self.stop_codon) > 0:
                stop_codons = {"[{0}:{1}]({2})".format(x[0] + 1, x[1], x[2]) for x in self.stop_codon}
                stop_codons = sorted(list(stop_codons))
                note = "internal stop codon at " + ",".join(stop_codons)
                feature.qualifiers.setdefault("note", []).append(note)
            if len(self.insertion) > 0 or len(self.deletion) > 0:
                note = "frameshifted"
                if len(self.insertion) > 0:
                    note += ", insertion at around " + ",".join(map(str, self.insertion))
                if len(self.deletion) > 0:
                    note += ", deletion at around " + ",".join(map(str, self.deletion))
                feature.qualifiers.setdefault("note", []).append(note)


class MBGDHit(Hit):
    def __init__(self, id_, mbgd_version, mbgd_tabid, clst_id, clst_descr, gene_symbol, gene_descr, e_value, score, identity, q_cov,
                 s_cov, flag, notes=None):
        self.id = id_
        self.mbgd_version = mbgd_version
        self.mbgd_tabid = mbgd_tabid
        self.clst_id = clst_id
        self.clst_descr = clst_descr
        self.gene_symbol = gene_symbol
        self.gene_descr = gene_descr
        self.score = score
        self.e_value = float(e_value)
        self.identity = float(identity)
        self.q_cov = q_cov
        self.s_cov = s_cov
        self.flag = flag
        if not notes:
            notes = []
        self.notes = notes

    def __str__(self):
        if self.clst_descr:
            return "MBGD: {{gene_id: '{x.id}', cluster_id: '{x.clst_id}', cluster_description: '{x.clst_descr}', gene: '{x.gene_symbol}', version: '{x.mbgd_version}', tabid: '{x.mbgd_tabid}', pid: '{x.identity:.1f}%', q_cov: '{x.q_cov:.1f}%', s_cov: '{x.s_cov:.1f}%', Eval: '{x.e_value:.1e}'}}".format(x=self)
        else:
            return "MBGD: {{gene_id: '{x.id}', cluster_id: '{x.clst_id}', gene_description: '{x.gene_descr}', version: '{x.mbgd_version}', tabid: '{x.mbgd_tabid}', pid: '{x.identity:.1f}%', q_cov: '{x.q_cov:.1f}%', s_cov: '{x.s_cov:.1f}%', Eval: '{x.e_value:.1e}'}}".format(x=self)

    def assign(self, feature, verbosity=2):
        self.assign_as_note(feature, verbosity)

    def assign_as_note(self, feature, verbosity=2):
        if verbosity >= 2:
            feature.qualifiers.setdefault("inference", []).append(self.get_inference())
            feature.qualifiers.setdefault("note", []).append(str(self))

    def get_inference(self):
        # inference: similar to AA sequence:MBGD:2016-01_default:317
        return "similar to AA sequence:MBGD:{0}_{1}:{2}".format(self.mbgd_version, self.mbgd_tabid, self.clst_id)


if __name__ == '__main__':
    pass
