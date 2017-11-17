#! /usr/bin/env python
# coding: UTF8

import os
from logging import getLogger
from copy import deepcopy
from .ddbj_submission import get_seq_rank

class GenBankSubmission(object):

    def __init__(self, genome, config):
        self.genome = genome
        self.output_dir = os.path.join(genome.workDir, "genbank")
        self.center_name = config.GENBANK_SUBMISSION.get("center_name", "")
        self.verbosity = config.GENBANK_SUBMISSION.get("output_verbosity", 1)
        self.enabled = config.GENBANK_SUBMISSION.get("enabled", True)
        self.metadata = {}
        self.logger = getLogger(__name__)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def get_file_prefix(self):
        bio_sample_id = self.metadata.get("biosample")
        strain = self.genome.strain
        if bio_sample_id and strain:
            return bio_sample_id + "_" + strain
        elif bio_sample_id:
            return bio_sample_id
        elif strain:
            return strain
        else:
            return "genbank"

    def create_submission_file(self):
        if self.enabled:
            if not self.center_name:
                self.logger.warning("'center_name' is not specified. 'my_center' is used tentatively.")
                self.center_name = "my_center"
            prefix = self.get_file_prefix()
            tbl_file = os.path.join(self.output_dir, prefix + ".tbl").replace(" ", "_").replace("(", "").replace(")", "")
            fsa_file = os.path.join(self.output_dir, prefix + ".fsa").replace(" ", "_").replace("(", "").replace(")", "")
            self.logger.info("Writing a GenBank tbl file to {}".format(tbl_file))
            self.logger.info("Writing a GenBank fsa file to {}".format(fsa_file))
            create_genbank_submission_file(self.genome, tbl_file, fsa_file, self.center_name, self.verbosity)
        else:
            self.logger.warning("'Generate GenBank Submission file' is disabled. Skip processing")


def location_to_str(location):
    start = str(int(location.start) + 1)
    end = str(int(location.end))
    if location.strand == 1:
        before = "<" if "<" in str(location.start) else ""
        after = ">" if ">" in str(location.end) else ""
        return before + start, after + end
    elif location.strand == -1:
        before = "<" if ">" in str(location.end) else ""
        after = ">" if "<" in str(location.start) else ""

        return before + end, after + start


def qualifier_to_tbl(qualifiers, key):
    ret = ""
    for value in qualifiers.get(key, []):
        ret += "\t\t\t{0}\t{1}\n".format(key, value)
    return ret


def feature_to_tbl(feature, center_name):
    ignored_qualifiers = ["locus_tag", "gene", "codon_start", "transl_table", "translation", "EC_number"]
    ret = ""
    start, end = location_to_str(feature.location)
    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
    if locus_tag:
        ret += "{0}\t{1}\t{2}\n".format(start, end, "gene")
        ret += "\t\t\t{0}\t{1}\n".format("locus_tag", locus_tag)
        ret += qualifier_to_tbl(feature.qualifiers, "gene")
    ret += "{0}\t{1}\t{2}\n".format(start, end, feature.type)
    for key in feature.qualifiers:
        if key not in ignored_qualifiers:
            ret += qualifier_to_tbl(feature.qualifiers, key)

    if feature.type == "CDS":
        if feature.qualifiers.get("codon_start", [1])[0] != 1:
            ret += qualifier_to_tbl(feature.qualifiers, "codon_start")
        ret += "\t\t\t{0}\tgnl|{1}|{2}\n".format("protein_id", center_name, locus_tag)
    return ret


def get_ff_definition(record, seq_rank):
    assert seq_rank in ["complete", "scaffold", "contig", "complete genome"]
    if seq_rank == "complete":
        seq_rank = "complete genome"
    organism = record.annotations.get("organism", "")
    if not organism:
        return ""
    strain = record.annotations.get("strain", "")
    seq_name = record.name
    ff_definition = "{organism} DNA, {seq_rank}: {seq_name}".format(organism=organism, seq_rank=seq_rank, seq_name=seq_name)
    if strain:
        ff_definition += ", strain: {strain}".format(strain=strain)
    return ff_definition


def create_fsa_modifiers(record):
    source_features = [f for f in record.features if f.type == "source"]
    if len(source_features) != 1:
        self.logger.error("Error. Source feature is empty or more than one source features are assigned.")
        exit(1)
    source_feature = source_features[0]
    modifiers = ["[{}={}]".format(key.replace("_", "-"), "; ".join(values)) for key, values 
                 in source_feature.qualifiers.items() if key != "mol_type"]
    topology = record.annotations.get("topology", "linear")
    if topology == "circular":
        modifiers.append("[topology=circular]")
    completedness = record.annotations.get("complete")
    if completedness:
        modifiers.append("[completedness=complete]")
        plasmid = record.annotations.get("plasmid")
        if plasmid:
            pass  # plasmid modifier is already added from source features.
        else:
            modifiers.append("[location=chromosome]")

    modifiers.append("[gcode=11]")
    modifiers = " ".join(modifiers)
    modifiers = modifiers.replace("[plasmid=", "[plasmid-name=")  # acceptable modifier is 'plasmid-name'
    return modifiers


def create_genbank_submission_file(genome, tbl_file, fsa_file, center_name, verbosity=2):
    # seq_rank = get_seq_rank(genome)
    R = deepcopy(list(genome.seq_records.values()))
    tbl_buffer = ""
    fsa_buffer = ""
    for record in R:
        tbl_buffer += ">Feature {}\n".format(record.name)
        for feature in record.features:
            if feature.type != "source":
                feature.assign_hit(verbosity=verbosity)
                tbl_buffer += feature_to_tbl(feature, center_name)

        record.description = create_fsa_modifiers(record)
        # ff_definition = get_ff_definition(record, seq_rank)
        # if ff_definition:
        #     record.description += " " + ff_definition
        fsa_buffer += record.format("fasta")
        record.description = ""
    with open(tbl_file, "w") as f:
        f.write(tbl_buffer)
    with open(fsa_file, "w") as f:
        f.write(fsa_buffer)



if __name__ == '__main__':
    pass
