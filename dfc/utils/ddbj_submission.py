#! /usr/bin/env python
# coding: UTF8

import sys
import os
from logging import getLogger
from copy import deepcopy
from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition
from dfc.genome import Genome
from dfc.utils.metadata_util import Metadata
from dfc import dfast_version
from Bio.SeqIO.InsdcIO import _insdc_location_string


class DDBJsubmission(object):

    def __init__(self, genome, config):
        self.genome = genome
        self.output_dir = os.path.join(genome.workDir, "ddbj")
        self.verbosity = config.DDBJ_SUBMISSION.get("output_verbosity", 2)
        self.enabled = config.DDBJ_SUBMISSION.get("enabled", True)
        self.logger = getLogger(__name__)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.metadata = self.load_metadata(config)

    def load_metadata(self, config):
        def _read_file(file_name):
            D = {}
            for line in open(file_name):
                key, value = line.strip("\n").split("\t")
                if value:
                    D[key] = value
            return D

        metadata_file = config.DDBJ_SUBMISSION.get("metadata_file")
        if metadata_file:
            if not os.path.exists(metadata_file):
                self.logger.error("Metadata file ({}) does not exist. Aborting...".format(metadata_file))
                exit(1)
            self.logger.info("Loading a metadata file from {}".format(metadata_file))

            return _read_file(metadata_file)
        else:
            return {}

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
            return "mss"

    def create_submission_file(self):
        if self.enabled:
            prefix = self.get_file_prefix()
            ann_file = os.path.join(self.output_dir, prefix + ".ann").replace(" ", "_").replace("(", "").replace(")", "")
            fasta_file = os.path.join(self.output_dir, prefix + ".fasta").replace(" ", "_").replace("(", "").replace(")", "")
            self.logger.info("Writing a DDBJ annotation file to {}".format(ann_file))
            self.logger.info("Writing a DDBJ sequence file to {}".format(fasta_file))
            create_ddbj_submission_file(self.genome, self.metadata, ann_file, fasta_file, self.verbosity)

            # experimental implementation to generate dfast_results.json
            from sys import version_info
            if version_info >= (3, 10):            
                try:
                    from dr_tools import drt_ann2json
                    json_out_dir = os.path.dirname(self.output_dir)
                    json_file = os.path.join(json_out_dir, "dfast_record.json")
                    drt_ann2json(ann_file, fasta_file, json_file, division="VRL")
                    self.logger.info(f"Converted MSS file into JSON. {ann_file} --> {json_file} [Experimental]")

                except ImportError:
                    self.logger.warning("[Experimental] DFAST Record JSON will not be generated. Please install dr_tools and use Python>=3.10.")
        else:
            self.logger.warning("'Generate DDBJ Submission file' is disabled. Skip processing")


# deprecated as replaced by _insdc_location_string.
def get_location_string(location):
    before = "<" if isinstance(location.start, BeforePosition) else ""
    after = ">" if isinstance(location.end, AfterPosition) else ""

    location_string = "{0}{1}..{2}{3}".format(before, int(location.start) + 1, after, int(location.end))
    if location.strand == -1:
        location_string = "complement({})".format(location_string)
    return location_string

def qualifier_to_table(qualifiers, key):
    ret = []
    for value in qualifiers.get(key, []):
        ret.append(["", "", "", key, str(value)])
    return ret


def feature_to_table(feature, rec_length):
    ret = []
    # location_string = get_location_string(feature.location)
    location_string = _insdc_location_string(feature.location, rec_length)
    for key in feature.qualifiers:
        if key == "translation":
            continue  # translation is not required for MSS
        if key == "EC_number":
            continue
        if key == "transl_except":
            value = feature.qualifiers[key][0]
            feature.qualifiers[key][0] = value.replace("complement(", "").replace("),", ",")
        ret += qualifier_to_table(feature.qualifiers, key)
    if not ret:
        ret.append(["", "", "", "", ""])  # add emtpty line
    ret[0][1] = feature.type
    ret[0][2] = location_string
    return ret


def source_feature_to_table(feature, record, seq_rank, rec_length, metadata):
    ret = []
    # location_string = get_location_string(feature.location)
    location_string = _insdc_location_string(feature.location, rec_length)
    for key in feature.qualifiers:
        ret += qualifier_to_table(feature.qualifiers, key)
    topology = record.annotations.get("topology", "linear")
    # organism = record.annotations.get("strain", "")
    strain = record.annotations.get("strain", "")
    # if strain:
    #     ret.append(["", "", "", "strain", strain])
    if seq_rank in ["contig", "scaffold"]:
        ret.append(["", "", "", "submitter_seqid", "@@[entry]@@"])
    # print(record.annotations)
    plasmid = record.annotations.get("plasmid")
    ret.append(create_ff_definiton(seq_rank, plasmid))
    ret += metadata.render_ex_source()
    ret[0][1] = feature.type
    ret[0][2] = location_string
    if topology == "circular":
        ret.insert(0, ["", "TOPOLOGY", "", "circular", ""])
    return ret


def create_ff_definiton(seq_rank, plasmid):
    assert seq_rank in ["complete", "scaffold", "contig"]
    if seq_rank == "complete":
        if plasmid:
            ff_definition = "@@[organism]@@ @@[strain]@@ plasmid @@[plasmid]@@ DNA, complete sequence"
        else:
            ff_definition = "@@[organism]@@ @@[strain]@@ DNA, complete genome"
        # if strain:
        #     ff_definition = "@@[organism]@@ @@[strain]@@ DNA, complete genome: @@[entry]@@"
        # else:
        #     ff_definition = "@@[organism]@@ DNA, complete genome: @@[entry]@@"
    else:
        ff_definition = "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"
        # if strain:
        #     ff_definition = "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"
        # else:
        #     ff_definition = "@@[organism]@@ DNA, @@[submitter_seqid]@@"
    return ["", "", "", "ff_definition", ff_definition]


def record_to_fasta(record):
    return ">{name}\n{seq}\n//\n".format(name=record.name, seq=str(record.seq))


def get_seq_rank(genome):
    if genome.complete:
        return "complete"
    else:
        features = [feature for feature in genome.features.values() if feature.type == "assembly_gap"]
        if len(features) > 0:
            return "scaffold"
        else:
            return "contig"

def create_ddbj_submission_file(genome, dict_metadata, ann_file, fasta_file, verbosity=2):
    R = deepcopy(list(genome.seq_records.values()))
    seq_rank = get_seq_rank(genome)  # complete, scaffold, or contig
    metadata = Metadata(dict_metadata)
    ann_buffer = metadata.render_common_entry(dfast_version=dfast_version, complete=genome.complete)
    fasta_buffer = ""
    for record in R:
        rec_length = len(record)
        entry_buffer = []
        for feature in record.features:
            feature.assign_hit(verbosity=verbosity)
            if feature.type == "source":
                entry_buffer += source_feature_to_table(feature, record, seq_rank, rec_length, metadata)
                # entry_buffer.append(create_ff_definiton(genome, seq_rank))
            else:
                entry_buffer += feature_to_table(feature, rec_length)
        entry_buffer[0][0] = record.name
        ann_buffer += entry_buffer
        fasta_buffer += record_to_fasta(record)
    with open(ann_file, "w") as f:
        for row in ann_buffer:
            f.write("\t".join(row) + "\n")
    with open(fasta_file, "w") as f:
        f.write(fasta_buffer)
