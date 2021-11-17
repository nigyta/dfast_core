#! /usr/bin/env python
# coding: UTF8

import re
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from .base_tools import StructuralAnnotationTool
from ..models.bio_feature import ExtendedFeature

from logging import getLogger
logger = getLogger(__name__)

class BlastFeatureN(StructuralAnnotationTool):
    """
    Import features from GFF file
    This class does not use external commands.

    Tool type: GFF_import

    """
    version = None
    TYPE = "feature"
    NAME = "BlastFeatureN"
    VERSION_CHECK_CMD = ["blastn", "-version"]
    VERSION_PATTERN = r"blastn: (.+)\+"
    SHELL = True

    def __init__(self, options=None, workDir="OUT"):
        super(BlastFeatureN, self).__init__(options, workDir)
        self.gbk_file_name = options.get("gbk_file_name")
        self.targets = options.get("targets", ["CDS"])
        self.targets_no_stranded = options.get("targets_non_stranded", ["stem_loop"])
        self.imported_features = {}
        self.logger.info("Following features will be imported from GFF: " + ", ".join(self.targets))
        self.annotated_features = {}

    def get_createdb_command(self, db_fasta_file, db_name, dbtype):
        return ["makeblastdb", "-dbtype", dbtype, "-in", db_fasta_file, "-out", db_name, "-hash_index", "-parse_seqids"]

    def get_aln_command(self, query_file, db_name, result_file):
        outfmt = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'"
        return ["blastn", "-task", "blastn", "-query", query_file, "-db", db_name, "-out", result_file,
                "-outfmt", outfmt, "-num_alignments 10 -max_hsps 5"]



    def run(self):
        def _get_query_genome_dict():
            query_genome_fasta = os.path.join(self.workDir, "input", "genome.fna")
            return SeqIO.to_dict(SeqIO.parse(query_genome_fasta, "fasta"))

        def _create_reference_fasta(dict_nuc_fasta, out_file_name):
            ret = [">{}\n{}\n".format(k, v) for k, v in dict_nuc_fasta.items()]
            ret = "".join(ret)
            with open(out_file_name, "w") as f:
                f.write(ret)

        def _parse_result():
            with open(self.outputFile) as f:
                for line in f:
                    qaccver, saccver, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = line.strip("\n").split("\t")
                    qstart, qend, sstart, send, qlen = int(qstart), int(qend), int(sstart), int(send), int(qlen)
                    feature_type = qaccver.split(".")[0]
                    if feature_type in self.targets_no_stranded:
                        left, right, strand = min(sstart, send), max(sstart, send), 1
                    elif sstart <= send:
                        left, right, strand = sstart, send, 1
                    else:
                        left, right, strand = send, sstart, -1
                    
                    # partial check
                    if strand == 1:
                        left_flag = "0" if qstart == 1 else "1"
                        right_flag = "0" if qend == qlen else "1"
                        partial_flag = left_flag + right_flag
                    else:
                        right_flag = "0" if qstart == 1 else "1"
                        left_flag = "0" if qend == qlen else "1"
                        partial_flag = left_flag + right_flag

                    location = self.getLocation(left, right, strand, partial_flag)
                    annotations = {}
                    # todo: add partial feature, codon_start qualifier
                    feature = ExtendedFeature(location=location, type=feature_type,
                                              # id=qaccver.format(self.__class__.__name__, i),
                                              id=qaccver,
                                              seq_id=saccver, annotations=annotations)
                    yield feature

        feature_num = 0
        dict_ref_features = {}
        dict_nuc_fasta = {}
        dict_prot_fasta = {}

        # prepare reference sequences
        for r in SeqIO.parse(self.gbk_file_name, "genbank"):
            for f in r.features:
                if f.type not in self.targets:
                    continue
                feature_num += 1
                feature_id = f.type + "." + str(feature_num)
                dict_ref_features[feature_id] = f
                if isinstance(f.location, CompoundLocation):
                    for i, loc_part in enumerate(f.location.parts, 1):
                        sub_feature_id = feature_id + "/" + str(i)
                        dict_nuc_fasta[sub_feature_id] = str(loc_part.extract(r).seq)
                else:
                    dict_nuc_fasta[feature_id] = str(f.extract(r).seq)
                if "translation" in f.qualifiers:
                    dict_prot_fasta[feature_id] = f.qualifiers["translation"][0]

        reference_nuc_fasta = os.path.join(self.workDir, "structural", "{0}.ref.nucl.fasta".format(self.__class__.__name__))
        _create_reference_fasta(dict_nuc_fasta, reference_nuc_fasta)
        reference_prot_fasta = os.path.join(self.workDir, "input", "ref.prot.fasta")
        _create_reference_fasta(dict_prot_fasta, reference_prot_fasta)

        # prepare BlastDB (Input genome is used as BlastDB)
        db_prefix = os.path.join(self.workDir, "structural", "{0}.genome.db".format(self.__class__.__name__))
        cmd = self.get_createdb_command(self.genomeFasta, db_prefix, dbtype="nucl")
        self.executeCommand(cmd, shell=self.__class__.SHELL)

        cmd = self.get_createdb_command(reference_prot_fasta, reference_prot_fasta, dbtype="prot")
        self.executeCommand(cmd, shell=self.__class__.SHELL)

        # run Blast alignment
        cmd = self.get_aln_command(reference_nuc_fasta, db_prefix, self.outputFile)
        self.executeCommand(cmd, shell=self.__class__.SHELL)



        dict_query_genome = _get_query_genome_dict()

        for feature in _parse_result():
            sequence = feature.seq_id
            self.annotated_features.setdefault(sequence, []).append(feature)

        join_features(self.annotated_features)
        add_translation_qualifier(self.annotated_features, dict_query_genome)
        transfer_qualifiers_from_reference(self.annotated_features, dict_ref_features)

    def getFeatures(self):

        return self.annotated_features



def join_features(dict_features):
    for seq_id in dict_features:
        features = dict_features[seq_id]
        dict_feature_by_id = {}
        for feature in features:
            feature_id = feature.id.split("/")[0]  # example of feature ID: CDS.1 CDS.2/1, CDS.2/2, CDS.3, ...  
            dict_feature_by_id.setdefault(feature_id, []).append(feature)
        tmp_features = []
        for feature_id, features_with_same_id in dict_feature_by_id.items():
            if len(features_with_same_id) == 1:
                tmp_features.append(features_with_same_id[0])
                # todo: add orphan feature
            elif len(features_with_same_id) > 1:
                compound_location = sum([f.location for f in features_with_same_id])
                # print(compound_location)
                updated_feature = features_with_same_id[0]
                updated_feature.location = compound_location
                updated_feature.id = feature_id
                tmp_features.append(updated_feature)
            else:
                # features_with_same_id is empty
                raise AssertionError
        dict_features[seq_id] = tmp_features


def add_translation_qualifier(dict_features, dict_query_genome):
    for seq_id, features in dict_features.items():
        seq_record = dict_query_genome[seq_id]
        for feature in features:
            # todo: add transl_table
            if feature.type == "CDS":
                # translation = feature.translate(seq_record.seq)
                print(len(feature.location))
                if len(feature.location) % 3 == 0:
                    translation = feature.translate(seq_record.seq)
                else:
                    translation = feature.translate(seq_record.seq, cds=False, to_stop=False)
                feature.qualifiers["translation"] = [translation]


def transfer_qualifiers_from_reference(dict_features, dict_ref_features):
    target_qualifiers = ["product", "gene", "function", "ribosomal_slippage"]
    for seq_id, features in dict_features.items():
        for feature in features:
            ref_feature = dict_ref_features.get(feature.id)
            if not ref_feature:
                logger.warning("No reference feature for {}".format(feature.id))
                continue
            for key, value in ref_feature.qualifiers.items():
                # todo: biopython currently cannot handle boolean qualifier such as ribosomal_slippage
                if key in target_qualifiers:
                    feature.qualifiers[key] = value
