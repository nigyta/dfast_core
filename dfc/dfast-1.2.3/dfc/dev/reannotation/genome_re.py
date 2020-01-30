#! /usr/bin/env python
# coding: UTF8

from collections import OrderedDict

import os
import math
from datetime import datetime
from logging import getLogger
import pickle
from copy import deepcopy
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, BeforePosition, AfterPosition
from Bio.Alphabet import IUPAC
from ...models.bio_feature import ExtendedFeature
from ...genome import Genome
logger = getLogger(__name__)


class GenomeReannotation(Genome):
    """
    Stores genome data as BioPython SeqRecord object

    TODO: sort and rename sequences
    TODO: force override option
    TODO: Minimum length
    """
    def __init__(self, config):

        super(GenomeReannotation, self).__init__(config)
        self.convert_features()

    def prepare_genome(self, config):
        query_genome_fasta = config.GENOME_FASTA
        use_original_name = config.GENOME_CONFIG.get("use_original_name", False)
        sort_sequence = config.GENOME_CONFIG.get("sort_sequence", True)
        minimum_length = config.GENOME_CONFIG.get("minimum_length", 0)
        R = [r for r in SeqIO.parse(open(query_genome_fasta), "genbank")]

        seq_dict = OrderedDict()

        if self.complete:
            logger.info("The query genome is treated as a complete genome with {} sequences.".format(len(R)))
            logger.info("'minimum_length' and 'sort_sequence' options will be ignored.")
            if len(R) != len(self.seq_topologies):
                logger.error("The numbers of sequences and seq_topologies do not match. Aborting...")
                logger.error("seq_topolgies: {}".format(self.seq_topologies))
                exit(1)
            elif len(R) != len(self.seq_types):
                logger.error("The numbers of sequences and seq_types do not match. Aborting...")
                logger.error("seq_types: {}".format(self.seq_types))
                exit(1)
            elif len(R) != len(self.seq_names) and not use_original_name:
                logger.error("The numbers of sequences and seq_names do not match. Aborting...")
                logger.error("seq_names: {}".format(self.seq_names))
                exit(1)
            if use_original_name:
                self.seq_names = []
            for i, r in enumerate(R):
                if use_original_name:
                    name = r.id
                    self.seq_names.append(name)
                else:
                    name = self.seq_names[i]
                    r.id = name
                if self.seq_types[i] == "plasmid":
                    logger.info("Sequence {num} (plasmid): {name} [{topology}]".format(num=(i + 1), name=name,
                                                                                      topology=self.seq_topologies[i]))
                else:
                    logger.info("Sequence {num}: {name}, [{topology}]".format(num=(i + 1), name=name,
                                                                             topology=self.seq_topologies[i]))
                r.name = name  # will be used in the LOCUS line in a gbk format
                r.description = ""  # will be used as a definition in a gbk file and as a fasta header. i.e. ">[seq.id] [seq.description]"
                seq_dict[r.id] = r
        else:
            logger.info("The query genome is treated as a draft genome with {} sequences.".format(len(R)))

            if sort_sequence:
                logger.info("Sequences are sorted by length (from longer to shorter).")
                R.sort(key=lambda x: -len(x))

            if minimum_length > 0:
                logger.info("Sequences shorter than {} will be eliminated.".format(minimum_length))
                R = [r for r in R if len(r) >= minimum_length]
                assert len(R) > 0
            if self.seq_names and not use_original_name:
                prefix = self.seq_names[0]  # First word is used as a seq prefix
                prefix = prefix.replace(" ", "_").replace(">", "_")
                prefix = prefix or "sequence"
            else:
                prefix = "sequence"
            digit = int(math.log10(len(R))) + 1
            if not use_original_name:
                logger.info("Sequences will be renamed as {0}, {1}...".format(prefix + str(1).zfill(digit), prefix + str(2).zfill(digit)))
            for i, r in enumerate(R, 1):
                if use_original_name:
                    name = r.id
                else:
                    name = prefix + str(i).zfill(digit)
                r.id = name
                r.name = name  # will be used in the LOCUS line in a gbk format
                r.description = ""  # will be used as a definition in a gbk file and as a fasta header. i.e. ">[seq.id] [seq.description]"
                seq_dict[r.id] = r

        # write a fasta file.
        output_directory = os.path.join(self.workDir, "input")
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        output_genome_fasta = os.path.join(output_directory, "genome.fna")
        with open(output_genome_fasta, "w") as f:
            SeqIO.write(R, f, "fasta")
        return output_genome_fasta, seq_dict

    def convert_features(self):
        def convert(seq_feature, feature_id, seq_id=None, primary_hit=None, secondary_hits=None, annotations=None):
            ext_feature = ExtendedFeature(location=seq_feature.location,
                                          type=seq_feature.type,
                                          location_operator=seq_feature.location_operator,
                                          strand=seq_feature.strand,
                                          id=feature_id,
                                          qualifiers=seq_feature.qualifiers,
                                          sub_features=None,
                                          ref=seq_feature.ref,
                                          ref_db=seq_feature.ref_db,
                                          primary_hit=primary_hit,
                                          seq_id=seq_id,
                                          secondary_hits=secondary_hits,
                                          annotations=annotations)
            if ext_feature.type == "CDS":
                note = ""
                if "product" in ext_feature.qualifiers:
                    note += "old_product: " + ext_feature.qualifiers["product"][0]
                if "protein_id" in ext_feature.qualifiers:
                    note += ", old_protein_id: " + ext_feature.qualifiers["protein_id"][0]
                    note = note.strip(", ")
                    del ext_feature.qualifiers["protein_id"]
                if "transl_table" not in ext_feature.qualifiers:
                    ext_feature.qualifiers["transl_table"] = [11]
                if "codon_start" not in ext_feature.qualifiers:
                    ext_feature.qualifiers["codon_start"] = [1]
                ext_feature.qualifiers.setdefault("note", []).append(note)
            if "old_locus_tag" in ext_feature.qualifiers:
                del ext_feature.qualifiers["old_locus_tag"]
            if "locus_tag" in ext_feature.qualifiers:
                ext_feature.qualifiers["old_locus_tag"] = ext_feature.qualifiers["locus_tag"]
            return ext_feature

        for r in self.seq_records.values():
            features = []
            seq_id = r.id
            for i, f in enumerate(r.features, start=1):
                if f.type == "gene":
                    continue
                feature_id = seq_id + "_" + str(i)
                ext_feature = convert(f, feature_id, seq_id=seq_id)
                features.append(ext_feature)
                self.features[feature_id] = ext_feature
            r.features = features


#     def add_source_information(self):
#         today = datetime.now().strftime('%d-%b-%Y').upper()
#         source = self.organism + " strain " + self.strain if self.strain else self.organism
#         if self.complete:
#             assert len(self.seq_types) == len(self.seq_topologies) == len(self.seq_names) == len(self.seq_records)
#             for i, record, seq_name, seq_type, seq_topology in zip(range(len(self.seq_records)), self.seq_records.values(),
#                                                                    self.seq_names, self.seq_types, self.seq_topologies):
#                 annotations = {
#                     "organism": self.organism, "source": source, "strain": self.strain,
#                     "complete": self.complete, "date": today, "topology": seq_topology,
#                     "data_file_division": "BCT",
#                     # "sequence_version": 1, "data_file_division": "BCT",
#                     # "taxonomy":['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Paphiopedilum'],
#                 }
#                 if seq_type == "plasmid" or seq_type == "p":
#                     annotations["plasmid"] = record.name
#                 logger.debug("DEBUG: " + str(annotations))
#                 record.annotations = annotations
#         else:
#             for i, record in enumerate(self.seq_records.values()):
#                 record.annotations = {
#                     "date": today, "data_file_division": "BCT", "topology": "linear",
#                     "organism": self.organism, "source": source, "strain": self.strain, "complete": self.complete
#                 }
# 
#     def load_metadata(self, config):
#         def _read_file(file_name):
#             D = {}
#             for line in open(file_name):
#                 key, value = line.strip("\n").split("\t")
#                 if value:
#                     D[key] = value
#             return D
# 
#         metadata_file = config.GENOME_SOURCE_INFORMATION.get("metadata_file")
#         if metadata_file:
#             if not os.path.exists(metadata_file):
#                 logger.error("Metadata file ({}) does not exist. Aborting...".format(metadata_file))
#                 exit(1)
#             logger.info("Loading a metadata file from {}".format(metadata_file))
# 
#             return _read_file(metadata_file)
#         else:
#             return {}
# 
#     def sort_features(self):
#         for record in self.seq_records.values():
#             record.features.sort(key=lambda x: x.location.start)
# 
#     def set_feature_dictionary(self):
#         self.features = OrderedDict()
#         for record in self.seq_records.values():
#             for feature in record.features:
#                 self.features[feature.id] = feature
# 
#     def add_source_features(self):
# 
#         for record in self.seq_records.values():
#             location = FeatureLocation(0, len(record), strand=1)
#             source_feature = ExtendedFeature(location=location, type="source", id=record.id, seq_id=record.id)
#             qualifiers = OrderedDict()
#             qualifiers["mol_type"] = [record.annotations.get("mol_type", "genomic DNA")]
#             qualifiers["organism"] = [record.annotations.get("organism", "")]
#             strain = record.annotations.get("strain")
#             if strain:
#                 qualifiers["strain"] = [strain]
#             plasmid = record.annotations.get("plasmid")
#             if plasmid:
#                 qualifiers["plasmid"] = [plasmid]
#             for key, values in self.additional_modifiers.items():
#                 for value in values:
#                     qualifiers.setdefault(key, []).append(value)
#             source_feature.qualifiers = qualifiers
#             if len(record.features) > 0 and record.features[0].type == "source":
#                 record.features[0] = source_feature
#             else:
#                 record.features.insert(0, source_feature)
# 
#     def to_genbank(self, file_name, verbosity=2):
#         logger.info("Writing a GenBank format file to {0}. (verbosity level={1})".format(file_name, verbosity))
#         R = deepcopy(list(self.seq_records.values()))
#         for record in R:
#             for feature in record.features:
#                 feature.assign_hit(verbosity=verbosity)
#         with open(file_name, "w") as f:
#             SeqIO.write(R, f, "genbank")
# 
#     def to_pickle(self, file_name):
#         with open(file_name, "wb") as f:
#             pickle.dump(self, f)
# 
# 
# allowed_qualifiers = ["bio_material",  "collected_by", "collection_date", "country", "culture_collection",
#                       "ecotype", "identified_by", "isolation_source", "host", "note",
#                       "serotype", "serovar", "sub_species", "sub_strain", "type_material", "variety",
#                       # "organism", "strain", "plasmid"
#                       ]
# 
# 
# def parse_additional_qualifier(qualifiers):
#     allowed_qualifiers_mod = [x.replace("_", "-") for x in allowed_qualifiers if "_" in x]
#     allowed_qualifiers_mod += allowed_qualifiers
#     ret = OrderedDict()
#     qualifiers = qualifiers.strip("; ").split(";")
#     for qualifier in qualifiers:
#         if qualifier.count("=") != 1:
#             continue
#         key, value = qualifier.strip().split("=")
#         key, value = key.strip(), value.strip()
#         if key not in allowed_qualifiers_mod:
#             continue
#         ret.setdefault(key, []).append(value)
#     return ret
# 
# if __name__ == '__main__':
#     from logging import getLogger, StreamHandler, FileHandler, DEBUG, WARN, INFO, Formatter
#     genome = Genome("test/genome.fna")
#     genome.toGenBank("test/genome.gbk")
#     genome.to_pickle("test/genome.pickle")