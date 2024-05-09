#! /usr/bin/env python
# coding: UTF8

from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, BeforePosition, AfterPosition
from logging import getLogger
from dfc.models.hit import ProteinHit, Hit, NuclHit, HmmHit, CddHit

class FeatureUtil(object):

    def __init__(self, genome, config):
        self.genome = genome
        self.logger = getLogger(__name__)
        cfg = config.FEATURE_ADJUSTMENT
        self.enable_remove_partial = cfg.get("remove_partial_features", True)

        self.enable_remove_overlapping = cfg.get("remove_overlapping_features", True)
        self.priority_for_overlapping_features = cfg.get("feature_type_priority", ["assembly_gap", "CRISPR", ("tmRNA", "tRNA", "rRNA"), "CDS"])

        self.enable_merge_cds = cfg.get("merge_cds", False)
        self.priority_for_merge_cds = cfg.get("tool_type_priority", {"MGA": 0, "Prodigal": 1})
        self.show_settings()

    def show_settings(self):
        if self.enable_remove_partial:
            self.logger.info("Remove_Partial_Feature is enabled.")

        if self.enable_remove_overlapping:
            self.logger.info("Remove_Overlapping_Feature is enabled. " +
                             "Priority: {}".format(self.priority_for_overlapping_features))

        if self.enable_merge_cds:
            self.logger.info("Merge_CDS is enabled. CDSs predicted from different tools are merged.\n" +
                             "[WARNING] Merge_CDS is a preliminary version. Tools with a lower number have higher priotity.\n" +
                             "Priority: {}".format(self.priority_for_merge_cds))

    def resolve_overlap(self):

        def _resolve_overlap(seq_record, threshold=10):
            removed = []
            masked_seq = str(seq_record.seq)
            for feature_types in self.priority_for_overlapping_features:
                if not isinstance(feature_types, tuple):
                    feature_types = (feature_types,)
                original_seq = masked_seq
                for feature in seq_record.features:
                    if feature.type not in feature_types:
                        continue
                    extracted_seq = str(feature.extract(original_seq))
                    #             print(masked_seq)
                    #             print(original_seq)
                    #             print(extracted_seq)
                    #             print(feature)
                    if 100.0 * extracted_seq.count("x") / len(extracted_seq) > threshold:
                        removed.append(feature.id)
                    else:
                        masked_seq = masked_seq[:feature.location.start] + "x" * len(extracted_seq) + masked_seq[
                                                                                                      feature.location.end:]

            seq_record.features = [feature for feature in seq_record.features if feature.id not in removed]
            for feature_id in removed:
                feature = self.genome.features[feature_id]
                self.logger.debug("Removing overlapping feature: {} at {}:{}".format(feature.id, feature.seq_id, str(feature.location)))
            return len(removed)

        count_removed = 0
        for seq_record in self.genome.seq_records.values():
            count_removed += _resolve_overlap(seq_record)
        self.genome.set_feature_dictionary()
        self.logger.info("Removed {} overlapping features.".format(count_removed))

    def remove_partial_features(self):
        exempted_features = ["CRISPR", "misc_feature"]
        removed = []

        for feature in self.genome.features.values():
            seq = self.genome.seq_records[feature.seq_id].seq
            if feature.location.start < 10 or feature.location.end > len(seq) - 10:
                # print(type(seq))
                # print(feature)
                # print(feature.annotations)
                # print(self.genome.seq_records)
                self.fix_partial_CDS(feature)
                # print(feature)
                # print(feature.annotations)
            if feature.type in exempted_features:
                continue
            start = feature.location.start
            end = feature.location.end
            if isinstance(start, BeforePosition) or isinstance(end, AfterPosition):
                self.logger.debug("Removing partial feature: {} at {}:{}".format(feature.id, feature.seq_id, str(feature.location)))
                # temporary disabled for dev
                removed.append(feature.id)
        for seq_record in self.genome.seq_records.values():
            seq_record.features = [feature for feature in seq_record.features if feature.id not in removed]
        self.genome.set_feature_dictionary()
        self.logger.info("Removed {} partial features.".format(len(removed)))

    def merge_cds(self):
        def _have_same_cds_end(one, other):
            """
            check if the two cds shares common stop position
            """
            if one.strand == 1:
                return one.seq_id == other.seq_id and one.strand == other.strand and one.location.end == other.location.end
            else:
                return one.seq_id == other.seq_id and one.strand == other.strand and one.location.start == other.location.start

        def _get_cds_start(f):
            if f.location.strand == 1:
                return f.location.start
            else:
                return f.location.end

        def _get_cds_end(f):
            if f.location.strand == 1:
                return f.location.end
            else:
                return f.location.start

        def _merge(feature_list):
            L = []
            for f1 in feature_list:
                if len(L) == 0:
                    L.append([f1])
                else:
                    for f2 in L[-1]:
                        if _have_same_cds_end(f1, f2):
                            if _get_cds_start(f1) == _get_cds_start(f2):
                                # f1.extended_attributes.setdefault("alt_features", []).append(f2)
                                break
                        else:
                            L.append([f1])
                            break
                    else:
                        L[-1].append(f1)
            return L

        def _get_priority(feature):
            return self.priority_for_merge_cds[feature.id.split("_")[0]]
            # return priority.get(feature.id.split("_")[0], 9)

        def _choose_representative_location(features):
            if len(features) == 1:
                return features[0]
            features.sort(key=lambda x: _get_priority(x))

            features_with_rbs = [x for x in features if x.annotations.get("rbs")]
            if len(features_with_rbs) > 0:
                return features_with_rbs[0]
            return features[0]

        # main part starts from here
        # extract cds features
        cds_features = [feature for feature in self.genome.features.values() if feature.type == "CDS"]
        non_cds_features = [feature for feature in self.genome.features.values() if feature.type != "CDS"]
        self.logger.info("Merging CDS features from different prediction tools. Start with {} CDSs.".format(len(cds_features)))

        # reset seq features
        for record in self.genome.seq_records.values():
            record.features = []

        # add representative features
        cnt = 0
        for features in _merge(cds_features):
            representative_feature = _choose_representative_location(features)
            seq_id = representative_feature.seq_id
            self.genome.seq_records[seq_id].features.append(representative_feature)
            cnt += 1

        # add non-cds features
        for feature in non_cds_features:
            seq_id = feature.seq_id
            self.genome.seq_records[seq_id].features.append(feature)

        # sort and reset dictionary
        self.genome.sort_features()
        self.genome.set_feature_dictionary()
        self.logger.info("Merged CDS features. {} CDSs in total.".format(cnt))

    def execute(self):
        # if self.enable_remove_partial:
        #     self.remove_partial_features()

        if self.enable_remove_overlapping:
            self.resolve_overlap()

        if self.enable_merge_cds:
            self.merge_cds()

    def execute_remove_partial(self):
        if self.enable_remove_partial:
            self.remove_partial_features()


    def fix_partial_CDS(self, feature, min_length=500):
        """
        fix pseudo-partial CDS predicted at the end of the contig, but is likely to be intact.
        The one that meet the following conditions will be retained
            - Starts at the coordinate 1 (+ strand) or at the end of contig (- strand)
            - 5'-end is missing and 3'-end is intact (+ strand) or the opposited case (- strand) 
            - length is multiple of 3 and > min_length
            - the CDS has significant protein hit or nucl hit.

        partial_flag 10, strand +, len>=500
        partial_flag 01, strand -, len>=500

        if the first 3 nuleotides are in [ATG, GTG, TTG], the CDS is fixed as intact,
        otherwise, it is annotated as misc_feature
        """
        def _fix_partial(feature, seq):
            self.logger.warning("Fixed partial CDS predicted at %s:%s", feature.seq_id, str(feature.location))
            if feature.strand == 1:
                feature.location = FeatureLocation(ExactPosition(0), feature.location.end, feature.strand)
            else:
                feature.location = FeatureLocation(feature.location.start, ExactPosition(len(seq)), feature.strand)
            feature.annotations["partial_flag"] = "00"
            if "partial" in feature.annotations:
                del feature.annotations["partial"]

        def _to_misc_feature(feature, hit):
            if not feature.type == "CDS":
                return  # Do nothing for features other than CDS
            self.logger.warning("Partial CDS predicted at %s:%s will be annotated as misc_feature", feature.seq_id, str(feature.location))
            feature.type = "misc_feature"
            if isinstance(hit, ProteinHit):
                note = f"partial CDS similar to {hit.id}:{hit.description}"
            elif isinstance(hit, NuclHit):
                note = "partial CDS, " + hit.model.info()
            else:
                note = "partial CDS"
            feature.qualifiers.setdefault("note", []).append(note)
            for key in ["product", "translation", "transl_table", "codon_start"]:
                if key in feature.qualifiers:
                    del feature.qualifiers[key]

        def _get_hit(feature):
            if feature.primary_hit:
                return feature.primary_hit
            elif feature.secondary_hits and isinstance(feature.secondary_hits[0], ProteinHit):
                return feature.secondary_hits[0]
            else:
                return None

        self.logger.debug("Fixing partial CDS: %s", feature.id)
        acceptable_codons = ["ATG", "GTG", "TTG"]
        seq = self.genome.seq_records[feature.seq_id].seq
        hit = _get_hit(feature)
        # self.logger.debug("Protein hit: %s", hit)
        if int(feature.location.start) == 0 and feature.strand == 1 and feature.annotations.get("partial_flag", "00") == "10" and len(feature) >= min_length:
            # case: fix left partial CDS (<1..##) to intact CDS
            first3 =  str(seq[:3]).upper()
            if first3 in acceptable_codons and len(feature) % 3 == 0 and hit:
                _fix_partial(feature, seq)
            elif isinstance(hit, ProteinHit) and hit.description != "hypothetical protein":
                _to_misc_feature(feature, hit)
            elif isinstance(hit, NuclHit):
                _to_misc_feature(feature, hit)
        elif int(feature.location.end) == len(seq) and feature.strand == -1 and feature.annotations.get("partial_flag", "00") == "01" and len(feature) >= min_length:
            # case: fix right partial CDS (##..>##)　to intact CDS
            first3 =  str(seq[-3:].reverse_complement()).upper()
            if first3 in acceptable_codons and len(feature) % 3 == 0 and hit:
                _fix_partial(feature, seq)
            elif isinstance(hit, ProteinHit) and hit.description != "hypothetical protein":
                _to_misc_feature(feature, hit)
            elif isinstance(hit, NuclHit):
                _to_misc_feature(feature, hit)
        elif isinstance(hit, NuclHit):
                _to_misc_feature(feature, hit)

if __name__ == '__main__':
    pass
