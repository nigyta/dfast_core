#! /usr/bin/env python
# coding: UTF8

import copy
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
        self.removed_by_rules = []  # descriptions of features removed by specific overlap rules
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
                    # misc_features converted from partial rRNAs still mask lower-priority features as rRNA
                    feature_type = "rRNA" if feature.annotations.get("putative_rrna") else feature.type
                    if feature_type not in feature_types:
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

            removed_features = [feature for feature in seq_record.features if feature.id in removed]
            seq_record.features = [feature for feature in seq_record.features if feature.id not in removed]
            return [describe_feature(feature, int(feature.location.start), int(feature.location.end), "overlapping")
                    for feature in removed_features]

        removed = []
        for seq_record in self.genome.seq_records.values():
            removed += _resolve_overlap(seq_record)
        self.genome.set_feature_dictionary()
        for description in removed:
            self.logger.debug("Removed overlapping feature: {}".format(description))
        if removed:
            self.logger.warning("Removed {} overlapping features not handled by specific overlap rules: {}".format(
                len(removed), "; ".join(removed)))
        else:
            self.logger.info("Removed 0 overlapping features.")

    def remove_partial_features(self):
        exempted_features = ["CRISPR", "misc_feature", "rRNA"]  # partial rRNAs are handled by adjust_rrna_features()
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

    def adjust_rrna_features(self):
        changed, converted = 0, 0
        for seq_record in self.genome.seq_records.values():
            gaps = [f for f in seq_record.features if f.type == "assembly_gap"]
            features = []
            for feature in seq_record.features:
                if feature.type != "rRNA":
                    features.append(feature)
                    continue
                location = str(feature.location)
                pieces, removed = adjust_rrna(feature, len(seq_record.seq), gaps)
                new_locations = [str(piece.location) for piece in pieces]
                if new_locations != [location]:
                    changed += 1
                    self.logger.debug("Adjusted rRNA {} at {}: {} -> {}".format(feature.id, feature.seq_id, location, ", ".join(new_locations) or "removed"))
                for description in removed:
                    self.logger.debug("Removed feature by overlap rule: {}".format(description))
                self.removed_by_rules += removed
                converted += sum(piece.type == "misc_feature" for piece in pieces)
                features += pieces
            seq_record.features = features
        self.genome.sort_features()  # split pieces may need reordering
        self.genome.set_feature_dictionary()
        self.logger.info("Adjusted locations of {} rRNA features at contig ends or assembly gaps. {} partial rRNA features are annotated as misc_feature.".format(changed, converted))

    def log_removed_by_rules(self):
        if self.removed_by_rules:
            self.logger.info("Removed {} features by specific overlap rules: {}".format(
                len(self.removed_by_rules), "; ".join(self.removed_by_rules)))
        else:
            self.logger.info("Removed 0 features by specific overlap rules.")

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

        # Overlap handling: specific rules run first, in feature type priority order
        # (assembly_gap > CRISPR > tRNA/tmRNA/rRNA > CDS), so that each rule sees the adjusted
        # higher-priority features. Add new rules (e.g. CDS vs tRNA) here in that order.
        # resolve_overlap() then runs as a fallback for overlaps that no specific rule handles.
        self.removed_by_rules = []
        self.adjust_rrna_features()
        self.log_removed_by_rules()

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
            # case: fix right partial CDS (##..>##) to intact CDS
            first3 =  str(seq[-3:].reverse_complement()).upper()
            if first3 in acceptable_codons and len(feature) % 3 == 0 and hit:
                _fix_partial(feature, seq)
            elif isinstance(hit, ProteinHit) and hit.description != "hypothetical protein":
                _to_misc_feature(feature, hit)
            elif isinstance(hit, NuclHit):
                _to_misc_feature(feature, hit)
        elif isinstance(hit, NuclHit):
                _to_misc_feature(feature, hit)


CONTIG_END_MARGIN = 10  # a partial rRNA ending within this many bases of a contig end is extended to the end
MIN_PIECE_LENGTH = 30  # pieces shorter than this after trimming or splitting by an overlap rule are removed
TRUNCATION_NOTES = {"contig end": "truncated at the contig end", "assembly gap": "truncated at an assembly gap"}


def adjust_rrna(feature, seq_len, gaps):
    """
    Fit an rRNA feature to the contig and the assembly gaps. The given feature is modified in place,
    and extra pieces are created when it is split.

    Barrnap partial hits ("aligned only N percent", i.e. < 80% of the expected length) are partial rRNAs:
    - They are converted to misc_feature ("putative rRNA, aligned only N percent ...").
    - An end less than CONTIG_END_MARGIN bases from a contig end is extended to that end and shown as partial (< or >).
    Other rRNAs stay rRNA and are never extended to a contig end.
    For all rRNAs:
    - An end overlapping an assembly gap is trimmed to the gap and shown as partial (< or >).
    - A gap inside the feature splits it into two pieces, each converted to misc_feature
      ("putative rRNA overlapping an assembly gap").
    - An end beyond the contig is clamped to the contig end and shown as partial.
    Truncated features get a note saying why.
    A feature lying entirely within gaps, and pieces shorter than MIN_PIECE_LENGTH after trimming
    or splitting, are removed.
    Returns (kept features, descriptions of removed features or pieces).
    """
    loc = feature.location
    notes = feature.qualifiers.get("note", [])
    aligned = [note for note in notes if note.startswith("aligned only")]
    start, end = int(loc.start), int(loc.end)
    # side value: None (exact) or the reason the side is partial
    left = "contig end" if isinstance(loc.start, BeforePosition) or start < 0 else None
    right = "contig end" if isinstance(loc.end, AfterPosition) or end > seq_len else None
    start, end = max(start, 0), min(end, seq_len)
    if aligned:
        if start < CONTIG_END_MARGIN:
            start, left = 0, "contig end"
        if seq_len - end < CONTIG_END_MARGIN:
            end, right = seq_len, "contig end"

    segments = [(start, end, left, right)]
    for gap in sorted(gaps, key=lambda g: int(g.location.start)):
        gap_start, gap_end = int(gap.location.start), int(gap.location.end)
        new_segments = []
        for s, e, lp, rp in segments:
            if gap_end <= s or e <= gap_start:
                new_segments.append((s, e, lp, rp))
                continue
            if s < gap_start:
                new_segments.append((s, gap_start, lp, "assembly gap"))
            if gap_end < e:
                new_segments.append((gap_end, e, "assembly gap", rp))
        segments = new_segments
    if not segments:
        return [], [describe_feature(feature, int(loc.start), int(loc.end), "inside an assembly gap")]

    split = len(segments) > 1
    original = (int(loc.start), int(loc.end))

    def _too_short(s, e):  # only pieces created by trimming or splitting are subject to MIN_PIECE_LENGTH
        return (split or (s, e) != original) and e - s < MIN_PIECE_LENGTH

    removed = [describe_feature(feature, s, e, "shorter than {} bp after trimming".format(MIN_PIECE_LENGTH))
               for s, e, _, _ in segments if _too_short(s, e)]
    segments = [seg for seg in segments if not _too_short(seg[0], seg[1])]
    other_notes = [note for note in notes if not note.startswith("aligned only")]
    product = feature.qualifiers.get("product", [""])[0]
    to_misc = split or bool(aligned)

    results = []
    for i, (s, e, lp, rp) in enumerate(segments):
        piece = feature if i == 0 else copy.deepcopy(feature)
        if i > 0:
            piece.id = "{}_{}".format(feature.id, i + 1)
        piece.location = FeatureLocation(BeforePosition(s) if lp else ExactPosition(s),
                                         AfterPosition(e) if rp else ExactPosition(e), loc.strand)
        truncation = [TRUNCATION_NOTES[reason] for reason in sorted({lp, rp} - {None})]
        if to_misc:
            prefix = "putative rRNA overlapping an assembly gap" if split else "putative rRNA"
            piece.qualifiers.pop("product", None)
            piece.qualifiers["note"] = ["{}, {}".format(prefix, aligned[0] if aligned else product)] + truncation + other_notes
            piece.type = "misc_feature"
            piece.annotations["putative_rrna"] = True
        elif truncation:
            piece.qualifiers["note"] = notes + truncation
        results.append(piece)
    return results, removed


def describe_feature(feature, start, end, reason):
    """One-line description of a (removed) feature for logs. start is 0-based, end is exclusive."""
    strand = "-" if feature.location.strand == -1 else "+"
    return "{} {} {}:{}..{}({}) [{}]".format(feature.id, feature.type, feature.seq_id, start + 1, end, strand, reason)


if __name__ == '__main__':
    pass
