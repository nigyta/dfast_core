#! /usr/bin/env python
# coding: UTF8

"""Map nucleotide references (CARD/VFDB) onto the genome with blastn and annotate the mapped regions."""

import json
import os
import re
from collections import Counter
from itertools import count
from Bio.SeqFeature import FeatureLocation
from .baseComponent import BaseAnnotationComponent
from .NuclSearch import models
from ..models.bio_feature import ExtendedFeature
from ..models.hit import NuclHit, PseudoGene
from ..utils.mapping_util import (BLAST_OUTFMT, read_blast_table, is_full_length, group_by_locus, is_valid_reference,
                                  resolve_region, classify_cds, same_locus)

PSEUDOGENE_ACTIONS = {"frameshift": "pseudogene: frameshift", "internal_stop": "pseudogene: internal stop codon",
                      "stop_lost": "pseudogene: stop codon lost", "no_start": "pseudogene: no start codon"}


def _location(region):
    return FeatureLocation(region.left0, region.right, strand=region.strand)


def _start(feature):
    """Coordinate of the start codon side: left end on the plus strand, right end on the minus strand."""
    return int(feature.location.start) if feature.location.strand == 1 else int(feature.location.end)


def new_cds(region, feature_id, transl_table):
    return ExtendedFeature(location=_location(region), type="CDS", id=feature_id, seq_id=region.seq_id,
                           qualifiers={"translation": [region.translation], "transl_table": [str(transl_table)],
                                       "codon_start": [1], "note": list(region.notes)},
                           annotations={"partial_flag": "00", "nucl_mapping": "rescued"})


def new_pseudogene(region, feature_id):
    return ExtendedFeature(location=_location(region), type="CDS", id=feature_id, seq_id=region.seq_id,
                           qualifiers={"pseudogene": ["unknown"], "note": list(region.notes)},
                           annotations={"partial_flag": "00", "nucl_mapping": "pseudogene"})


def relocate_cds(cds, region, ref_label):
    """Move the start codon of a same-gene CDS to that of the mapped region; the CDS keeps its id."""
    old_start = _start(cds)
    old_length = len(cds.location)
    cds.location = _location(region)
    shift = abs(_start(cds) - old_start)
    direction = "upstream" if len(cds.location) > old_length else "downstream"
    cds.qualifiers["translation"] = [region.translation]
    cds.qualifiers["codon_start"] = [1]
    cds.qualifiers.setdefault("note", []).extend(region.notes)
    cds.qualifiers["note"].append(f"start codon moved {shift} bp {direction} by mapping of {ref_label}")
    cds.annotations["partial_flag"] = "00"
    cds.annotations.pop("partial", None)
    cds.annotations["nucl_mapping"] = "start corrected"


class NuclMapping(BaseAnnotationComponent):

    def __init__(self, genome, options, workDir, CPU):
        options["aligner"] = "blastn"
        super(NuclMapping, self).__init__(genome, options, workDir, CPU)
        self.pipeline_dir = workDir
        self.database = options.get("database", "")
        self.db_name = options.get("db_name", "")
        self.pident_cutoff = options.get("pident_cutoff", 90)
        self.evalue_cutoff = options.get("evalue_cutoff", 1e-10)
        self.overlap_tolerance = options.get("overlap_tolerance", 30)
        self.transl_table = options.get("transl_table", 11)
        self.blast_task = options.get("blast_task", "dc-megablast")
        self.pseudogene_pident_cutoff = options.get("pseudogene_pident_cutoff", 95)
        self.aligner_options = options.get("aligner_options", "")
        self.result_file = os.path.join(self.workDir, "mapping.out")
        self.hits_json = {}
        self.get_db_version()
        self.load_references()

    def get_db_version(self):
        version_file = self.database + ".version"
        if os.path.exists(version_file):
            db_version = open(version_file).read().strip()
            self.logger.info(f"Database version: {self.db_name} {db_version}")
        else:
            self.logger.warning(f"Database version could not be found: {self.db_name} UNDETERMINED")

    def load_references(self):
        self.references = models[self.db_name].load_from_tsv_file(self.database + ".nucl.ref")
        self.ref_validity = {ref_id: is_valid_reference(model.nucl_seq, self.transl_table)
                             for ref_id, model in self.references.items()}

    def genome_fasta(self):
        return os.path.join(self.pipeline_dir, "input", "genome.fna")

    def create_db(self):
        db_prefix = os.path.join(self.workDir, "genome")
        self.executeCommand(self.aligner.format_db_command(self.genome_fasta(), db_prefix))
        return db_prefix

    def blast_command(self, db_prefix):
        return ["blastn", "-task", self.blast_task, "-query", self.database + ".nucl.fasta", "-db", db_prefix,
                "-out", self.result_file, "-outfmt", "'" + BLAST_OUTFMT + "'",
                "-evalue", str(self.evalue_cutoff), "-perc_identity", str(self.pident_cutoff),
                "-max_target_seqs", "1000", "-max_hsps", "50", "-num_threads", str(self.CPU)] \
            + self.aligner_options.split()

    def select_hits(self):
        hits = []
        for hit in read_blast_table(self.result_file):
            if not is_full_length(hit, self.pident_cutoff):
                continue
            if hit.seq_id not in self.genome.seq_records:
                self.logger.warning(f"Mapping hit {hit.ref_id} refers to an unknown sequence ({hit.seq_id}); ignored.")
                continue
            hits.append(hit)
        return group_by_locus(hits)

    def apply(self, groups):
        """Rescue, relocate or call pseudogenes from the full-length mappings, one group (locus) at a time."""
        map_ids = [int(m.group(1)) for m in map(re.compile(r"MAP_(\d+)").fullmatch, self.genome.features) if m]
        self.map_ids = count(max(map_ids, default=0) + 1)
        self.counts = Counter()
        # Resolved before acting, at the final location of a region extended to the next stop codon.
        resolved = [r for r in map(self.resolve_group, groups) if r]
        # A CDS that is the same gene of some group is never removed as overlapping by another group (operon
        # neighbours overlapping by more than overlap_tolerance).
        self.same_genes = set()
        for region, _ in resolved:
            features = self.genome.seq_records[region.seq_id].features
            self.same_genes.update(cds.id for cds in features if cds.type == "CDS"
                                   and classify_cds(cds, region, self.overlap_tolerance) == "same_gene")
        for region, group in resolved:
            self.apply_group(region, group)
        self.genome.sort_features()
        self.genome.set_feature_dictionary()
        c = self.counts
        self.logger.info(f"{self.db_name} mapping: {c['rescued']} rescued, {c['start corrected']} start corrected, "
                         f"{c['pseudogene']} pseudogenes, {c['annotated']} annotated, {c['removed']} CDSs removed")

    def rank(self, region, hit):
        """Members with a valid reference resolving to an intact CDS first, then by identity: a longer, less identical
        reference must not turn an intact gene into a pseudogene. Invalid references rank after valid ones."""
        valid = self.ref_validity[hit.ref_id][0]
        intact = region is not None and valid and region.kind == "valid"
        return region is not None, intact, valid, hit.pident, hit.bitscore

    def resolve_group(self, group):
        """(region of the primary hit, members with the primary first), or None if no member region resolves."""
        record = self.genome.seq_records[group[0].seq_id]
        ranked = sorted(((resolve_region(record.seq, h, self.transl_table, self.ref_validity[h.ref_id][1]), h)
                         for h in group), key=lambda rh: self.rank(*rh), reverse=True)
        region = ranked[0][0]
        return None if region is None else (region, [h for _, h in ranked])

    def apply_group(self, region, group):
        record = self.genome.seq_records[region.seq_id]
        primary = group[0]
        valid = self.ref_validity[primary.ref_id][0]
        location = (region.left0, region.right)
        for feature in record.features:  # a feature already mapped at this locus (e.g. by CARD before VFDB) wins
            feature_location = (int(feature.location.start), int(feature.location.end))
            if feature.annotations.get("nucl_mapping") and feature.location.strand == region.strand \
                    and same_locus(feature_location, location):
                if feature_location == location:
                    self.annotate(feature, group, secondary=True)
                else:
                    self.logger.debug(f"Mapping of {primary.ref_id} at {location} overlaps {feature.id}; skipped.")
                return

        classes = {cds.id: classify_cds(cds, region, self.overlap_tolerance) for cds in record.features
                   if cds.type == "CDS" and not cds.annotations.get("nucl_mapping")}
        same = next((f for f in record.features if classes.get(f.id) == "same_gene"), None)
        defect = region.kind != "valid" and not (region.kind == "no_start" and same is not None)
        distant = defect and primary.pident < self.pseudogene_pident_cutoff  # too distant to call a pseudogene
        if not valid or distant:  # annotate only
            if valid:
                self.logger.debug(f"{region.kind} in mapping of {primary.ref_id} ({primary.pident}%) at {location} "
                                  f"not called as a pseudogene: identity below {self.pseudogene_pident_cutoff}%")
            if same is not None:
                self.annotate(same, group)
            return

        feature = None
        removed = {cds_id for cds_id, c in classes.items() if c == "overlapping" and cds_id not in self.same_genes}
        if region.kind == "valid" and same is None:
            action, feature = "rescued", new_cds(region, self.new_id(), self.transl_table)
            record.features.append(feature)
        elif region.kind == "valid" and (int(same.location.start), int(same.location.end)) == location:
            self.annotate(same, group, notes=region.notes)
        elif region.kind == "valid":
            ref_label = f"{self.db_name}:{self.references[primary.ref_id].to_dict()['gene']}"
            self.logger.debug(f"{same.id} relocated from {same.location} to {location} by mapping of {ref_label}")
            action, feature = "start corrected", same
            relocate_cds(same, region, ref_label)
        elif region.kind == "no_start" and same is not None:  # the reference start is not usable
            self.annotate(same, group, notes=[n for n in region.notes if n != "no start codon"])
        else:
            action, feature = PSEUDOGENE_ACTIONS[region.kind], new_pseudogene(region, self.new_id())
            feature.secondary_hits.append(PseudoGene(primary.ref_id, region.stop_codons, region.indels))
            record.features.append(feature)
            if same is not None:
                removed.add(same.id)  # replaced by the pseudogene
        self.remove_cds(record, removed)
        if feature is not None:
            self.set_hits(feature, group, action)

    def new_id(self):
        return f"MAP_{next(self.map_ids)}"

    def nucl_hits(self, group):
        """NuclHits of a group: the primary hit first, then the others flagged as alternative hits."""
        return [NuclHit(h.ref_id, self.references[h.ref_id], self.db_name, h.evalue, h.bitscore, h.pident,
                        100.0, 100.0, "alternative hit" if i else "") for i, h in enumerate(group)]

    def annotate(self, feature, group, secondary=False, notes=()):
        """Attach the hits to a feature without changing it, as NuclSearch does."""
        hit, *others = self.nucl_hits(group)
        feature.qualifiers.setdefault("note", []).extend(notes)
        if secondary or feature.primary_hit:
            feature.secondary_hits.append(hit)
        else:
            feature.primary_hit = hit
        feature.secondary_hits.extend(others)
        self.hits_json[feature.id] = hit.to_dict()
        self.counts["annotated"] += 1

    def set_hits(self, feature, group, action):
        """Make the mapping the primary hit of a feature it created or changed."""
        hit, *others = self.nucl_hits(group)
        feature.primary_hit = hit
        feature.secondary_hits.extend(others)
        d = hit.to_dict()
        d["model"]["note"] = action + "; " + d["model"]["note"]
        self.hits_json[feature.id] = d
        self.counts[action.split(":")[0]] += 1
        if action != "start corrected":
            self.logger.debug(f"{feature.id} ({action}) added at {feature.location} by mapping of {hit.id}")

    def remove_cds(self, record, cds_ids):
        for feature in record.features:
            if feature.id in cds_ids:
                self.logger.debug(f"{feature.id} at {feature.location} removed: overlaps a mapped gene")
        record.features[:] = [f for f in record.features if f.id not in cds_ids]
        self.counts["removed"] += len(cds_ids)

    def write_hits(self):
        hits_json_file = os.path.join(self.workDir, "nucl_hits.json")
        with open(hits_json_file, "w") as f:
            json.dump(self.hits_json, f, indent=4)

    def run(self):
        db_prefix = self.create_db()
        self.executeCommand(self.blast_command(db_prefix), shell=True)
        self.apply(self.select_hits())
        self.write_hits()
