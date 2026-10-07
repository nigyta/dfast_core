#! /usr/bin/env python
# coding: UTF8

import csv
import glob
import json
import os
from .baseComponent import BaseAnnotationComponent
from ..models.hit import DefenseHit
from ..tools.defense_finder import DefenseFinder as DefenseFinderTool


def read_genes_file(file_name):
    """Read *_defense_finder_genes.tsv. Returns a list of dicts (one per gene in a system)."""
    with open(file_name) as f:
        return list(csv.DictReader(f, delimiter="\t"))


class DefenseFinder(BaseAnnotationComponent):
    """Detects anti-phage defense systems (and anti-defense systems with '-a') among the CDSs.

    Proteins are given in the gembase format, '<replicon>_<number>', so that a system never spans two sequences.
    Replicon names are generated (R0001, R0002, ...) because macsyfinder trims trailing 'i'/'b' from them.
    """
    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        super(DefenseFinder, self).__init__(genome, options, workDir, CPU)
        self.tool = DefenseFinderTool()
        self.models_dir = options.get("models_dir", "")
        self.cmd_options = options.get("cmd_options", "")
        if not os.path.exists(self.models_dir):
            self.logger.error("DefenseFinder models are not found in '{0}'. "
                              "Run 'dfast_file_downloader.py --defensefinder' to download them.".format(self.models_dir))
            exit(1)
        self.query_file = os.path.join(self.workDir, "query.faa")
        self.query_ids = {}  # gembase id -> feature id

    def prepare_query(self):
        with open(self.query_file, "w") as f:
            for i, record in enumerate(self.genome.seq_records.values(), 1):
                cds_features = [feature for feature in record.features if feature.type == "CDS"]
                for j, feature in enumerate(cds_features, 1):
                    gembase_id = "R{0:04d}_{1:05d}".format(i, j)
                    self.query_ids[gembase_id] = feature.id
                    f.write(">{0}\n{1}\n".format(gembase_id, feature.qualifiers.get("translation", [""])[0]))

    def set_results(self):
        systems = {}
        for genes_file in glob.glob(os.path.join(self.workDir, "*_defense_finder_genes.tsv")):
            for row in read_genes_file(genes_file):
                feature_id = self.query_ids.get(row["hit_id"])
                if feature_id is None:
                    continue
                hit = DefenseHit(row["type"], row["subtype"], row["activity"], row["gene_name"], row["sys_id"])
                self.genome.features[feature_id].secondary_hits.append(hit)
                system = systems.setdefault(row["sys_id"], {
                    "type": row["type"], "subtype": row["subtype"], "activity": row["activity"], "features": []})
                system["features"].append(feature_id)
        # read by utils.summarize_amr
        with open(os.path.join(self.workDir, "systems.json"), "w") as f:
            json.dump(systems, f, indent=2)
        self.logger.info("DefenseFinder found {0} systems.".format(len(systems)))

    def run(self):
        self.prepare_query()
        self.logger.info("DefenseFinder will be performed using {0} CPUs.".format(self.CPU))
        cmd = self.tool.get_command(self.query_file, self.workDir, self.models_dir, self.CPU, self.cmd_options)
        self.executeCommand(cmd)
        self.set_results()
