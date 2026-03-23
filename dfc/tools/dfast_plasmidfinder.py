#! /usr/bin/env python
# coding: UTF8

import os
from .base_tools import ContigAnnotationTool
from ..models.bio_feature import ExtendedFeature
import json

class Plasmidfinder(ContigAnnotationTool):
    """
    Plasmidfinder

    Tool type: contig annotation
    URL: https://bitbucket.org/genomicepidemiology/plasmidfinder
    REF: Carattoli et al., 2014

    """
    version = None
    TYPE = "source"
    NAME = "Plasmidfinder"
    VERSION_CHECK_CMD = ["plasmidfinder.py", "-v"]
    VERSION_PATTERN = r"(\d+\.\d+\.\d+)"


    def __init__(self, options=None, workDir="OUT"):
        if options is None:
            options = {}
        super(Plasmidfinder, self).__init__(options, workDir)
        self.cmd_options = options.get("cmd_options", "")
        self.db_path = options.get("db_path", "")
        self.method_path = options.get("method_path", "")
        self.output_directory = os.path.join(self.workDir, "contig_annotation", "plasmidfinder")
        self.result_file = os.path.join(self.output_directory, "results.json")
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

    def getCommand(self):
        cmd = ["plasmidfinder.py", "--infile", self.genomeFasta,
               "--outputPath", self.output_directory,
               "--tmp_dir", self.output_directory,
               "--databasePath", self.db_path,
               "-j", self.result_file]
        if self.method_path:
            cmd += ["-mp", self.method_path]
        return cmd

    def getResult(self):
        D = {}  # result: key=contig_name, value=list of note strings
        report = {}  # for AMR report
        data = json.load(open(self.result_file))
        seq_regions = data.get("seq_regions", {})
        for key, hit in seq_regions.items():
            plasmid_type = hit["name"]
            identity = str(hit["identity"])
            alignment_length = hit["alignment_length"]
            ref_length = hit.get("ref_seq_lenght", hit.get("ref_seq_length", 0))
            coverage = f"{alignment_length}/{ref_length}"
            accession = hit["ref_acc"]
            contig_name = hit["query_id"]
            start_pos = hit["query_start_pos"]
            end_pos = hit["query_end_pos"]
            positions_in_contig = f"{contig_name}:{start_pos}..{end_pos}"

            # Extract organism from ref_database (e.g. "PlasmidFinder-2.2.0:enterobacteriales")
            ref_database = hit.get("ref_database", "")
            organism = ref_database.split(":")[-1] if ":" in ref_database else ""

            description = f"Possibly derived from {organism} {plasmid_type} type plasmid."
            hit_info = f"Similar to {accession}, {positions_in_contig}, identity:{identity}%, aligned:{coverage}"

            list_note = D.setdefault(contig_name, [])
            list_note.append(f"PlasmidFinder: {plasmid_type}")
            list_note.append(description)
            list_note.append(hit_info)
            report.setdefault(contig_name, []).append(
                f"PlasmidFinder: {plasmid_type}, Description: {description}, HitInfo: {hit_info}")
        return D, report

