#! /usr/bin/env python
# coding: UTF8

import os
from .base_tools import ContigAnnotationTool
# from Bio import SeqIO
# from Bio.SeqFeature import SeqFeature, FeatureLocation
from ..models.bio_feature import ExtendedFeature
import json

class Plasmidfinder(ContigAnnotationTool):
    """
    Plasmidfinder

    Tool type: contig annotation
    URL: https://github.com/tseemann/barrnap
    REF:

    """
    version = None
    TYPE = "source"
    NAME = "Plasmidfinder"
    VERSION_CHECK_CMD = None
    VERSION_PATTERN = None
    # SHELL = True


    def __init__(self, options=None, workDir="OUT"):
        """
        """
        if options is None:
            options = {}
        super(Plasmidfinder, self).__init__(options, workDir)
        self.cmd_options = options.get("cmd_options", "")
        self.db_path = options.get("db_path", "")
        self.output_directory = os.path.join(self.workDir, "contig_annotation", "plasmidfinder")
        self.result_file = os.path.join(self.output_directory, "data.json")
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

    def setVersion(self):
        """
        Plasmidfinder has no version number.
        """
        version = "-"
        self.__class__.version = version
        self.logger.info("{self.__class__.NAME} initialized. (Version {self.__class__.version})".format(self=self))


    def getCommand(self):
        """plasmidfinder.py -i testseq_AF250878.fasta -o OUT/ -tmp tmpOUT/ -p ${DB_ROOT}/plasmidfinder_db"""
        cmd = ["plasmidfinder.py", "--infile", self.genomeFasta, "--outputPath", self.output_directory, "--tmp_dir", self.output_directory, "--databasePath", self.db_path]
        return cmd

    def getResult(self):
        """to be implemented"""
        D = {}  # result
        data = json.load(open(self.result_file))
        pf_results = data.get("plasmidfinder", {})
        for organism, results in pf_results.get("results", {}).items():
            # print(organism)
            for database, hits in results.items():
                # print(database)
                if isinstance(hits, dict):
                    for hit_id, hit in hits.items():
                        identity = str(hit["identity"])
                        coverage = f'{hit["HSP_length"]}/{hit["template_length"]}'
                        plasmid_type = hit["plasmid"]
                        contig_name = hit["contig_name"]
                        accession = hit["accession"]
                        positions_in_contig = f'{contig_name}:{hit["positions_in_contig"]}'
                        # print(hit_id)
                        # example
                        # Possibly derived from Enterobacteriales IncFIA(HI1) type plasmid.
                        # Similar to AF250878, sequence1:157357..157744, identity:100.0%, aligned:388/388
                        description = f"Possibly derived from {organism} {plasmid_type} type plasmid."
                        hit_info = f"Similar to {accession}, {positions_in_contig}, identity:{identity}%, aligned:{coverage}"
                        # print("\t".join([description, hit_info]))

                        list_note = D.setdefault(contig_name, [])
                        list_note.append(f"PlasmidFinder: {plasmid_type}")
                        list_note.append(description)
                        list_note.append(hit_info)
        # print(D)
        return D

