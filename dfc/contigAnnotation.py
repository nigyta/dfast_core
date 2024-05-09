#! /usr/bin/env python
# coding: UTF8

import os
import shutil
from logging import getLogger
from concurrent import futures
from Bio.SeqFeature import ExactPosition, FeatureLocation
from .tools.dfast_plasmidfinder import Plasmidfinder

TOOLS = {
    "Plasmidfinder": Plasmidfinder
}


class ContigAnnotation(object):
    def __init__(self, genome, config):
        self.genome = genome

        # For compatibility for older config file without "CONTIG_ANNOTATION"
        if hasattr(config, "CONTIG_ANNOTATION"):
            self.configs = config.CONTIG_ANNOTATION
        else:
            self.configs = []

        self.workDir = config.WORK_DIR
        self.child_dir = os.path.join(self.workDir, "contig_annotation")
        if not os.path.exists(self.child_dir):
            os.makedirs(self.child_dir)

        self.CPU = config.CPU
        self.logger = getLogger(__name__)

        self.tools = []
        self.logger.info("Initializing contig annotation tools... ")

        for tool_config in self.configs:
            Tool = TOOLS.get(tool_config["tool_name"])
            if Tool:
                if tool_config.get("enabled"):
                    tool = Tool(options=tool_config.get("options", {}), workDir=self.workDir)
                    self.tools.append(tool)
            else:
                self.logger.error(
                    "{0} is not registered in {1}. \nProcess aborted due to an error.".format(tool_config["tool_name"],
                                                                                              __name__))
                exit(1)

    def execute(self):
        """
        Run contig annotation tools in parallel using multi threading
        :return:
        """
        self.logger.info("Start contig annotation process using {self.CPU} CPUs".format(self=self))

        with futures.ThreadPoolExecutor(max_workers=self.CPU) as executor:
            mappings = {}
            for tool in self.tools:
                mappings[executor.submit(tool.run)] = tool.__class__.NAME

            for future in futures.as_completed(mappings):
                target = mappings[future]
                result = future.result()  # result will be None
                msg = '{target} finished.'.format(target=target)
                self.logger.debug(msg)
        return self.collectResults()

    def collectResults(self):
        """
        Collect Results from each tool
        """
        dict_source_notes = {}
        dict_contig_annotation_report = {}
        for tool in self.tools:
            result, dict_report = tool.getResult()   # Dict: key=sequence_id, value=List of Note
            for seq_id, list_note in result.items():
                tmp_note_list = dict_source_notes.setdefault(seq_id, [])
                tmp_note_list += list_note
            for seq_id, report in dict_report.items():
                dict_contig_annotation_report.setdefault(seq_id, []).extend(report)
        return dict_source_notes, dict_contig_annotation_report

    def cleanup(self):
        shutil.rmtree(self.child_dir)
