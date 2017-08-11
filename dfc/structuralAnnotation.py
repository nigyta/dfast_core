#! /usr/bin/env python
# coding: UTF8

import os
import shutil
from logging import getLogger
from concurrent import futures
from Bio.SeqFeature import ExactPosition
from .tools.mga import MGA
from .tools.barrnap import Barrnap
from .tools.aragorn import Aragorn
from .tools.gap import GAP
from .tools.CRT import CRT
# from .tools.glimmer import Glimmer
from .tools.prodigal import Prodigal
from .tools.tRNAscan import tRNAscan
from .tools.rnammer import RNAmmer

TOOLS = {
    "GAP": GAP,
    "MGA": MGA,
    "Barrnap": Barrnap,
    "Aragorn": Aragorn,
    "CRT": CRT,
    # "Glimmer": Glimmer,
    "Prodigal": Prodigal,
    "tRNAscan": tRNAscan,
    "RNAmmer": RNAmmer,
}


class StructuralAnnotation(object):
    def __init__(self, genome, config):
        self.genome = genome
        self.configs = config.STRUCTURAL_ANNOTATION

        self.workDir = config.WORK_DIR
        self.child_dir = os.path.join(self.workDir, "structural")
        if not os.path.exists(self.child_dir):
            os.makedirs(self.child_dir)

        self.CPU = config.CPU
        self.logger = getLogger(__name__)

        self.tools = []
        self.logger.info("Initializing structural annotation tools... ")

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
                exit()

    def execute(self):
        """
        Run structural annotation tools in parallel using multi threading
        :return:
        """
        self.logger.info("Start structural annotation process using {self.CPU} CPUs".format(self=self))

        with futures.ThreadPoolExecutor(max_workers=self.CPU) as executor:
            mappings = {}
            for tool in self.tools:
                mappings[executor.submit(tool.run)] = tool.__class__.NAME

            for future in futures.as_completed(mappings):
                target = mappings[future]
                result = future.result()  # result will be None
                msg = '{target} finished.'.format(target=target)
                self.logger.debug(msg)
        self.setFeatures()

    def setFeatures(self):
        """
        Collect annotated features from each tool, and set them to the Genome object.
        :param genome:
        :return:
        """
        for tool in self.tools:
            dict_features = tool.getFeatures()

            detected_features = 0
            for seqID in dict_features.keys():
                record = self.genome.seq_records[seqID]
                features = dict_features[seqID]
                detected_features += len(features)
                if tool.TYPE == "CDS":
                    for feature in features:
                        nucseq = feature.location.extract(record.seq)
                        offset = feature.qualifiers.get("codon_start", [1])[0] - 1
                        rightOffset = -1 * ((len(nucseq) - offset) % 3)
                        if rightOffset == 0:
                            if isinstance(feature.location.start, ExactPosition) and isinstance(feature.location.end, ExactPosition):
                                translation = nucseq[offset:].translate(table=tool.transl_table, cds=True)
                            else:
                                translation = nucseq[offset:].translate(table=tool.transl_table, stop_symbol="")
                        else:
                            translation = nucseq[offset:rightOffset].translate(table=tool.transl_table, stop_symbol="")
                        feature.qualifiers["translation"] = [str(translation)]
                record.features += features
            self.logger.info("{num} {tool.TYPE} features were detected by {tool.__class__.NAME}.".format(
                num=detected_features, tool=tool))
        self.genome.sort_features()
        self.genome.set_feature_dictionary()

    def cleanup(self):
        shutil.rmtree(self.child_dir)
