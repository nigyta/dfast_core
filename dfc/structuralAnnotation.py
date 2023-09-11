#! /usr/bin/env python
# coding: UTF8

import os
import shutil
from logging import getLogger
from concurrent import futures
from Bio.SeqFeature import ExactPosition, FeatureLocation
from Bio.Data.CodonTable import TranslationError
from .tools.mga import MGA
from .tools.barrnap import Barrnap
from .tools.aragorn import Aragorn
from .tools.gap import GAP
from .tools.CRT import CRT
# from .tools.glimmer import Glimmer
from .tools.prodigal import Prodigal
from .tools.genemarkS2 import GeneMarkS2
from .tools.tRNAscan import tRNAscan
from .tools.rnammer import RNAmmer
from .tools.gff_importer import GFFimporter

TOOLS = {
    "GAP": GAP,
    "MGA": MGA,
    "Barrnap": Barrnap,
    "Aragorn": Aragorn,
    "CRT": CRT,
    # "Glimmer": Glimmer,
    "Prodigal": Prodigal,
    "GeneMarkS2": GeneMarkS2,
    "tRNAscan": tRNAscan,
    "RNAmmer": RNAmmer,
    "GFF_import": GFFimporter,
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
                exit(1)

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
        
        def _get_translation(feature, seq):
            nucseq = feature.location.extract(seq)
            offset = feature.qualifiers.get("codon_start", [1])[0] - 1
            right_offset = -1 * ((len(nucseq) - offset) % 3)
            if hasattr(tool, "transl_table"):
                transl_table = tool.transl_table
            else:
                transl_table = feature.qualifiers.get("transl_table", [11])[0]
            if transl_table == 4:
                start_codons = ["TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "GTG"]  # and ATG  for transl_table 4
            else:
                start_codons = ["TTG", "CTG", "ATT", "ATC", "ATA", "GTG"]  # and ATG   for transl_table 11

            if right_offset == 0:
                if isinstance(feature.location.start, ExactPosition) and isinstance(feature.location.end,
                                                                                    ExactPosition):
                    try:
                        translation = nucseq[offset:].translate(table=transl_table, cds=True)
                    except TranslationError:
                        translation = nucseq[offset:].translate(table=transl_table, to_stop=True)
                        if len(translation) * 3 != len(nucseq[offset:]):
                            self.logger.warning(
                                "Translation error in {}. In-frame stop codon exists. Translation was terminated at the first in-frame stop codon.".format(
                                    feature.id))
                            before = str(feature.location)
                            if feature.location.strand == 1:
                                start = feature.location.start
                                end = ExactPosition(start + offset + len(translation) * 3 + 3)  # Trailing +3 for stop-codon
                                feature.location = FeatureLocation(start, end, 1)
                            else:
                                end = feature.location.end
                                start = ExactPosition(end - offset - len(translation) * 3 - 3)  # Trailing -3 for stop-codon
                                feature.location = FeatureLocation(start, end, -1)
                            after = str(feature.location)
                            self.logger.warning("CDS[{}] was fixed from {} to {}.".format(feature.id, before, after))
                else:
                    translation = nucseq[offset:].translate(table=transl_table, to_stop=True)
            else:
                translation = nucseq[offset:right_offset].translate(table=transl_table)  # , stop_symbol="")
            translation = str(translation)
            first_codon = str(nucseq[offset:offset + 3]).upper()
            if first_codon in start_codons:
                translation = "M" + translation[1:]
            return translation

        for tool in self.tools:
            dict_features = tool.getFeatures()

            detected_features = 0
            for seqID in dict_features.keys():
                record = self.genome.seq_records[seqID]
                features = dict_features[seqID]
                detected_features += len(features)
                if tool.TYPE == "CDS":
                    for feature in features:
                        translation = _get_translation(feature, record.seq)
                        # nucseq = feature.location.extract(record.seq)
                        # offset = feature.qualifiers.get("codon_start", [1])[0] - 1
                        # right_offset = -1 * ((len(nucseq) - offset) % 3)
                        # if right_offset == 0:
                        #     if isinstance(feature.location.start, ExactPosition) and isinstance(feature.location.end, ExactPosition):
                        #         try:
                        #             translation = nucseq[offset:].translate(table=tool.transl_table, cds=True)
                        #         except TranslationError:
                        #             translation = nucseq[offset:].translate(table=tool.transl_table, to_stop=True)
                        #             if len(translation) * 3 != len(nucseq[offset:]):
                        #                 self.logger.warning("Translation error. In-frame stop codon exists. Translation was terminated at the first in-frame stop codon. [{}]".format(feature.id))
                        #                 before = str(feature.location)
                        #                 if feature.location.strand == 1:
                        #                     start = feature.location.start
                        #                     end = ExactPosition(start + offset + len(translation) * 3)
                        #                     feature.location = FeatureLocation(start, end, 1)
                        #                 else:
                        #                     end = feature.location.end
                        #                     start = ExactPosition(end - offset - len(translation) * 3)
                        #                     feature.location = FeatureLocation(start, end, -1)
                        #                 after = str(feature.location)
                        #                 self.logger.warning("CDS is fixed from {} to {}. [{}]".format(before, after, feature.id))
                        #     else:
                        #         translation = nucseq[offset:].translate(table=tool.transl_table, to_stop=True)
                        # else:
                        #     translation = nucseq[offset:right_offset].translate(table=tool.transl_table)  # , stop_symbol="")
                        feature.qualifiers["translation"] = [translation]
                elif tool.TYPE == "GFF":  # TODO: preliminary implementation (redundant implementation of the code above) 
                    for feature in features:
                        if feature.type == "CDS":
                            translation = _get_translation(feature, record.seq)
                            # nucseq = feature.location.extract(record.seq)
                            # offset = feature.qualifiers.get("codon_start", [1])[0] - 1
                            # transl_table = feature.qualifiers.get("transl_table", [11])[0]
                            # right_offset = -1 * ((len(nucseq) - offset) % 3)
                            # if right_offset == 0:
                            #     try:
                            #         translation = nucseq[offset:].translate(table=transl_table, cds=True)
                            #     except TranslationError as e:
                            #         print(e)
                            #         self.logger.warning("Warning: Translation error. In-frame stop codon may exist. {}".format(feature.id))
                            #         translation = nucseq[offset:].translate(table=transl_table, to_stop=True)
                            # else:
                            #     translation = nucseq[offset:right_offset].translate(table=transl_table)  # , stop_symbol="")
                            feature.qualifiers["translation"] = [translation]
                record.features += features
            self.logger.info("{num} {tool.TYPE} features were detected by {tool.__class__.NAME}.".format(
                num=detected_features, tool=tool))
        self.genome.sort_features()
        self.genome.set_feature_dictionary()

    def cleanup(self):
        shutil.rmtree(self.child_dir)
