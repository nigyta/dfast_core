#! /usr/bin/env python
# coding: UTF8

import re
from .base_tools import StructuralAnnotationTool
from ..models.bio_feature import ExtendedFeature


class GFFimporter(StructuralAnnotationTool):
    """
    Import features from GFF file
    This class does not use external commands.

    Tool type: GFF_import

    """
    version = None
    TYPE = "GFF"
    NAME = "GFFimporter"
    VERSION_CHECK_CMD = None
    VERSION_PATTERN = None

    def __init__(self, options=None, workDir="OUT"):
        super(GFFimporter, self).__init__(options, workDir)
        self.gff_file_name = options.get("gff_file_name")
        self.targets = options.get("targets", ["CDS"])
        self.imported_features = {}
        self.logger.info("Following features will be imported from GFF: " + ", ".join(self.targets))
    def setVersion(self):
        """
        GFFimporter is not an external program.
        """
        version = "1.0"
        self.__class__.version = version
        self.logger.info("{self.__class__.NAME} initialized. (Version {self.__class__.version})".format(self=self))

    def run(self):
        """Example of standard GFF format.
            sequence1	GeneMark.hmm	gene	337	2799	.	+	.	ID=1
            sequence1	GeneMark.hmm	CDS	337	2799	310.21	+	0	Parent=1; start_score=-0.30
            sequence1	GeneMark.hmm	gene	5683	6459	.	-	.	ID=5
            sequence1	GeneMark.hmm	CDS	5683	6459	96.29	-	0	Parent=5; start_score=2.63
        """

        def _parse_gff():
            i = 0
            with open(self.gff_file_name) as f:
                for line in f:
                    if line.startswith("##FASTA"):
                        break
                    elif line.startswith("#"):
                        continue
                    sequence_id, tool_name, feature_type, left, right, score, strand, _, qualifiers = line.strip("\n").split("\t")
                    tool_name = tool_name.replace(" ", "-")
                    if feature_type not in self.targets:
                        continue
                    i += 1
                    qualifiers = dict([x.split("=") for x in qualifiers.strip(";").split(";")])

                    gff_id = qualifiers.get("ID", "{0}_{1}".format(tool_name, i))
                    location = self.getLocation(left, right, strand)
                    feature = ExtendedFeature(location=location, type=feature_type, id=gff_id,
                                              seq_id=sequence_id, annotations={})

                    if feature_type == "CDS":
                        feature.qualifiers = {
                            # "product": [qualifiers.get("product", "hypothetical protein")],
                            "product": ["hypothetical protein"],
                            "transl_table": [qualifiers.get("transl_table", 11)],
                            "codon_start": [qualifiers.get("codon_start", 1)]
                        }
                    elif feature_type == "rRNA":
                        feature.qualifiers = {
                            "product": [qualifiers.get("product", "unknown rRNA")]
                        }
                    elif feature_type == "tRNA":
                        feature.qualifiers = {
                            "product": [qualifiers.get("product", "unknown tRNA")]
                        }
                    else:
                        feature.qualifiers = {
                            "product": [qualifiers.get("product", "unknown product")]
                        }
                    # if gff_id:
                    #     feature.qualifiers["note"] = ["gff_gene_id:" + gff_id]
                    yield feature

        for feature in _parse_gff():
            seq_id = feature.seq_id
            self.imported_features.setdefault(seq_id, []).append(feature)

    def getFeatures(self):
        return self.imported_features


