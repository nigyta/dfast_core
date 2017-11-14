#! /usr/bin/env python
# coding: UTF8

import re
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from .base_tools import StructuralAnnotationTool
from ..models.bio_feature import ExtendedFeature


class GAP(StructuralAnnotationTool):
    """
    GAP annotator
    Find assembly gaps in the genome
    This class does not use external commands.
    It runs in the same way as other Tool-classes

    Tool type: assembly_gap
    URL:
    REF:

    """
    version = None
    TYPE = "assembly_gap"
    NAME = "GAPannotator"
    VERSION_CHECK_CMD = None
    VERSION_PATTERN = None

    def __init__(self, options=None, workDir="OUT"):
        super(GAP, self).__init__(options, workDir)
        self.gap_type = options.get("gap_type", "within scaffold")
        self.linkage_evidence = options.get("linkage_evidence", "paired-ends")
        self.len_cutoff = options.get("len_cutoff", 5)
        if self.len_cutoff < 1:
            self.logger.error("len_cutoff ({}) must be larger than 0. Aborted...".format(self.len_cutoff))
            exit(1)

    def setVersion(self):
        """
        GAP is not an external program.
        """
        version = "1.0"
        self.__class__.version = version
        self.logger.info("{self.__class__.NAME} initialized. (Version {self.__class__.version})".format(self=self))

    def run(self):
        def _generateGapFeature(genomeFasta):
            i = 0
            for record in SeqIO.parse(open(self.genomeFasta), "fasta", IUPAC.ambiguous_dna):
                startPosition = 0
                seq = str(record.seq).upper()
                pat = "(" + "N" * self.len_cutoff + "+)"
                for fragment in re.split(pat, seq):
                    endPosition = startPosition + len(fragment)
                    if fragment.startswith("N"):
                        i += 1
                        # qualifiers = {"estimated_length": [len(fragment)], "gap_type": [self.gap_type],
                        #               "linkage_evidence": [self.linkage_evidence]}
                        qualifiers = {"estimated_length": ["known"], "gap_type": [self.gap_type],
                                      "linkage_evidence": [self.linkage_evidence]}
                        location = FeatureLocation(startPosition, endPosition, strand=1)
                        feature = ExtendedFeature(location, type="assembly_gap", qualifiers=qualifiers,
                                                  id="{0}_{1}".format(self.__class__.__name__, i), seq_id=record.id)

                        assert str(feature.extract(record).seq).upper() == fragment
                        yield record.id, feature
                    startPosition = endPosition

        self.gapFeatures = {}
        for seqID, gapFeature in _generateGapFeature(self.genomeFasta):
            self.gapFeatures.setdefault(seqID, []).append(gapFeature)

    def getFeatures(self):
        return self.gapFeatures


