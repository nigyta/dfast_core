#! /usr/bin/env python
# coding: UTF8

from collections import Counter
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, BeforePosition, AfterPosition
from .base_tools import StructuralAnnotationTool
from ..models.bio_feature import ExtendedFeature

MIN_VERSION = (0, 1, 4)


class Diced(StructuralAnnotationTool):
    """
    Diced (CRISPR detection, Rust implementation of the MinCED/CRT approach, run in-process as a Python package)

    Tool type: CRISPR prediction
    URL: https://github.com/althonos/diced
    The default parameters are those of CRT, the tool Diced replaced in DFAST.
    """

    version = None
    TYPE = "CRISPR (later will be replaced by repeat_region)"
    NAME = "Diced"
    DEFAULT_OPTIONS = {"min_repeat_length": 19, "max_repeat_length": 38, "min_spacer_length": 19, "max_spacer_length": 48}

    def __init__(self, options=None, workDir="OUT"):
        super(Diced, self).__init__(options, workDir)
        self.params = {k: self.options.get(k, v) for k, v in self.DEFAULT_OPTIONS.items()}
        self.crisprs = {}  # {seq_id: [(start, end, most common repeat, seq length), ...]}, 0-based, end-exclusive

    def setVersion(self):
        self.logger.info("Checking {0} version... ".format(self.__class__.NAME))
        try:
            import diced
            version = diced.__version__
            ok = tuple(int(x) for x in version.split(".")[:3]) >= MIN_VERSION
        except (ImportError, ValueError):
            ok = False
        if not ok:
            self.logger.error("Diced >= {0} is required but was not found. Install it with "
                              "\"pip install 'diced>=0.1.4'\". Aborted...".format(".".join(map(str, MIN_VERSION))))
            exit(1)
        self.__class__.version = version
        self.logger.info("{0} initialized. (Version {1})".format(self.__class__.NAME, version))

    def run(self):
        import diced
        self.crisprs = {}
        for record in SeqIO.parse(self.genomeFasta, "fasta"):
            seq = str(record.seq).upper()
            for crispr in diced.scan(seq, **self.params):
                # Do not touch crispr.spacers; iterating it panics at the end in diced 0.1.4.
                repeat = Counter(str(x) for x in crispr.repeats).most_common(1)[0][0]
                self.crisprs.setdefault(record.id, []).append((crispr.start, crispr.end, repeat, len(seq)))

    def getFeatures(self):
        D = {}
        i = 0
        for seq_id, crisprs in self.crisprs.items():
            for start, end, repeat, length in crisprs:
                i += 1
                start_position = BeforePosition(0) if start == 0 else start
                end_position = AfterPosition(length) if end == length else end
                feature = ExtendedFeature(location=FeatureLocation(start_position, end_position, strand=1),
                                          type="CRISPR", seq_id=seq_id)
                feature.id = "{0}_{1}".format(self.__class__.__name__, i)
                feature.qualifiers["inference"] = \
                    ["COORDINATES:alignment:{0}:{1}".format(self.__class__.NAME, self.__class__.version)]
                feature.qualifiers["rpt_family"] = ["CRISPR"]
                feature.qualifiers["rpt_type"] = ["direct"]
                feature.qualifiers["rpt_unit_seq"] = [repeat.lower()]
                D.setdefault(seq_id, []).append(feature)
        return D
