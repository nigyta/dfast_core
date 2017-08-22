#! /usr/bin/env python
# coding: UTF8

from logging import getLogger
import math


class LocusTagGenerator(object):
    def __init__(self, genome, config):
        self.logger = getLogger(__name__)
        self.genome = genome
        self.locus_tag_prefix = config.LOCUS_TAG_SETTINGS.get("locus_tag_prefix", "LOCUS")
        self.step = config.LOCUS_TAG_SETTINGS.get("step", 1)
        self.separate = config.LOCUS_TAG_SETTINGS.get("use_separate_tags", True)
        self.symbols = config.LOCUS_TAG_SETTINGS.get("symbols", {"CDS": "", "rRNA": "r", "tRNA": "t"})
        if self.locus_tag_prefix:
            self.enabled = True
            self.logger.info("Locus_tag settings: locus_tag_prefix={self.locus_tag_prefix} and step={self.step}.".format(self=self))
            if self.separate:
                examples = ", ".join([key + ": " + self.locus_tag_prefix + "_" + value + "000xx" for key, value in self.symbols.items()])
                self.logger.info("Locus_tags are assigned separately to each feature type. e.g. " + examples + ".")
        else:
            self.enabled = False
            self.logger.info("locus_tag_prefix is empty. Locus_tags will not be assigned.")

    def execute(self):
        if not self.enabled:
            return
        features = ", ".join(list(self.symbols.keys()))
        self.logger.info("Assigning locus_tags to " + features)

        counts = {key: 0 for key in self.symbols.keys()}
        count = 0
        digit = 1 if len(self.genome.features) == 0 else int(math.log10(len(self.genome.features) * self.step)) + 1
        for feature in self.genome.features.values():
            type_ = feature.type
            if type_ in self.symbols:
                counts[type_] += 1
                count += 1
                if self.separate:
                    locus_tag = self.locus_tag_prefix + "_" + self.symbols[type_] + str(self.step * counts[type_]).zfill(digit)
                else:
                    locus_tag = self.locus_tag_prefix + "_" + str(self.step * count).zfill(digit)
                feature.qualifiers["locus_tag"] = [locus_tag]
            else:
                if "locus_tag" in feature.qualifiers:
                    del feature.qualifiers["locus_tag"]

if __name__ == '__main__':
    pass
