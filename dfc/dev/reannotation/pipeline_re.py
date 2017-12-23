#! /usr/bin/env python
# coding: UTF8

import os
import shutil
from datetime import datetime
from logging import FileHandler, Formatter
from dfc.utils.path_util import create_output_directory
from dfc.functionalAnnotation import FunctionalAnnotation
from .genome_re import GenomeReannotation
from dfc.structuralAnnotation import StructuralAnnotation
from dfc.utils.feature_util import FeatureUtil
from dfc.utils.locus_tag_generator import LocusTagGenerator
from dfc.utils.format_converter import write_results
from dfc.utils.ddbj_submission import DDBJsubmission
from dfc.utils.genbank_submission import GenBankSubmission
from dfc.utils.genome_stat import GenomeStat
from dfc import dfast_version


class Pipeline():
    def __init__(self, config, logger):

        self.config = config
        self.logger = logger
        self.start_time = datetime.now()
        # set/create output directory
        create_output_directory(config)
        log_file = os.path.join(config.WORK_DIR, "application.log")
        log_file_handler = FileHandler(log_file, mode="w")
        if config.DEBUG:
            log_file_handler.setFormatter(Formatter("%(asctime)s %(levelname)8s %(message)s [%(module)s %(lineno)d]",
                                                    datefmt='%Y/%m/%d %H:%M:%S'))
        else:
            log_file_handler.setFormatter(Formatter("%(asctime)s %(message)s", datefmt='%Y/%m/%d %H:%M:%S'))
        self.logger.addHandler(log_file_handler)

        self.logger.info("DFAST-REANNOTATION pipeline started. (version {})".format(dfast_version))
        self.logger.warn("[WARNING] Results of the pipeline are not guaranteed as this is a prototype implementation.")
        self.logger.info("Results will be generated into '{}'.".format(config.WORK_DIR))

        # initializing query genome.
        self.genome = GenomeReannotation(config)

        self.ltg = LocusTagGenerator(self.genome, config)
        self.fu = FeatureUtil(self.genome, config)
        self.sa = StructuralAnnotation(self.genome, config)
        self.fa = FunctionalAnnotation(self.genome, config)
        self.ddbj = DDBJsubmission(self.genome, config)
        self.genbank = GenBankSubmission(self.genome, config)

    def execute(self):
        # self.sa.execute()  # execute structural annotation
        # self.fu.execute()  # feature adjustment: sort, remove_partial, (merge)
        self.fa.execute()  # functional annotation
        self.ltg.execute()  # assigning locus_tags
        self.genome.add_source_features()  # set source feature

        # writing result files.
        write_results(self.genome, self.config)
        GenomeStat.execute(self.genome)
        
        self.logger.warn("INSDC submission files will be created, although the file format may be invalid.")
        self.ddbj.create_submission_file()
        self.genbank.create_submission_file()

        if self.config.DEBUG:
            self.genome.to_pickle(os.path.join(self.config.WORK_DIR, "genome.pickle"))
        
        end_time = datetime.now()
        running_time = end_time - self.start_time
        running_time = running_time.total_seconds()
        h, remainder = divmod(running_time, 3600)
        m, s = divmod(remainder, 60)
        self.logger.info("DFAST pipeline completed!")
        self.logger.info("Total running time: {0:.0f}h{1:.0f}m{2:.0f}s".format(h, m, s))

    def cleanup(self):
        self.sa.cleanup()
        self.fa.cleanup()
        shutil.rmtree(os.path.join(self.config.WORK_DIR, "input"))