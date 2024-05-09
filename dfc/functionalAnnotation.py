#! /usr/bin/env python
# coding: UTF8

import os
from logging import getLogger
from concurrent import futures

from .components.baseComponent import BaseAnnotationComponent
from .components.DBsearch import DBsearch
from .components.OrthoSearch import OrthoSearch
from .components.CDDsearch import CDDsearch
from .components.HMMscan import HMMscan
from .components.PseudoGeneDetection import PseudoGeneDetection
from .components.BlastSearch import BlastSearch
from .components.MBGDsearch import MBGDsearch
from .components.DnaAfinder import DnaAfinder
from .components.NuclSearch import NuclSearch

COMPONENTS = {
    "BASE": BaseAnnotationComponent,
    "DBsearch": DBsearch,
    "OrthoSearch": OrthoSearch,
    "CDDsearch": CDDsearch,
    "HMMscan": HMMscan,
    "PseudoGeneDetection": PseudoGeneDetection,
    "BlastSearch": BlastSearch,
    "MBGDsearch": MBGDsearch,
    "DnaAfinder": DnaAfinder,
    "NuclSearch": NuclSearch,
}


class FunctionalAnnotation(object):
    def __init__(self, genome, config):
        self.genome = genome
        self.configs = config.FUNCTIONAL_ANNOTATION
        self.workDir = config.WORK_DIR
        self.CPU = config.CPU
        self.logger = getLogger(__name__)

        self.components = []
        self.logger.info("Initializing annotation components... ")

        for component_config in self.configs:
            Component = COMPONENTS.get(component_config["component_name"])
            if Component:
                if component_config.get("enabled"):
                    self.logger.info("{Component.__name__} initialized.".format(Component=Component))
                    component = Component(self.genome, options=component_config.get("options", {}), workDir=self.workDir, CPU=self.CPU)
                    if isinstance(component, PseudoGeneDetection):
                        component.set_components(self.get_components())
                    self.components.append(component)
            else:
                self.logger.error("{0} is not registered in {1}. \nProcess aborted due to an error.".format(component_config["component_name"], __name__))
                exit(1)

    def get_components(self):
        '''
        Currently, this method is only for PseudoGeneDetection.
        Returns a list that contains all reference data
        :return:
        '''
        return [component for component in self.components if hasattr(component, "references")]


    def execute(self):
        """
        Run functional annotation tools parallelly using multi threading
        :return:
        """
        self.logger.info("Start executing functional annotation components.")
        for component in self.components:
            component.run()

    def cleanup(self):
        for component in self.components:
            component.cleanup()