#! /usr/bin/env python
# coding: UTF8

import os
import subprocess
import shutil
from logging import getLogger
from concurrent import futures

from ..tools.ghostz import Ghostz
from ..tools.ghostx import Ghostx
from ..tools.blastp import Blastp
from ..tools.diamond import Diamond
from ..tools.blastn import Blastn

ALIGNERS = {
    "ghostz": Ghostz,
    "ghostx": Ghostx,
    "blastp": Blastp,
    "diamond": Diamond,
    "blastn": Blastn
 }

class BaseAnnotationComponent(object):

    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        self.logger = getLogger(__name__)
        self.__class__.instances += 1
        self.genome = genome
        self.options = options
        self.workDir = os.path.join(workDir, self.__class__.__name__)
        if self.__class__.instances > 1:
            self.workDir += "_" + str(self.__class__.instances)
        if not os.path.exists(self.workDir):
            os.makedirs(self.workDir)
        cpu = self.options.get("cpu")
        self.CPU = cpu if cpu else CPU  # if component-specific cpu number is not specified, global CPU number is used.  
        self.skipAnnotatedFeatures = self.options.get("skipAnnotatedFeatures", False)
        self.logger.info("Setting {0} options. {1}".format(self.__class__.__name__, str(self.options)))
        aligner = self.options.get("aligner")
        if aligner:
            if aligner not in ALIGNERS:
                self.logger.error("Aligner ({}) is not registered in {}".format(aligner, __file__))
                exit(1)
            aligner_option = self.options.get("aligner_option", {})
            self.aligner = ALIGNERS[aligner](aligner_option)
        self.query_files = {}
        self.query_sequences = {}
        self.commands = []

        
    def splitList(self, L, n):
        """
        Generic use function to split list
        """
        for i in range(n):
            yield L[i::n]

    def prepareQueries(self, type_="protein", split_query=True):
        """
        todo: implement other types (cds, extended_cds, ...)

        :param type_:
        :param split_query:
        :return:
        """
        cds_features = [feature for feature in self.genome.features.values() if feature.type == "CDS"]
        if self.skipAnnotatedFeatures:
            cds_features = [feature for feature in cds_features if feature.primary_hit is None]

        split_num = self.CPU if split_query else 1
        i = 0
        for features in self.splitList(cds_features, split_num):
            fasta_buffer = ""
            for feature in features:
                ID = feature.id
                if type_ == "protein":
                    sequence = feature.qualifiers.get("translation", [""])[0]
                elif type_ == "nucleotide":  # get CDS
                    seq = self.genome.seq_records[feature.seq_id].seq
                    sequence = str(feature.extract(seq))
                else:
                    self.logger.error("Query type must be protein or nuclaotide.")
                    raise AssertionError

                fasta_buffer += ">{0}\n{1}\n".format(ID, sequence)
                self.query_sequences[ID] = sequence
            if not fasta_buffer:
                self.logger.debug("Query file is empty. Skip processing.")
                continue
            fasta_file = os.path.join(self.workDir, "query{0}.fasta".format(i))
            self.query_files[i] = fasta_file
            i += 1
            with open(fasta_file, "w") as f:
                f.write(fasta_buffer)

    @staticmethod
    def get_coverage(start, end, length):
        return 100.0 * (int(end) - int(start) + 1) / int(length)

    def executeCommand(self, cmd, shell=False):
        """
        Any output to standard error will be handled as an error.
        Some tools write logs to standard error. In that case, logs must be redirected to standard output.
        :param cmd:
        :param shell:
        :return:
        """
        self.logger.debug('Running command "{0}" ({1})'.format(" ".join(cmd), self.__class__.__name__))
        if shell:
            cmd = " ".join(cmd)
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, shell=shell)
        out, err = p.communicate()
        if p.returncode != 0 and err:
            self.logger.error("Process aborted due to an error in {self.__class__.__name__}.".format(self=self))
            self.logger.error(err.decode("utf8"))
            exit()
        return out, err

    def createCommands(self):
        """
        This is an example
        implement this

        """
        for i in range(self.CPU):
            cmd = ["echo", "query{0}=".format(i), self.query_files["query{0}".format(i)]]
            self.commands.append(cmd)

    def executeCommands(self, shell=False, process_name=None, verbose=True):
        if process_name is None:
            process_name = self.__class__.__name__
        if verbose:
            self.logger.info("{0} will be performed using {1} CPUs.".format(process_name, self.CPU))
        with futures.ThreadPoolExecutor(max_workers=self.CPU) as executor:
            mappings = {}
            num = len(self.commands)
            for i, command in enumerate(self.commands):

                mappings[executor.submit(self.executeCommand, command, shell)] = i + 1

            for future in futures.as_completed(mappings):
                finished_command = mappings[future]
                result = future.result()
                if verbose:
                    self.logger.info("{0} done {1}/{2}.".format(process_name, finished_command, num))
        self.commands = []

    def run(self):
        self.prepareQueries()
        self.createCommands()
        self.executeCommands(shell=True)

    def cleanup(self):
        shutil.rmtree(self.workDir)


if __name__ == '__main__':
    pass
