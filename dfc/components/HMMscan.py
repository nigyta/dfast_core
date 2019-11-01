#! /usr/bin/env python
# coding: UTF8


import os
from .baseComponent import BaseAnnotationComponent
from ..models.hit import HmmHit
from ..tools.hmmer import Hmmer_hmmscan
from ..utils.reffile_util import check_hmm_db_file

class HMMscan(BaseAnnotationComponent):
    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        super(HMMscan, self).__init__(genome, options, workDir, CPU)
        self.hmmer = Hmmer_hmmscan(options=options)
        self.database = options.get("database", "")
        self.db_name = options.get("db_name", "")
        check_hmm_db_file(self.database)

    def createCommands(self):
        for i, query in self.query_files.items():
            # example of ghostz ["ghostz", "aln", "-i", queryFileName, "-d", dbFileName, "-o", outFileName, "-b", "1"]
            result_file = os.path.join(self.workDir, "result{0}.out".format(i))
            cmd = self.hmmer.get_command(query, self.database, result_file)
            self.commands.append(cmd)

    def parseResult(self, fileName):
        '''
        File format example, tab-separated.
            #                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
            # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
            #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
        '''
        for line in open(fileName):
            if line.startswith("#"):
                continue
            field = line.strip("\n").split()
            query_id, accession, name, description, evalue, score, bias = field[2], field[1], field[0], " ".join(field[18:]), field[4], field[5], field[6]
            description = description.replace('"', "'")
            hmm = HmmHit(accession, name, description, evalue, score, bias, self.db_name)
            yield query_id, hmm

    def set_results(self):

        self.logger.info("Collecting HMMscan results.")
        hitDict = {}
        for i in range(len(self.query_files)):
            result_file = os.path.join(self.workDir, "result{0}.out".format(i))
            for query_id, hmm in self.parseResult(result_file):
                hitDict.setdefault(query_id, []).append(hmm)

        for query_id, hmms in hitDict.items():
            # sort by score in descending order. The best hit comes first.
            hmms.sort(key=lambda x: x.score, reverse=True)
            best_hmmhit = hmms[0]
            feature = self.genome.features[query_id]

            feature.secondary_hits.append(best_hmmhit)

    def run(self):
        self.prepareQueries()
        self.createCommands()
        self.executeCommands(shell=False)
        self.set_results()

if __name__ == '__main__':
    pass
