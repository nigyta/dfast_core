#! /usr/bin/env python
# coding: UTF8


import os
from .baseComponent import BaseAnnotationComponent
from ..tools.hmmer import Hmmer_hmmsearch


class DnaAfinder(BaseAnnotationComponent):
    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        CPU = 1
        super(DnaAfinder, self).__init__(genome, options, workDir, CPU)
        self.hmmer = Hmmer_hmmsearch(options=options)
        self.hmm_profile = options.get("hmm_profile", "")
        self.offset = int(options.get("offset", 0))

    def createCommands(self):
        result_file = os.path.join(self.workDir, "result0.out")
        cmd = self.hmmer.get_command(self.query_files[0], self.hmm_profile, result_file)
        self.commands.append(cmd)

    def getBestHit(self, fileName):
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
            query_id, evalue, score = field[0], field[4], field[5]
            return query_id, evalue, score
        else:
            return None, None, None

    def set_results(self):
        def _rotate(seq_record, start, strand, offset=0):
            seq = seq_record.seq
            length = len(seq)
            if strand == 1:
                origin = start - offset
                rotated = seq[origin:] + seq[:origin]
                self.logger.warn("'{0}' will be flipped/rotated so that the base position {1} comes first.".format(seq_record.id, origin + 1))
            else:
                origin = start + offset
                if origin > length:
                    origin = origin % length
                rotated = seq[origin:] + seq[:origin]
                rotated = rotated.reverse_complement()
                self.logger.warn("'{0}' will be flipped/rotated so that the base position {1} comes first.".format(seq_record.id, origin))
            return rotated

        self.logger.info("Checking HMMsearch results.")
        result_file = os.path.join(self.workDir, "result0.out")
        query_id, evalue, score = self.getBestHit(result_file)
        if query_id is not None:
            dnaA_gene = self.genome.features[query_id]
            # print(query_id, evalue, score)
            # print(dnaA_gene.seq_id)
            strand = dnaA_gene.location.strand
            if strand == 1:
                start = int(dnaA_gene.location.start)
                strand_symbol = "+" if strand == 1 else "-"
                self.logger.warn("DnaA gene was found at {4}:{0}[{1}]. (E-val: {2}, score: {3})".format(
                    start + 1, strand_symbol, evalue, score, dnaA_gene.seq_id))
            else:
                start = int(dnaA_gene.location.end)
                strand_symbol = "-" if strand == 1 else "-"

                self.logger.warn("DnaA gene was found at {4}:{0}[{1}]. (E-val: {2}, score: {3})".format(
                    start, strand_symbol, evalue, score, dnaA_gene.seq_id))
            seq_record = self.genome.seq_records[dnaA_gene.seq_id]
            seq_record.seq = _rotate(seq_record, start, strand, self.offset)
        else:
            self.logger.warn("DnaA gene was not found. The sequence file will be generated as is.")
            # print(seq_record)
            # print(seq_record.seq)


    def run(self):
        self.prepareQueries()
        self.createCommands()
        self.executeCommands(shell=False)
        self.set_results()

if __name__ == '__main__':
    pass
