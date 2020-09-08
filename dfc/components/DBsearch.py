#! /usr/bin/env python
# coding: UTF8

import os
from .baseComponent import BaseAnnotationComponent
from ..models.protein import Protein
from ..models.hit import ProteinHit
from ..utils.reffile_util import check_db_file, read_db_attributes

class DBsearch(BaseAnnotationComponent):
    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        super(DBsearch, self).__init__(genome, options, workDir, CPU)
        self.evalue_cutoff = options.get("evalue_cutoff", 1e-5)
        self.qcov_cutoff = options.get("qcov_cutoff", 70)
        self.scov_cutoff = options.get("scov_cutoff", 70)
        self.pident_cutoff = options.get("pident_cutoff", 0)
        self.db_name = options.get("db_name", "")
        database = options.get("database", "")
        if database.endswith(".ref"):
            database = database[:-4]
        self.database = database
        aligner_name = self.options.get("aligner")
        check_db_file(self.database, aligner_name)

        self.references = self.read_reference()
        """
        todo: check if db file exists.
        """

    def read_reference(self):
        reference_file = self.database + ".ref"
        D = Protein.read_from_dfast_reference(reference_file)
        db_attributes = read_db_attributes(reference_file)
        if db_attributes:
            str_attributes = " [" + ", ".join([key + "=" + value for key, value in db_attributes.items()]) + "]"
        else:
            str_attributes = ""
        self.logger.info("Reference DB loaded. {0} sequences.{1}".format(len(D), str_attributes))
        return D
    

    def createCommands(self):
        for i, query in self.query_files.items():
            # example of ghostz ["ghostz", "aln", "-i", queryFileName, "-d", dbFileName, "-o", outFileName, "-b", "1"]
            # query_file = self.querie_files["query{0}".format(i)]
            result_file = os.path.join(self.workDir, "alignment{0}.out".format(i))

            cmd = self.aligner.get_command(query, self.database, result_file)
            self.commands.append(cmd)

    def parse_result(self, file_name):
        """
        to parse alignment result.
        By default, alignment result is tabular data.
        q_id, s_id, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore

        :param file_name:
        :return:
        """

        with open(file_name) as f:
            data = f.readlines()
            for row in data:
                q_id, s_id, pident, length, mismatch, gapopen, qstart, qend, \
                sstart, send, evalue, bitscore = row.strip("\n").split("\t")
                protein = self.references[s_id]
                evalue = float(evalue)
                if evalue > self.evalue_cutoff:
                    # filtering high e-value hit
                    continue
                qlen = len(self.query_sequences[q_id])
                slen = len(protein.sequence)
                q_cov = DBsearch.get_coverage(qstart, qend, qlen)
                s_cov = DBsearch.get_coverage(sstart, send, slen)
                # ProteinHit: id_, description, gene, ec_number, source_db, organism, db_name, e_value, score, identity, q_cov, s_cov, notes=[]
                flag = ""
                hit = ProteinHit(s_id, protein.description, protein.gene, protein.ec_number, protein.source_db,
                                 protein.organism, self.db_name, evalue, bitscore, pident, q_cov, s_cov,
                                 flag, notes=[])
                feature = self.genome.features[q_id]
                if hit.identity < self.pident_cutoff:
                    self.logger.debug("Low identity hit. {}-{} ({:.1f}%)".format(q_id, s_id, hit.identity))
                    hit.flag = "low identity"
                    feature.secondary_hits.append(hit)
                elif q_cov < self.qcov_cutoff or s_cov < self.scov_cutoff:
                    self.logger.debug("Partial hit. {}-{} (q_cov:{:.1f}%, s_cov{:.1f}%)".format(q_id, s_id, q_cov, s_cov))
                    hit.flag = "partial hit"
                    feature.secondary_hits.append(hit)
                elif feature.primary_hit:
                    feature.secondary_hits.append(hit)
                else:
                    feature.primary_hit = hit

    def set_results(self):
        for i in range(len(self.query_files)):
            result_file = os.path.join(self.workDir, "alignment{0}.out".format(i))
            self.parse_result(result_file)

    def run(self):
        # if not self.enabled:
        #     self.logger.warning("[Warning] {} is disabled. Skip execution.".format(self.__class__.__name__))
        self.prepareQueries()
        self.createCommands()
        self.executeCommands(shell=True)
        self.set_results()


if __name__ == '__main__':
    pass
