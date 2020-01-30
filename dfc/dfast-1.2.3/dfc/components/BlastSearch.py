#! /usr/bin/env python
# coding: UTF8

import os
from .baseComponent import BaseAnnotationComponent
from ..models.hit import ProteinHit
from ..tools.blastdbcmd import Blastdbcmd

from ..utils.ref_util import fasta_parsers  # keys:  auto/ncbi/uniprot/plain
from ..models.protein import Protein


class BlastSearch(BaseAnnotationComponent):
    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        super(BlastSearch, self).__init__(genome, options, workDir, CPU)
        self.evalue_cutoff = options.get("evalue_cutoff", 1e-5)
        self.qcov_cutoff = options.get("qcov_cutoff", 70)
        self.scov_cutoff = options.get("scov_cutoff", 70)
        self.pident_cutoff = options.get("pident_cutoff", 0)

        self.database = options.get("database", "")
        self.db_name = options.get("db_name", "")
        aligner_name = options.get("aligner")
        assert aligner_name == "blastp"
        self.blastdbcmd = Blastdbcmd()
        # check_db_file(self.database, aligner_name)
        self.dbtype = options.get("dbtype")
        self.parser = fasta_parsers[self.dbtype]
        self.references = {}
        """
        todo: check if db file exists.
        """

    def createCommands(self):
        for i, query in self.query_files.items():
            result_file = os.path.join(self.workDir, "alignment{0}.out".format(i))

            cmd = self.aligner.get_extended_command(query, self.database, result_file)
            self.commands.append(cmd)

    def parse_result(self, file_name):
        """
        to parse alignment result.
        -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle'

        :param file_name:
        :return:
        """
        sid_list = []
        with open(file_name) as f:
            data = f.readlines()
            for row in data:
                q_id, s_id, pident, length, mismatch, gapopen, qstart, qend, \
                sstart, send, evalue, bitscore, qlen, slen, stitle = row.strip("\n").split("\t")
                # protein = self.references[s_id]
                evalue = float(evalue)
                if evalue > self.evalue_cutoff:
                    # filtering high e-value hit
                    continue
                q_cov = BlastSearch.get_coverage(qstart, qend, qlen)
                s_cov = BlastSearch.get_coverage(sstart, send, slen)
                # ProteinHit: id_, description, gene, ec_number, source_db, organism, db_name, e_value, score, identity, q_cov, s_cov, notes=[]
                flag = ""
                s_id_, product, organism, gene, source_db, ec_number = self.parser(s_id, stitle)
                hit = ProteinHit(s_id, product, gene, ec_number, source_db,
                                 organism, self.db_name, evalue, bitscore, pident, q_cov, s_cov,
                                 flag, notes=[])
                sid_list.append(s_id)
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
        return sid_list

    def set_results(self):
        sid_list = []
        for i in range(len(self.query_files)):
            result_file = os.path.join(self.workDir, "alignment{0}.out".format(i))
            sid_list += self.parse_result(result_file)
        sid_file = os.path.join(self.workDir, "sid_list.txt")
        with open(sid_file, "w") as f:
            f.write("\n".join(sid_list))

    def set_ref_info(self):
        """
        This method is for later components (e.g. PseudoGeneDetection)
        Create a fasta file containing hit proteins by using blastdbcmd.
        then, set reference info from the fasta file.
        """
        entry_file = os.path.join(self.workDir, "sid_list.txt")
        output_fasta = os.path.join(self.workDir, "reference.fasta")
        error_log_file = os.path.join(self.workDir, "blastdbcmd.err")
        cmd = self.blastdbcmd.get_command(entry_file, self.database, output_fasta, error_log_file)
        self.commands = [cmd]
        self.executeCommands(shell=True)

        if os.path.exists(output_fasta) and os.path.getsize(output_fasta) > 0:
            self.references = Protein.read_from_fasta(output_fasta, parser_type=self.dbtype)
            self.logger.debug("{} sequences were added to the reference.".format(len(self.references)))
        else:
            self.logger.warning("Fasta file for hit proteins does not exist or is empty. blastdbcmd might have failed.")

    def run(self):
        self.prepareQueries()
        self.createCommands()
        self.executeCommands(shell=True)
        self.set_results()
        self.set_ref_info()

if __name__ == '__main__':
    pass
