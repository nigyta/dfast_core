#! /usr/bin/env python
# coding: UTF8

import json
import os
from .baseComponent import BaseAnnotationComponent
from ..models.protein import Protein
from ..models.hit import ProteinHit, NuclHit
from ..utils.reffile_util import check_db_file, read_db_attributes
from ..tools.blastdbcmd import Blastdbcmd
from ..models.nucref import PLASMID_DB  # , NucRef, NucRefBase
from ..models.card import CARD
from ..models.vfdb import VFDB
# from ..models.plasmid_db import PLASMID_DB


models = {
    "PLASMID_DB": PLASMID_DB,
    "CARD": CARD,
    "VFDB": VFDB
}

class NuclSearch(BaseAnnotationComponent):
    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        super(NuclSearch, self).__init__(genome, options, workDir, CPU)
        self.evalue_cutoff = options.get("evalue_cutoff", 1e-5)
        self.qcov_cutoff = options.get("qcov_cutoff", 90)
        self.scov_cutoff = options.get("scov_cutoff", 90)
        self.pident_cutoff = options.get("pident_cutoff", 90)

        self.database = options.get("database", "")  
        self.db_name = options.get("db_name", "")
        self.get_db_version()
        self.ref_model = models[self.db_name]
        aligner_name = "blastn"  # options.get("aligner")
        assert aligner_name == "blastn"
        self.blastdbcmd = Blastdbcmd()
        # check_db_file(self.database, aligner_name)

        # self.dbtype = options.get("dbtype")
        # self.parser = fasta_parsers[self.dbtype]

        # self.references = {}
        """
        todo: check if db file exists.
        """

    def get_db_version(self):
        versoin_file = self.database + ".version"
        if os.path.exists(versoin_file):
            db_version = open(versoin_file).read().strip()
            self.logger.info(f"Database version: {self.db_name} {db_version}")
        else:
            db_version = "UNDETERMINED"
            self.logger.warning(f"Database version could not be found: {self.db_name} {db_version}")

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

        ref_tsv_file = self.database + ".nucl.ref"
        ref_models = self.ref_model.load_from_tsv_file(ref_tsv_file)
        sid_list = []
        hits = {}

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
                q_cov = NuclSearch.get_coverage(qstart, qend, qlen)
                s_cov = NuclSearch.get_coverage(sstart, send, slen)
                if not (q_cov > self.qcov_cutoff and s_cov > self.scov_cutoff):
                    continue
                # ProteinHit: id_, description, gene, ec_number, source_db, organism, db_name, e_value, score, identity, q_cov, s_cov, notes=[]
                flag = ""
                
                # to_be_implemented
                # s_id_, product, organism, gene, source_db, ec_number = self.parser(s_id, stitle)
                s_id, product, organism, gene, source_db, ec_number = s_id, s_id.split("|")[-1], "-", "-", "-", "-"
                model = ref_models[s_id]
                hit = NuclHit(s_id, model, self.db_name, evalue, bitscore, pident, q_cov, s_cov,
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
                hits[feature.id] = hit.to_dict()
        # sid_list will be used for PsuedoGeneDetection (not implemented yet)
        # hits will be used for json output and to make a summary file for AMR detection
       
        return sid_list, hits

    def set_results(self):
        sid_list = []
        hits = {}
        for i in range(len(self.query_files)):
            result_file = os.path.join(self.workDir, "alignment{0}.out".format(i))
            tmp_sid_list, tmp_hits = self.parse_result(result_file)
            sid_list += tmp_sid_list
            hits.update(tmp_hits)
        sid_file = os.path.join(self.workDir, "sid_list.txt")
        hits_json_file = os.path.join(self.workDir, "nucl_hits.json")
        json.dump(hits, open(hits_json_file, "w"), indent=4)


    # def set_ref_info(self):
    #     """
    #     This method is for later components (e.g. PseudoGeneDetection)
    #     Create a fasta file containing hit proteins by using blastdbcmd.
    #     then, set reference info from the fasta file.
    #     """
    #     entry_file = os.path.join(self.workDir, "sid_list.txt")
    #     output_fasta = os.path.join(self.workDir, "reference.fasta")
    #     error_log_file = os.path.join(self.workDir, "blastdbcmd.err")
    #     cmd = self.blastdbcmd.get_command(entry_file, self.database, output_fasta, error_log_file)
    #     self.commands = [cmd]
    #     self.executeCommands(shell=True)

    #     if os.path.exists(output_fasta) and os.path.getsize(output_fasta) > 0:
    #         self.references = Protein.read_from_fasta(output_fasta, parser_type=self.dbtype)
    #         self.logger.debug("{} sequences were added to the reference.".format(len(self.references)))
    #     else:
    #         self.logger.warning("Fasta file for hit proteins does not exist or is empty. blastdbcmd might have failed.")

    def run(self):
        # if not self.enabled:
        #     self.logger.warning("[Warning] {} is disabled. Skip execution.".format(self.__class__.__name__))
        self.prepareQueries(type_="nucleotide")
        self.createCommands()
        self.executeCommands(shell=True)
        self.set_results()


if __name__ == '__main__':
    pass
