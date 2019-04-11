#! /usr/bin/env python
# coding: UTF8

import os
from Bio import SeqIO
from .DBsearch import DBsearch
from .baseComponent import BaseAnnotationComponent
from ..models.hit import MBGDHit
from ..tools.blastdbcmd import Blastdbcmd

from ..utils.ref_util import fasta_parsers  # keys:  auto/ncbi/uniprot/plain
from ..models.protein import Protein


class MBGDreference():
    def __init__(self, seq_id, clst_id, clst_descr, gene_symbol, gene_descr, sequence):
        self.seq_id = seq_id
        self.clst_id = clst_id
        self.gene_descr = gene_descr
        self.gene_symbol = gene_symbol
        self.clst_descr = clst_descr
        self.sequence = sequence


def MBGD_fasta_reader(fasta_file_name, mbgd_definition_file):
    cluster_info = MBGDcluster.read_cluster_info(mbgd_definition_file)
    # example) hin:HI1272 [cluster=1] ABC transporter ATP-binding protein
    R = list(SeqIO.parse(fasta_file_name, "fasta"))
    D = {}
    for r in R:
        seq_id, clst_id, gene_descr = r.description.split(" ", 2)
        clst_id = clst_id.split("=")[-1].strip("]")
        cluster = cluster_info.get(clst_id)
        if cluster:
            clst_descr, gene_symbol = cluster.descr, cluster.gene
        else:
            clst_descr, gene_symbol = "", ""
        ref = MBGDreference(seq_id, clst_id, clst_descr, gene_symbol, gene_descr, str(r.seq))
        D[seq_id] = ref
    return D


class MBGDsearch(BaseAnnotationComponent):
    # instances = 0

    def __init__(self, genome, options, workDir, CPU):
        super(MBGDsearch, self).__init__(genome, options, workDir, CPU)
        self.evalue_cutoff = options.get("evalue_cutoff", 1e-5)
        self.qcov_cutoff = options.get("qcov_cutoff", 70)
        self.scov_cutoff = options.get("scov_cutoff", 70)
        self.pid_cutoff = options.get("pid_cutoff", 40)
        self.mbgd_version = options.get("mbgd_version", "unknown_version")
        self.mbgd_tabid = options.get("mbgd_tabid", "unknown_tabid")
        self.mbgd_definition_file = options.get("mbgd_definition_file", "mbgd_2015-01_default.tab")
        self.database = options.get("database", "")
        aligner_name = options.get("aligner")
        assert aligner_name == "ghostx"
        self.blastdbcmd = Blastdbcmd()
        # check_db_file(self.database, aligner_name)
        # self.dbtype = options.get("mbgd")
        # assert self.dbtype == "mbgd"
        # 
        # self.parser = fasta_parsers[self.dbtype]
        self.commands = []
        self.references = {}
        """
        todo: check if db file exists.
        """

    def createCommands(self):
        for i, query in self.query_files.items():
            # example of ghostz ["ghostz", "aln", "-i", queryFileName, "-d", dbFileName, "-o", outFileName, "-b", "1"]
            # query_file = self.querie_files["query{0}".format(i)]
            result_file = os.path.join(self.workDir, "alignment{0}.out".format(i))

            cmd = self.aligner.get_command(query, self.database, result_file)
            self.commands.append(cmd)


    def set_results(self):
        for i in range(len(self.query_files)):
            result_file = os.path.join(self.workDir, "alignment{0}.out".format(i))
            self.parse_result(result_file)

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
                s_id = s_id.split(" ")[0]
                ref = self.references[s_id]
                evalue = float(evalue)
                pident = float(pident)
                if evalue > self.evalue_cutoff or pident < self.pid_cutoff:
                    # filtering high e-value hit
                    continue
                qlen = len(self.query_sequences[q_id])
                slen = len(ref.sequence)
                q_cov = self.get_coverage(qstart, qend, qlen)
                s_cov = self.get_coverage(sstart, send, slen)
                if q_cov < self.qcov_cutoff or s_cov < self.scov_cutoff:
                    continue
                hit = MBGDHit(s_id, self.mbgd_version, self.mbgd_tabid, ref.clst_id, ref.clst_descr, ref.gene_symbol,
                              ref.gene_descr, evalue, bitscore, pident, q_cov, s_cov, flag="", notes=None)
                feature = self.genome.features[q_id]
                feature.secondary_hits.append(hit)

    def set_ref_info(self):
        """
        This method is to collect s_id from hit results,
        and use BlastDBcmd to create a fasta file of hit proteins.
        Then, set reference dict from the fasta file.
        """
        result_files = os.path.join(self.workDir, "alignment*.out")
        entry_file = os.path.join(self.workDir, "sid_list.txt")
        cmd = ['cut -f2 {0} | cut -f1 -d" " > {1}'.format(result_files, entry_file)]
        self.commands = [cmd]
        self.executeCommands(shell=True, process_name="collecting seqids")

        reference_fasta = os.path.join(self.workDir, "reference.fasta")
        error_log_file = os.path.join(self.workDir, "blastdbcmd.err")
        cmd = self.blastdbcmd.get_command(entry_file, self.database, reference_fasta, error_log_file)
        self.commands = [cmd]
        self.executeCommands(shell=True, process_name="BLASTDBCMD")

        if os.path.exists(reference_fasta) and os.path.getsize(reference_fasta) > 0:
            self.references = MBGD_fasta_reader(reference_fasta, self.mbgd_definition_file)
            self.logger.debug("{} sequences were added to the reference.".format(len(self.references)))
        else:
            self.logger.warning("Fasta file for hit proteins does not exist or is empty. blastdbcmd might have failed.")

    def run(self):
        self.prepareQueries()
        self.createCommands() # inherited from DBsearch
        self.executeCommands(shell=True)
        self.set_ref_info()
        self.set_results()  # inherited from DBsearch


class MBGDcluster():
    def __init__(self, id_, gene, descr):
        self.id = id_
        self.gene = gene
        self.descr = descr

    def __repr__(self):
        return "{0.id} {0.gene} {0.descr}".format(self)

    @staticmethod
    def read_cluster_info(descr_file_name):
        # descr_file_name = "mbgd_2015-01_default.tab"
        D = {}
        f = open(descr_file_name)
        for line in f:
            if line.startswith("Cluster"):
                break
        for line in f:
            cols = line.strip("\n").split("\t")
            cluster_id, gene, descr = cols[0], cols[2], cols[7]
            cluster=MBGDcluster(cluster_id, gene, descr)
            D[cluster_id] = cluster
        return D

if __name__ == '__main__':
    pass
