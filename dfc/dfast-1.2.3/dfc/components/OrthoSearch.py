#! /usr/bin/env python
# coding: UTF8

import os
from collections import namedtuple
from .baseComponent import BaseAnnotationComponent
from ..models.protein import Protein
from ..models.hit import ProteinHit

Alignment = namedtuple("Alignment",
                       ["q_id", "s_id", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore"])


class OrthoSearch(BaseAnnotationComponent):
    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        super(OrthoSearch, self).__init__(genome, options, workDir, CPU)
        self.evalue_cutoff = options.get("evalue_cutoff", 1e-5)
        self.qcov_cutoff = options.get("qcov_cutoff", 70)
        self.scov_cutoff = options.get("scov_cutoff", 70)

        self.reference_files = options.get("references", "")
        self.references = {}
        self.results = {}
        self.all_hits = {}

        self.ref_dir = os.path.join(self.workDir, "refs")
        if not os.path.exists(self.ref_dir):
            os.makedirs(self.ref_dir)
        self.result_dir = os.path.join(self.workDir, "result")
        if not os.path.exists(self.result_dir):
            os.makedirs(self.result_dir)


    def prepare_references(self):
        for i, file_name in enumerate(self.reference_files):
            D = Protein.read_reference(file_name)  # file format will be automatically inferred.
            output_fasta_file = os.path.join(self.ref_dir, "reference{0}.faa".format(i))
            Protein.write_as_fasta(D, output_fasta_file)
            self.references.update(D)

    def format_dbs(self):
        for i in range(len(self.reference_files)):
            fasta_file = os.path.join(self.ref_dir, "reference{0}.faa".format(i))
            assert os.path.exists(fasta_file)
            db_file = os.path.join(self.ref_dir, "reference{0}".format(i))
            self.commands.append(self.aligner.format_db_command(fasta_file, db_file))

        # query file
        query_fasta_file = self.query_files[0]
        assert os.path.exists(query_fasta_file)
        query_db = os.path.join(self.workDir, "query0")
        self.commands.append(self.aligner.format_db_command(query_fasta_file, query_db))
        self.executeCommands(shell=True, process_name="Format DB")

    def createCommands(self):
        query_fasta_file = self.query_files[0]
        query_db = os.path.join(self.workDir, "query0")

        for i in range(len(self.reference_files)):
            ref_fasta_file = os.path.join(self.ref_dir, "reference{0}.faa".format(i))
            ref_db = os.path.join(self.ref_dir, "reference{0}".format(i))
            result_file_f = os.path.join(self.result_dir, "query0-reference{0}.out".format(i))
            result_file_r = os.path.join(self.result_dir, "reference{0}-query0.out".format(i))

            cmd = self.aligner.get_command(query_file=query_fasta_file, db_name=ref_db, result_file=result_file_f)
            cmd.append("2> /dev/null")  # Drop warning for O containing query sequence
            self.commands.append(cmd)
            cmd = self.aligner.get_command(query_file=ref_fasta_file, db_name=query_db, result_file=result_file_r)
            cmd.append("2> /dev/null")  # Drop warning for O containing query sequence
            self.commands.append(cmd)

        # self-self alignment
        result_file_self = os.path.join(self.result_dir, "query0-query0.out")
        cmd = self.aligner.get_self_alignment_command(query_file=query_fasta_file, db_name=query_db, result_file=result_file_self)
        cmd.append("2> /dev/null")  # Drop warning for O containing query sequence
        self.commands.append(cmd)

    def read_inner_hits(self):
        """
        key (query_id_1, query_id_2)
        value bit_score
        """
        result_file_self = os.path.join(self.result_dir, "query0-query0.out")
        D = {}
        with open(result_file_self) as f:
            data = f.readlines()
            for row in data:
                q_id, s_id, pident, length, mismatch, gapopen, qstart, qend, \
                sstart, send, evalue, bitscore = row.strip("\n").split("\t")
                evalue = float(evalue)
                if evalue > self.evalue_cutoff:
                    continue
                D[(q_id, s_id)] = Alignment(q_id, s_id, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)
        return D

    def read_hits(self, qry, target):
        result_file_self = os.path.join(self.result_dir, "{0}-{1}.out".format(qry, target))
        D = {}
        with open(result_file_self) as f:
            data = f.readlines()
            for row in data:
                q_id, s_id, pident, length, mismatch, gapopen, qstart, qend, \
                sstart, send, evalue, bitscore = row.strip("\n").split("\t")
                evalue = float(evalue)
                if evalue > self.evalue_cutoff:
                    continue
                D[q_id] = Alignment(q_id, s_id, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)
        return D

    def find_orthologues(self):

        # ProteinHit: id_, description, gene, ec_number, source_db, organism, db_name, e_value, score, identity, q_cov, s_cov, notes=[]
        def _get_protein_hit(alignment, flag):
            protein = self.references[alignment.s_id]
            qlen = len(self.query_sequences[alignment.q_id])
            slen = len(protein.sequence)
            q_cov = OrthoSearch.get_coverage(alignment.qstart, alignment.qend, qlen)
            s_cov = OrthoSearch.get_coverage(alignment.sstart, alignment.send, slen)
            return ProteinHit(alignment.s_id, protein.description, protein.gene, protein.ec_number, protein.source_db,
                             protein.organism, self.__class__.__name__, alignment.evalue,
                             alignment.bitscore, alignment.pident, q_cov, s_cov,
                             flag, notes=[])

        def _filter_hit(hit):
            return hit.q_cov >= self.qcov_cutoff and hit.s_cov >= self.scov_cutoff

        inner_hits = self.read_inner_hits()
        for i in range(len(self.reference_files)):
            forward_hits = self.read_hits("query0", "reference{0}".format(i))
            reverse_hits = self.read_hits("reference{0}".format(i), "query0")

            for q_id, forward_hit in forward_hits.items():
                hit = _get_protein_hit(forward_hit, flag="")
                self.all_hits.setdefault(q_id, []).append(hit)
                reverse_hit = reverse_hits.get(forward_hit.s_id)
                if reverse_hit is None:
                    continue
                else:
                    if q_id == reverse_hit.s_id:
                        # RBH
                        hit.flag = "RBH"
                        if _filter_hit(hit):
                            self.results.setdefault(q_id, []).append(hit)
                        else:
                            self.all_hits.setdefault(q_id, []).append(hit)

                    else:
                        outer_score = float(forward_hit.bitscore)
                        inner_hit = inner_hits.get((q_id, reverse_hit.s_id))
                        if inner_hit is None:
                            continue
                        else:
                            inner_score = float(inner_hit.bitscore)
                            if inner_score >= outer_score:
                                hit.flag = "inparalogue"
                                if _filter_hit(hit):
                                    self.results.setdefault(q_id, []).append(hit)

    def set_results(self):
        hypothetical_proteins = ["hypothetical protein", "conserved protein", "uncharacterized protein", "conserved hypothetical protein"]
        for q_id, hits in self.results.items():
            # print(q_id, list(map(str, hits)))
            hits.sort(key=lambda x: -1 * float(x.score))
            best_hit = hits[0]
            feature = self.genome.features[q_id]
            # print(best_hit)
            if best_hit.description.lower() in hypothetical_proteins:
                # hypothetical protein will be added as a secondary hit.
                feature.secondary_hits.append(best_hit)
            elif best_hit.q_cov < self.qcov_cutoff or best_hit.s_cov < self.scov_cutoff:
                best_hit.flag += "Partial hit, "
                feature.secondary_hits.append(best_hit)
            elif feature.primary_hit is None:
                feature.primary_hit = best_hit

            else:
                best_hit.flag = "Secondary hit, " + best_hit.flag
                feature.secondary_hits.append(best_hit)

        # add partial or alternative hit to secondary_hits
        for q_id, hits in self.all_hits.items():
            if q_id in self.results:
                continue
            else:
                hits.sort(key=lambda x: -1 * float(x.score))
                best_hit = hits[0]
                feature = self.genome.features[q_id]
                if best_hit.q_cov < self.qcov_cutoff or best_hit.s_cov < self.scov_cutoff:
                    best_hit.flag = "Partial"
                feature.secondary_hits.append(best_hit)

    def run(self):
        self.prepareQueries(split_query=False)
        self.prepare_references()
        self.format_dbs()
        self.createCommands()
        self.executeCommands(shell=True)
        self.find_orthologues()
        self.set_results()


if __name__ == '__main__':
    pass
