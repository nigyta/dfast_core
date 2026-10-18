#! /usr/bin/env python
# coding: UTF8

"""Map nucleotide references (CARD/VFDB) onto the genome with blastn and annotate the mapped regions."""

import os
from .baseComponent import BaseAnnotationComponent
from .NuclSearch import models
from ..utils.mapping_util import BLAST_OUTFMT, read_blast_table, is_full_length, group_by_locus, is_valid_reference


class NuclMapping(BaseAnnotationComponent):
    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        options["aligner"] = "blastn"
        super(NuclMapping, self).__init__(genome, options, workDir, CPU)
        self.pipeline_dir = workDir
        self.database = options.get("database", "")
        self.db_name = options.get("db_name", "")
        self.pident_cutoff = options.get("pident_cutoff", 90)
        self.evalue_cutoff = options.get("evalue_cutoff", 1e-10)
        self.overlap_tolerance = options.get("overlap_tolerance", 30)
        self.transl_table = options.get("transl_table", 11)
        self.aligner_options = options.get("aligner_options", "")
        self.result_file = os.path.join(self.workDir, "mapping.out")
        self.hits_json = {}
        self.get_db_version()
        self.load_references()

    def get_db_version(self):
        version_file = self.database + ".version"
        if os.path.exists(version_file):
            db_version = open(version_file).read().strip()
            self.logger.info(f"Database version: {self.db_name} {db_version}")
        else:
            self.logger.warning(f"Database version could not be found: {self.db_name} UNDETERMINED")

    def load_references(self):
        self.references = models[self.db_name].load_from_tsv_file(self.database + ".nucl.ref")
        self.ref_validity = {ref_id: is_valid_reference(model.nucl_seq, self.transl_table)
                             for ref_id, model in self.references.items()}

    def genome_fasta(self):
        return os.path.join(self.pipeline_dir, "input", "genome.fna")

    def create_db(self):
        db_prefix = os.path.join(self.workDir, "genome")
        self.executeCommand(self.aligner.format_db_command(self.genome_fasta(), db_prefix))
        return db_prefix

    def blast_command(self, db_prefix):
        return ["blastn", "-task", "dc-megablast", "-query", self.database + ".nucl.fasta", "-db", db_prefix,
                "-out", self.result_file, "-outfmt", "'" + BLAST_OUTFMT + "'",
                "-evalue", str(self.evalue_cutoff), "-perc_identity", str(self.pident_cutoff),
                "-max_target_seqs", "1000", "-max_hsps", "50", "-num_threads", str(self.CPU)] \
            + self.aligner_options.split()

    def select_hits(self):
        hits = []
        for hit in read_blast_table(self.result_file):
            if not is_full_length(hit, self.pident_cutoff):
                continue
            if hit.seq_id not in self.genome.seq_records:
                self.logger.warning(f"Mapping hit {hit.ref_id} refers to an unknown sequence ({hit.seq_id}); ignored.")
                continue
            hits.append(hit)
        return group_by_locus(hits)

    def apply(self, groups):
        pass  # Task 6

    def write_hits(self):
        pass  # Task 7

    def run(self):
        db_prefix = self.create_db()
        self.executeCommand(self.blast_command(db_prefix), shell=True)
        self.apply(self.select_hits())
        self.write_hits()
