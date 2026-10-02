#! /usr/bin/env python
# coding: UTF8


import os
from .baseComponent import BaseAnnotationComponent
from ..models.hit import HmmHit
from ..tools.hmmer import Hmmer_hmmscan
from ..utils.reffile_util import check_hmm_db_file

# Family types of the NCBI HMM collection, from the most to the least specific (as used by PGAP for naming).
FAMILY_TYPE_PRIORITY = ["exception", "equivalog", "subfamily", "equivalog_domain", "subfamily_domain", "domain"]


def read_hmm_attributes(file_name):
    """Read hmm_PGAP.tsv of the NCBI HMM collection: {ncbi_accession: {family_type, for_naming, product_name, ...}}"""
    attributes = {}
    with open(file_name) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            attributes[fields[0]] = {"family_type": fields[6], "for_naming": fields[8], "product_name": fields[10],
                                     "gene_symbol": fields[11], "ec_number": fields[13]}
    return attributes


def select_naming_hit(hmms):
    """
    The hit used to name the protein: an HMM marked for naming whose product is not a hypothetical protein,
    preferring the most specific family type, then the highest score. None if there is no such hit.
    """
    candidates = [hmm for hmm in hmms if hmm.attributes and hmm.attributes["for_naming"] == "Y"
                  and hmm.attributes["product_name"] and "hypothetical protein" not in hmm.attributes["product_name"].lower()]
    if not candidates:
        return None
    rank = {family_type: i for i, family_type in enumerate(FAMILY_TYPE_PRIORITY)}
    return min(candidates, key=lambda hmm: (rank.get(hmm.attributes["family_type"], len(rank)), -hmm.score))


class HMMscan(BaseAnnotationComponent):
    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        super(HMMscan, self).__init__(genome, options, workDir, CPU)
        self.hmmer = Hmmer_hmmscan(options=options)
        self.database = options.get("database", "")
        self.db_name = options.get("db_name", "")
        check_hmm_db_file(self.database)
        # Attribute table of the NCBI HMM collection (hmm_PGAP.tsv). When given, hits are used to name proteins.
        attribute_file = options.get("attributes", "")
        if attribute_file and not os.path.exists(attribute_file):
            self.logger.error("HMM attribute file ({}) does not exist. Aborting...".format(attribute_file))
            exit(1)
        self.attributes = read_hmm_attributes(attribute_file) if attribute_file else {}

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
            hmm = HmmHit(accession, name, description, evalue, score, bias, self.db_name, self.attributes.get(accession))
            yield query_id, hmm

    def set_results(self):

        self.logger.info("Collecting HMMscan results.")
        hitDict = {}
        for i in range(len(self.query_files)):
            result_file = os.path.join(self.workDir, "result{0}.out".format(i))
            for query_id, hmm in self.parseResult(result_file):
                hitDict.setdefault(query_id, []).append(hmm)

        named = 0
        for query_id, hmms in hitDict.items():
            # sort by score in descending order. The best hit comes first.
            hmms.sort(key=lambda x: x.score, reverse=True)
            best_hmmhit = hmms[0]
            feature = self.genome.features[query_id]
            naming_hit = select_naming_hit(hmms) if feature.primary_hit is None else None
            if naming_hit:
                feature.primary_hit = naming_hit  # gives product, gene, EC_number and inference
                named += 1
                self.logger.debug("Named {} by {}: {}".format(query_id, naming_hit.accession, naming_hit.attributes["product_name"]))
            if best_hmmhit is not naming_hit:
                feature.secondary_hits.append(best_hmmhit)
        if self.attributes:
            self.logger.info("{} CDSs were named by {} HMMs.".format(named, self.db_name))

    def run(self):
        self.prepareQueries()
        self.createCommands()
        self.executeCommands(shell=False)
        self.set_results()

if __name__ == '__main__':
    pass
