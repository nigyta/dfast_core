#! /usr/bin/env python
# coding: UTF8


# ! /usr/bin/env python
# coding: UTF8

import os
from .baseComponent import BaseAnnotationComponent
from ..models.cdd_model import CDDmodel
from ..models.hit import CddHit
from ..tools.rpsblast import RPSblast
from ..tools.rpsbproc import Rpsbproc


class CDDsearch(BaseAnnotationComponent):
    instances = 0

    def __init__(self, genome, options, workDir, CPU):
        super(CDDsearch, self).__init__(genome, options, workDir, CPU)
        self.rpsblast = RPSblast(options=options)
        self.rpsbproc = Rpsbproc(options=options)
        self.database = options.get("database", "")

    def create_rpsblast_commands(self):
        for i, query in self.query_files.items():
            result_file = os.path.join(self.workDir, "alignment{0}.xml".format(i))

            cmd = self.rpsblast.get_command(query, self.database, result_file)
            self.commands.append(cmd)

    def create_rpsbproc_commands(self):
        for i in range(len(self.query_files)):
            rpsblast_result_file = os.path.join(self.workDir, "alignment{0}.xml".format(i))
            result_file = os.path.join(self.workDir, "rpsbproc{0}.out".format(i))

            cmd = self.rpsbproc.get_command(rpsblast_result_file, result_file)
            self.commands.append(cmd)

    def parse_result(self):

        def _read_rpsbproc_result(file_name):
            '''
            The file structure of rpsbproc is like below.

            #QUERY  <query-id>  <seq-type>  <seq-length>    <definition-line>
            #DOMAINS
            #<session-ordinal>  <query-id[readingframe]>    <hit-type>  <PSSM-ID>   <from>  <to>    <E-Value>   <bitscore>  <accession> <short-name>    <incomplete>    <superfamily PSSM-ID>
            #more such lines......
            #ENDDOMAINS
            #SITES
            #<session-ordinal>  <query-id[readingframe]>    <annot-type>    <title> <residue(coordinates)>  <complete-size> <mapped-size>   <source-domain>
            #more such lines......
            #ENDSITES
            #MOTIFS
            #<session-ordinal>  <query-id[readingframe]>    <title> <from>  <to>    <source-domain>
            #more such lines......
            #ENDMOTIFS
            #ENDQUERY   <query-id>
            '''

            for line in open(file_name):
                line = line.strip()
                if line == "":
                    return  # Bug fix for Python 3.7
                    # raise StopIteration  # the end of data
                elif line.startswith("END") or line.startswith("#") or line.startswith("SESSION") or line.startswith(
                        "DATA"):
                    continue
                elif line.startswith("QUERY"):
                    # parse query
                    fields = line.split("\t")
                    query_id = fields[4].split()[0]

                elif line.startswith("DOMAINS") or line.startswith("SITES") or line.startswith("MOTIFS"):
                    result_type = line[:-1].title()  # remove trailing 'S'

                else:
                    # data line
                    fields = line.split("\t")
                    assert int(fields[0])  # the first column must be a number
                    if result_type == "Domain":
                        num, _q_id, hit_type, pssm_id, from_, to_, evalue, score, accession, short_name, incomplete, superfam_pssm_id = fields
                        yield query_id, result_type, hit_type, pssm_id, from_, to_, evalue, score, accession, short_name, incomplete
                    elif result_type == "Site":
                        pass
                    elif result_type == "Motif":
                        pass
                        # This is the end of "_read_rpsbproc_result"

        cdd_definition_file = os.path.join(self.rpsbproc.rpsbproc_data, "cddid.tbl")
        cdd_definitions = CDDmodel.read_cdd_definition(cdd_definition_file)

        for i in range(len(self.query_files)):
            result_file = os.path.join(self.workDir, "rpsbproc{0}.out".format(i))
            for (query_id, result_type, hit_type, pssm_id, from_, to_, evalue, score, accession,
                 short_name, incomplete) in _read_rpsbproc_result(result_file):

                cdd = cdd_definitions.get(pssm_id)
                cdd_description = cdd.description.replace('"', "'")
                cdd_hit = CddHit(result_type, hit_type, pssm_id, from_, to_, evalue, score, accession, short_name, incomplete, cdd_description)
                feature = self.genome.features[query_id]
                feature.secondary_hits.append(cdd_hit)

    def run(self):
        self.prepareQueries()
        self.create_rpsblast_commands()
        self.executeCommands(shell=True, process_name="RPSblast")
        self.create_rpsbproc_commands()
        self.executeCommands(shell=True, process_name="rpsbproc")
        self.parse_result()
        # self.set_results()


if __name__ == '__main__':
    pass
