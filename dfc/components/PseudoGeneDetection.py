#! /usr/bin/env python
# coding: UTF8

import os
import re
from collections import namedtuple
from Bio.SeqFeature import FeatureLocation
from .baseComponent import BaseAnnotationComponent
from ..utils.ddbj_submission import get_location_string
from ..models.protein import Protein
from ..models.hit import ProteinHit, PseudoGene
from ..tools.last import Lastal, Lastdb
from ..tools.blastp import Blastp

Extended_CDS = namedtuple("Extended_CDS", ["id", "seq", "cds_start", "cds_end", "abs_start", "abs_end", "strand"])
Alignment = namedtuple("Alignment", ["id", "start", "aligned_length", "strand", "total_length", "alignment"])


class PseudoGeneDetection(BaseAnnotationComponent):
    instances = 0
    PAT_SCORE = re.compile(r"a score=(\d+) EG2=(.+) E=(.+)")

    def __init__(self, genome, options, workDir, CPU):
        super(PseudoGeneDetection, self).__init__(genome, options, workDir, CPU)
        self.extension = options.get("extension", 300)
        self.scov_cutoff = options.get("scov_cutoff", 85)
        self.min_extension = 50
        self.candidates = []  # candidates for pseudogenes. should be tuple (query_cds.id, protein.id)
        # self.database = options.get("database", "")
        self.lastdb = Lastdb(options)
        self.lastal = Lastal(options)
        self.components = []
        self.references = {}

        self.trans_except = []
        self.blastp = Blastp()


    def prepare_queries(self):

        def _get_nucseq(feature, seq_record, extension):
            start = feature.location.start - extension
            cds_start = extension
            cds_end = cds_start + len(feature.location)
            if start < 0:
                offset = 0 - start
                start = start + offset
                cds_start = cds_start - offset
                cds_end = cds_end - offset
            end = feature.location.end + extension
            if end > len(seq_record):
                offset = end - len(seq_record)
                end = end - offset
            if feature.location.strand == -1:
                # swap cds_start and cd_end
                length = end - start
                cds_start, cds_end = length - cds_end, length - cds_start
            extended_location = FeatureLocation(start=start, end=end, strand=feature.strand)
            nucseq = extended_location.extract(seq_record)

            # for debug
            # print(feature.id, start, end, feature.location.strand, cds_start, cds_end, str(nucseq.seq), str(nucseq.seq)[cds_start:cds_end])
            # print(str(nucseq.seq)[cds_start:cds_end])
            # print(str(feature.extract(seq_record).seq))

            assert str(nucseq.seq)[cds_start:cds_end] == str(feature.extract(seq_record).seq)
            return Extended_CDS(feature.id, str(nucseq.seq), int(cds_start), int(cds_end),
                                int(start), int(end), int(feature.location.strand))

        def _generate_partial_hit(features):
            """
            :param features: list(iterator)  of features. (self.features.values())
            :return: feature_object,
            """
            for feature in features:
                if feature.primary_hit and isinstance(feature.primary_hit, ProteinHit):
                    if feature.primary_hit.s_cov < self.scov_cutoff:
                        yield feature, feature.primary_hit.id
                for hit in feature.secondary_hits:
                    if isinstance(hit, ProteinHit) and hit.s_cov < self.scov_cutoff:
                        yield feature, hit.id

        # cds_features = [feature for feature in self.genome.features.values() if feature.type == "CDS"]
        candidate_features = []
        for feature, hit_id in _generate_partial_hit(self.genome.features.values()):
            self.candidates.append((feature.id, hit_id))
            ext_cds = _get_nucseq(feature, self.genome.seq_records[feature.seq_id], extension=self.extension)
            str_fasta = ">{0}\n{1}\n".format(ext_cds.id, ext_cds.seq)
            candidate_features.append(str_fasta)
            self.query_sequences[ext_cds.id] = ext_cds

        # Currently, does not create splitted queries. split_num is set to 1.
        split_num = 1
        for i, list_fasta in enumerate(self.splitList(candidate_features, split_num)):
            fasta_buffer = "".join(list_fasta)
            # for str_fasta in list_fasta:
            #     fasta_buffer += str_fasta
            fasta_file = os.path.join(self.workDir, "query{0}.fna".format(i))
            self.query_files[i] = fasta_file
            with open(fasta_file, "w") as f:
                f.write(fasta_buffer)

        # write query info (for debug)
        buffer = ""
        query_info_file = os.path.join(self.workDir, "query.info")
        for query in self.query_sequences.values():
            buffer += "\t".join(map(str, query)) + "\n"
        with open(query_info_file, "w") as f:
            f.write(buffer)
        # print(self.candidates)

    def set_components(self, list_components):
        self.components = list_components

    def set_references(self):
        for component in self.components:
            self.references.update(component.references)

        reference_ids = {ref_id for feature_id, ref_id in self.candidates}

        fasta_buffer = ""
        for ref_id in reference_ids:
            ref = self.references.get(ref_id)
            if ref is None:
                continue
            assert isinstance(ref, Protein)
            fasta_buffer += ref.to_fasta()
        fasta_file = os.path.join(self.workDir, "reference.faa")
        with open(fasta_file, "w") as f:
            f.write(fasta_buffer)

    def format_db(self):
        fasta_file = os.path.join(self.workDir, "reference.faa")
        assert os.path.exists(fasta_file)
        db_name = os.path.join(self.workDir, "reference")
        self.commands.append(self.lastdb.get_command(fasta_file, db_name))
        self.executeCommands(shell=True, process_name="Last format DB")

    def align(self):
        db_name = os.path.join(self.workDir, "reference")
        for i in range(len(self.query_files)):
            fasta_file = os.path.join(self.workDir, "query{0}.fna".format(i))
            result_file = os.path.join(self.workDir, "lastal{0}.out".format(i))
            assert os.path.exists(fasta_file)
            self.commands.append(self.lastal.get_command(fasta_file, db_name, result_file))
        self.executeCommands(shell=True, process_name="Last alignment")

    def read_last_result(self, last_result_file):

        def _parse_alignment(line):
            _, id_, start, aligned_length, strand, total_length, alignment = line.strip("\n").split()
            alignment = Alignment(id_, int(start), int(aligned_length), strand, int(total_length), alignment)
            return alignment

        def _get_score(line):
            m = PseudoGeneDetection.PAT_SCORE.search(line)
            assert m is not None
            score = int(m.group(1))
            eg2 = float(m.group(2))
            evalue = float(m.group(3))
            return score, eg2, evalue

        with open(last_result_file) as f:
            for line in f:
                if line.startswith("a "):
                    score, eg2, evalue = _get_score(line)
                    line = next(f)
                    ref_alignment = _parse_alignment(line)
                    line = next(f)
                    query_alignment = _parse_alignment(line)
                    yield score, eg2, evalue, ref_alignment, query_alignment

    def scan_alignment(self, alignment, ref_alignment):

        def _to_abs(pos, query):
            """convert relative position to absolute position"""
            if query.strand == -1:
                return query.abs_end - pos
            else:
                return pos + query.abs_start

        query = self.query_sequences[alignment.id]  # query is extended_cds object

        # debug
        # print("[debug] cds_start:", query.cds_start, "cds_end:", query.cds_end, "strand:", query.strand, "abs:",
        #       query.abs_start, "-", query.abs_end)
        # print("[debug] alignment_start:", alignment.start, "alignment_length:", alignment.aligned_length,
        #       "total_length", alignment.total_length)

        D = {"insertion": set(), "deletion": set(), "stop_codon": set(), "transl_except": None}
        if (alignment.start > query.cds_start - self.min_extension) and (
                        alignment.start + alignment.aligned_length < query.cds_end + self.min_extension):
            return D  # check if aligned region spans the extended regions.

        pos = alignment.start
        for char, ref_char in zip(alignment.alignment, ref_alignment.alignment):
            if char == "-":
                pass
            elif char == "/":
                # deletion
                D["deletion"].add(_to_abs(pos, query))
                # print("/ deletion", "at", pos)
            elif char == "\\":
                # insertion
                D["insertion"].add(_to_abs(pos, query))
                # print("\\ insertion", "at", pos)
                pos += 1
            elif char == "*":
                # print("stop codon", pos, "(strand:{})".format(query.strand))
                # print("codon", query.seq[pos:pos + 3])
                abs_start = _to_abs(pos, query)
                abs_end = _to_abs(pos + 3, query)
                if query.strand == -1:
                    abs_start, abs_end = abs_end, abs_start
                # print("abs_position", abs_start, abs_end)
                strand = "-" if query.strand == -1 else "+"
                if ref_char == "U" or ref_char == "O":
                    aa = "selenocysteine" if ref_char == "U" else "pyrrolysine"
                    stop_codon_location = FeatureLocation(start=abs_start, end=abs_end, strand=query.strand)
                    D["transl_except"] = (abs_start, abs_end, strand, aa, stop_codon_location)
                else:
                    D["stop_codon"].add((abs_start, abs_end, strand))
                pos += 3
            else:
                pos += 3
        return D

    def find_pseudo(self):
        debug_out = ""
        all_results = {}
        for i in range(len(self.query_files)):
            result_file = os.path.join(self.workDir, "lastal{0}.out".format(i))
            assert os.path.exists(result_file)
            for score, eg2, evalue, ref_alignment, query_alignment in self.read_last_result(result_file):
                if not (query_alignment.id, ref_alignment.id) in self.candidates:
                    continue
                tmp_dict = all_results.setdefault(
                    query_alignment.id, {"insertion": set(), "deletion": set(), "stop_codon": set(), "transl_except": None})
                D = self.scan_alignment(query_alignment, ref_alignment)

                query = self.query_sequences[query_alignment.id]

                location = self.genome.features[query_alignment.id].location
                strand = "-" if location.strand == -1 else "+"
                debug_out += "QRY:{0} [{1}:{2}]({3}) REF:{4}\n".format(query_alignment.id,
                                                    location.start + 1, location.end, strand, ref_alignment.id)
                debug_out += "[Relative coordinate] cds: {start}-{end}, ".format(start=(query.cds_start + 1), end=query.cds_end) + \
                             "aligned: ({start}-{end})\n".format(start=(query_alignment.start + 1), end=(query_alignment.start + query_alignment.aligned_length))
                debug_out += "REF: " + ref_alignment.alignment + "\n"
                debug_out += "QRY: " + query_alignment.alignment + "\n"

                tmp_dict["stop_codon"] = tmp_dict["stop_codon"].union(D["stop_codon"])
                tmp_dict["insertion"] = tmp_dict["insertion"].union(D["insertion"])
                tmp_dict["deletion"] = tmp_dict["deletion"].union(D["deletion"])

                note = ""
                if D["transl_except"]:
                    abs_start, abs_end, strand, aa, stop_codon_location = D["transl_except"]
                    self.logger.debug("QRY:{0} [{1}..{2}]({3}) is a candidate of transl_except for {4}.".format(
                        query_alignment.id, abs_start, abs_end, strand, aa))
                    self.trans_except.append((query_alignment.id, stop_codon_location))  # region = N-term or C-term
                if D["stop_codon"]:
                    note = "Stop codon at " + ", ".join(
                        ["[{0}:{1}]({2})".format(x[0]+1, x[1], x[2]) for x in list(D["stop_codon"])]) + "."
                if D["insertion"]:
                    note = "Frameshift due to the insertion at around " + ", ".join(map(str, D["insertion"])) + "."
                if D["deletion"]:
                    note = "Frameshift due to the deletion at around " + ", ".join(map(str, D["deletion"])) + "."
                if note:
                    self.logger.debug("QRY:{0} [{1}:{2}]({3}) {4}".format(query_alignment.id,
                                                    location.start + 1, location.end, strand, note))
                    debug_out += note + "\n"
                debug_out += "\n"

        # set all results:
        fs_count = 0
        pseudo_count = 0
        for feature_id, result in all_results.items():
            feature = self.genome.features[feature_id]
            pseudo_gene = PseudoGene(result["stop_codon"], result["insertion"], result["deletion"])
            feature.secondary_hits.append(pseudo_gene)
            if len(result["stop_codon"]) > 0:
                pseudo_count += 1
            if len(result["insertion"]) > 0 or len(result["deletion"]) > 0:
                fs_count += 1
 
        self.logger.info("{} CDS features were marked as possible pseudo due to internal stop codons.".format(pseudo_count))
        self.logger.info("{} CDS features were marked as possible pseudo due to frameshift.".format(fs_count))

        final_result_file = os.path.join(self.workDir, "final_result.out")
        with open(final_result_file, "w") as f:
            f.write(debug_out)

    def find_transl_except(self):

        def _get_candidates(feature_id, stop_codon_location):
            feature = self.genome.features[feature_id]
            feature_list = list(self.genome.features.keys())
            my_idx = feature_list.index(feature_id)
            strand = feature.strand
            # self.logger.debug("CDS:{0} Location:{1}, Stop codon:{2}".format(feature_id, feature.location, stop_codon_location))
            if strand == 1:
                if stop_codon_location.end < feature.location.start:
                    other_idx = my_idx - 1
                elif feature.location.end - 3 <= stop_codon_location.start:
                    other_idx = my_idx + 1
                else:
                    self.logger.warning("[Warning in transl_except annotation] Stop codon was found outside of the CDS. Skip processing."
                      + "CDS:{0} Location:{1}, Stop codon:{2}".format(feature_id, feature.location, stop_codon_location))
                    return None
            if strand == -1:
                if stop_codon_location.end <= feature.location.start + 3:
                    other_idx = my_idx - 1
                elif feature.location.end < stop_codon_location.start:
                    other_idx = my_idx + 1
                else:
                    self.logger.warning("Stop codon was found outside of the CDS Skip processing."
                      + "CDS:{0} Location:{1}, Stop codon:{2}".format(feature_id, feature.location, stop_codon_location))
                    return None
            return feature_list[min(my_idx, other_idx)], feature_list[max(my_idx, other_idx)]

        def _concatenate_features(left_id, right_id):
            '''
            Let the N-term part 'parent' and C-term part 'child'.
            Two features are concatenated based on the parent feature.
            Child features will be removed in the downstream process.

            Returns parent's feature_id and child's feature_id
            
            If consistency check fails, returns None.
            '''

            left_feature = self.genome.features[left_id]
            right_feature = self.genome.features[right_id]
            if left_feature.type != "CDS" or right_feature.type != "CDS" or left_feature.seq_id != right_feature.seq_id:
                return None
            if left_feature.strand == right_feature.strand == 1:
                parent, child = left_feature, right_feature
                stop_codon_location = FeatureLocation(start=parent.location.end - 3, end=parent.location.end, strand=1)
            elif left_feature.strand == right_feature.strand == -1:
                parent, child = right_feature, left_feature
                stop_codon_location = FeatureLocation(start=parent.location.start, end=parent.location.start + 3, strand=-1)
            else:
                return None

            concatenated_location = FeatureLocation(start=left_feature.location.start, end=right_feature.location.end,
                                                    strand=left_feature.strand)

            seq_id = parent.seq_id
            whole_seq = self.genome.seq_records[seq_id]
            # annotations = parent.annotations.copy()
            # qualifiers = parent.qualifiers.copy()
            transl_table = parent.qualifiers.get("trasl_table", [11])[0]  # if not available, use translation table 11.

            extracted_seq = concatenated_location.extract(whole_seq.seq)
            translated_seq = str(extracted_seq.translate(table=transl_table).rstrip("*"))
            if translated_seq.count("*") == 1:
                stop_codon_pos = translated_seq.index("*") + 1
            else:
                return None  # Only one stop codon must be included in aa.

            stop_codon = str(stop_codon_location.extract(whole_seq.seq))
            if stop_codon.upper() == "TGA":  # opal > Selenocysteine, Sec, U
                # /transl_except=(pos:complement(5272379..5272381),aa:Sec)
                transl_except = "(pos:{},aa:Sec)".format(get_location_string(stop_codon_location))
                translated_seq = translated_seq.replace("*", "U")
                note_value = "codon on position {} is selenocysteine opal codon.".format(stop_codon_pos)
            elif stop_codon.upper() == "TAG":  # amber > pyrrolysine, Pyl, O
                # /transl_except=(pos:213..215,aa:Pyl ) 
                transl_except = "(pos:{},aa:Pyl)".format(get_location_string(stop_codon_location))
                translated_seq = translated_seq.replace("*", "O")
                note_value = "codon on position {} is pyrrolysine amber codon.".format(stop_codon_pos)
            else:
                return None  # stop codon must be either of TGA/TAG


            # todo: Change this to hit obejct
            parent.location = concatenated_location
            parent.qualifiers["translation"] = [translated_seq]
            parent.qualifiers["transl_except"] = [transl_except]
            parent.qualifiers.setdefault("note", []).append(note_value)
            parent.primary_hit, parent.secondary_hits = None, []

            # print("left", left_feature.location)
            # print("right", right_feature.location)
            # print(concatenated_location)
            return parent.id, child.id

        # ----- main part starts from here -----

        # candidate is a feature_id pair of two neighboring CDS features
        candidates = [_get_candidates(feature_id, stop_codon_location) for feature_id, stop_codon_location in self.trans_except]
        candidates = {x for x in candidates if x is not None}  # remove None or redundancy

        proteins_with_transl_except, removed = [], []
        for left_id, right_id in candidates:
            result = _concatenate_features(left_id, right_id)
            if result:
                proteins_with_transl_except.append(result[0])  # parent_id
                removed.append(result[1])  # child_id
        if proteins_with_transl_except:
            self.logger.info("Internal stop codons in the following {0} CDSs are translated to selenosysteine/pyrrolysine, which will be annotated with 'transl_except'. {1}".format(len(proteins_with_transl_except), proteins_with_transl_except))

        if removed:  # Removing child features
            self.logger.info("Removed {0} CDSs. {1}".format(len(removed), removed))
            for seq_record in self.genome.seq_records.values():
                seq_record.features = [feature for feature in seq_record.features if feature.id not in removed]
            self.genome.set_feature_dictionary()

        return proteins_with_transl_except

    def execute_blast(self, targets):

        def _prepare_query(targets, query_fasta_file):
            cds_features = [self.genome.features[target] for target in targets]
            fasta_buffer = ""
            for feature in cds_features:
                sequence = feature.qualifiers.get("translation", [""])[0]
                if sequence:
                    fasta_buffer += ">{0}\n{1}\n".format(feature.id, sequence)
                else:
                    self.logger.error("Translated sequence is empty. {0}".format(feature.id))
                    raise AssertionError
            with open(query_fasta_file, "w") as f:
                f.write(fasta_buffer)

        def _parse_result(blast_result_file):

            for row in open(blast_result_file):
                q_id, s_id, pident, length, mismatch, gapopen, qstart, qend, \
                    sstart, send, evalue, bitscore = row.strip("\n").split("\t")
                protein = self.references[s_id]
                evalue = float(evalue)

                qlen = len(self.genome.features[q_id].qualifiers["translation"][0])
                slen = len(protein.sequence)
                q_cov = PseudoGeneDetection.get_coverage(qstart, qend, qlen)
                s_cov = PseudoGeneDetection.get_coverage(sstart, send, slen)
                flag = ""
                hit = ProteinHit(s_id, protein.description, protein.gene, protein.ec_number, protein.source_db,
                                 protein.organism, self.__class__.__name__, evalue, bitscore, pident, q_cov, s_cov,
                                 flag, notes=[])
                feature = self.genome.features[q_id]
                if q_cov < self.scov_cutoff or s_cov < self.scov_cutoff:
                    hit.flag = "Partial hit"
                    feature.secondary_hits.append(hit)
                elif feature.primary_hit:
                    feature.secondary_hits.append(hit)
                else:
                    feature.primary_hit = hit

        query_fasta_file = os.path.join(self.workDir, "proteins_with_transl_except.faa")
        _prepare_query(targets, query_fasta_file)
        reference_fasta_file = os.path.join(self.workDir, "reference.faa")
        db_name = os.path.join(self.workDir, "reference")
        blast_result_file = os.path.join(self.workDir, "blast_result.out")

        self.executeCommand(self.blastp.format_db_command(reference_fasta_file, db_name))  # format db
        blast_cmd = self.blastp.get_command(query_fasta_file, db_name, blast_result_file)
        blast_cmd.append("2> /dev/null")  # Drop warning for O containing query sequence
        self.executeCommand(blast_cmd, shell=True)  # blastp
        _parse_result(blast_result_file)

    def run(self):
        self.prepare_queries()
        self.set_references()
        self.format_db()
        self.align()
        self.find_pseudo()
        proteins_with_transl_except = self.find_transl_except()
        if len(proteins_with_transl_except) > 0:
            self.execute_blast(targets=proteins_with_transl_except)
    
if __name__ == '__main__':
    pass
