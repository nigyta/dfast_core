import logging
from collections import OrderedDict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
import dfc.genome as genome_mod
from dfc.components.NuclMapping import NuclMapping, new_cds, new_pseudogene, relocate_cds
from dfc.models.bio_feature import ExtendedFeature
from dfc.models.card import CARD
import json
from dfc.utils.mapping_util import MappedRegion
from dfc.utils.summarize_amr import summarize_amr
from dfc.utils.summarize_pseudo import summarize_pseudo
from dfc.models.hit import NuclHit, PseudoGene
from tests.test_mapping_util import _hit, GENE
from dfc.functionalAnnotation import COMPONENTS

class FakeGenome:
    def __init__(self, records):
        self.seq_records = OrderedDict((r.id, r) for r in records)
        self.features = OrderedDict()
        self.set_feature_dictionary()
    sort_features = genome_mod.Genome.sort_features
    set_feature_dictionary = genome_mod.Genome.set_feature_dictionary

def _component(tmp_path, genome, **options):
    comp = NuclMapping.__new__(NuclMapping)          # skip __init__ (blastn version check, reference loading)
    comp.genome, comp.workDir, comp.pipeline_dir = genome, str(tmp_path), str(tmp_path)
    comp.logger = logging.getLogger("test")
    comp.db_name, comp.pident_cutoff, comp.overlap_tolerance, comp.transl_table = "CARD", 90, 30, 11
    comp.pseudogene_pident_cutoff, comp.overlap_fraction, comp.truncated_fraction = 95, 0.2, 0.7
    comp.inference = INFERENCE
    comp.result_file = str(tmp_path / "mapping.out")
    comp.nucl_refs, comp.ref_validity, comp.hits_json = {}, {}, {}
    for k, v in options.items(): setattr(comp, k, v)
    return comp

INFERENCE = "COORDINATES:alignment:blastn:2.17.0"
MGA_INFERENCE = "COORDINATES:ab initio prediction:MetaGeneAnnotator"

def test_registered():
    assert COMPONENTS["NuclMapping"] is NuclMapping

def test_select_hits_keeps_full_length_hits_on_known_contigs(tmp_path, caplog):
    genome = FakeGenome([SeqRecord(Seq("A" * 5000), id="sequence1")])
    comp = _component(tmp_path, genome)
    rows = ["ARO:1\tsequence1\t99.0\t300\t3\t0\t1\t300\t1001\t1300\t0.0\t500\t300\t5000\t300",
            "ARO:2\tsequence1\t99.0\t200\t2\t0\t50\t250\t2001\t2200\t0.0\t300\t300\t5000\t200",   # partial
            "ARO:3\tunknown\t99.0\t300\t3\t0\t1\t300\t1001\t1300\t0.0\t500\t300\t5000\t300",       # unknown contig
            "ARO:4\tunknown\t99.0\t300\t3\t0\t1\t300\t2001\t2300\t0.0\t500\t300\t5000\t300"]
    (tmp_path / "mapping.out").write_text("\n".join(rows) + "\n")
    groups = comp.select_hits()
    assert [[h.ref_id for h in g] for g in groups] == [["ARO:1"]]
    assert len([r for r in caplog.records if r.levelno == logging.WARNING]) == 1   # warned once per sequence

def test_blast_command(tmp_path):
    comp = _component(tmp_path, FakeGenome([]), database="/db/card/CARD", evalue_cutoff=1e-10, aligner_options="", CPU=4,
                      blast_task="megablast")
    cmd = comp.blast_command("/w/genome")
    assert cmd[:3] == ["blastn", "-task", "megablast"] and "/db/card/CARD.nucl.fasta" in cmd
    assert cmd[cmd.index("-perc_identity") + 1] == "90" and cmd[cmd.index("-num_threads") + 1] == "4"


def _setup(tmp_path, gene, cds=(), validity=(True, True), btop=None, ref_len=None):
    """One contig: 100 C + gene + 100 C, with predicted CDSs given as (id, start0, end, strand)."""
    rec = SeqRecord(Seq("C" * 100 + gene + "C" * 100), id="sequence1")
    rec.features = [ExtendedFeature(location=FeatureLocation(s, e, strand=st), type="CDS", id=i, seq_id="sequence1",
                                    qualifiers={"translation": ["M"], "transl_table": ["11"], "codon_start": [1],
                                                "inference": [MGA_INFERENCE]},
                                    annotations={"partial_flag": "00"}) for i, s, e, st in cds]
    comp = _component(tmp_path, FakeGenome([rec]))
    comp.nucl_refs = {"ARO:1": CARD("ARO:1", "OXA-48", "OXA beta-lactamase", "P1", "", "N1", gene, "E. coli", "")}
    comp.ref_validity = {"ARO:1": validity}
    n = ref_len or len(gene)
    send = 100 + len(gene) if btop else 100 + n   # genome-side end of the aligned gene
    hit = _hit(qlen=n, qend=n, sstart=101, send=send, slen=len(rec), btop=btop or str(n))
    return comp, rec, hit

def _cds_features(rec):
    return {f.id: f for f in rec.features if f.type == "CDS"}

def test_rescue_adds_cds_with_primary_hit(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    comp.apply([[hit]])
    f = _cds_features(rec)["MAP_1"]
    assert (int(f.location.start), int(f.location.end), f.location.strand) == (100, 100 + len(GENE), 1)
    assert f.qualifiers["translation"] == ["M" + "A" * 20] and f.primary_hit.id == "ARO:1"
    assert comp.hits_json["MAP_1"]["model"]["note"].startswith("rescued; ")
    assert comp.genome.features["MAP_1"] is f and f.qualifiers["inference"] == [INFERENCE]

def test_annotate_keeps_identical_cds(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 100, 100 + len(GENE), 1)])
    comp.apply([[hit]])
    assert list(_cds_features(rec)) == ["MGA_1"]
    f = rec.features[0]
    assert f.primary_hit.id == "ARO:1" and f.annotations["nucl_mapping"] == "annotated"
    assert not comp.hits_json["MGA_1"]["model"]["note"].startswith(("rescued", "start corrected", "pseudogene"))

def test_start_correction_relocates_same_gene(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 130, 100 + len(GENE), 1)])
    comp.apply([[hit]])
    f = _cds_features(rec)["MGA_1"]
    assert int(f.location.start) == 100 and f.qualifiers["translation"] == ["M" + "A" * 20]
    assert "start codon moved 30 bp upstream by mapping of CARD:OXA-48" in f.qualifiers["note"]
    assert f.annotations["nucl_mapping"] == "start corrected" and f.primary_hit.id == "ARO:1"
    assert comp.hits_json["MGA_1"]["model"]["note"].startswith("start corrected; ")
    assert f.qualifiers["inference"] == [INFERENCE]          # replaces the ab initio COORDINATES inference

def test_frameshift_makes_pseudogene_and_removes_fragments(tmp_path):
    gene = GENE[:10] + GENE[11:]   # 1-bp deletion at position 111 of the genome
    comp, rec, hit = _setup(tmp_path, gene, ref_len=len(GENE), btop="10A-" + str(len(GENE) - 11),
                            cds=[("MGA_1", 100, 140, 1), ("MGA_2", 140, 100 + len(gene), 1)])
    comp.apply([[hit]])
    assert list(_cds_features(rec)) == ["MAP_1"]
    f = rec.features[0]
    assert f.qualifiers["pseudogene"] == ["unknown"] and "translation" not in f.qualifiers
    assert not {"transl_table", "codon_start"} & set(f.qualifiers)   # DDBJ ANN3150: not allowed with /pseudogene
    assert "frameshift at around 111" in f.qualifiers["note"]
    assert any(isinstance(h, PseudoGene) for h in f.secondary_hits) and f.primary_hit.id == "ARO:1"
    assert comp.hits_json["MAP_1"]["model"]["note"].startswith("pseudogene: frameshift; ")
    assert f.qualifiers["inference"] == [INFERENCE]
    f.assign_hit()
    assert not [n for n in f.qualifiers["note"] if "frameshifted" in n or "possible pseudo" in n]

def test_distant_reference_does_not_call_pseudogene(tmp_path):
    gene = GENE[:10] + GENE[11:]
    cds = [("MGA_1", 100, 140, 1), ("MGA_2", 140, 100 + len(gene), 1)]
    comp, rec, hit = _setup(tmp_path, gene, ref_len=len(GENE), btop="10A-" + str(len(GENE) - 11), cds=cds)
    hit.pident = 93.6                                       # a single cross-strain reference (aslA)
    comp.apply([[hit]])
    f = _cds_features(rec)
    assert list(f) == ["MGA_1", "MGA_2"] and f["MGA_2"].primary_hit.id == "ARO:1"

def test_unresolved_region_is_skipped(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE[:-3] + "GCT", validity=(True, False))   # no stop within 100 bp
    comp.apply([[hit]])
    assert _cds_features(rec) == {}

def _add_reference(comp, ref_id, nucl_seq, validity=(True, True)):
    comp.nucl_refs[ref_id] = CARD(ref_id, "OXA-1", "OXA beta-lactamase", "P2", "", "N2", nucl_seq, "E. coli", "")
    comp.ref_validity[ref_id] = validity

def test_identical_intact_reference_beats_longer_frameshifted_one(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 130, 100 + len(GENE), 1)])
    hit.pident = 100.0
    longer = GENE[:10] + "A" + GENE[10:]                    # 98.8%, one gap against the genome (hlyE)
    _add_reference(comp, "ARO:2", longer)
    other = _hit(ref_id="ARO:2", pident=98.8, bitscore=600.0, qlen=len(longer), qend=len(longer), sstart=101,
                 send=100 + len(GENE), slen=len(rec), btop="10A-" + str(len(GENE) - 10))
    comp.apply([[other, hit]])
    f = _cds_features(rec)
    assert list(f) == ["MGA_1"] and int(f["MGA_1"].location.start) == 100
    assert f["MGA_1"].primary_hit.id == "ARO:1" and [h.id for h in f["MGA_1"].secondary_hits] == ["ARO:2"]

def test_intact_reference_beats_longer_one_with_internal_stop(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE + "GCTGCTTAA", ref_len=len(GENE))
    hit.pident = 100.0
    longer = GENE[:-3] + "TACGCTGCTTAA"                     # 3 codons longer, 99.9% (cusC)
    _add_reference(comp, "ARO:2", longer)
    other = _hit(ref_id="ARO:2", pident=99.9, bitscore=600.0, qlen=len(longer), qend=len(longer), sstart=101,
                 send=100 + len(longer), slen=len(rec), btop=str(len(longer)))
    comp.apply([[other, hit]])
    f = _cds_features(rec)
    assert list(f) == ["MAP_1"] and "pseudogene" not in f["MAP_1"].qualifiers
    assert (int(f["MAP_1"].location.start), int(f["MAP_1"].location.end)) == (100, 100 + len(GENE))
    assert f["MAP_1"].primary_hit.id == "ARO:1"

def test_other_strand_cds_is_kept_neighbour_kept(tmp_path):
    right = 100 + len(GENE)
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 110, right - 10, -1), ("MGA_2", right - 20, right + 100, 1)])
    comp.apply([[hit]])
    assert sorted(_cds_features(rec)) == ["MAP_1", "MGA_1", "MGA_2"] and comp.counts["removed"] == 0


def test_same_strand_cds_removed_only_when_overlap_is_20_percent(tmp_path):
    right = 100 + len(GENE)
    # MGA_1: 35 bp of 300 bp (12%) overlap -> kept; MGA_2: inside the gene -> removed
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", right - 35, right + 265, 1), ("MGA_2", 110, 150, 1)])
    comp.apply([[hit]])
    assert sorted(_cds_features(rec)) == ["MAP_1", "MGA_1"] and comp.counts["removed"] == 1


def test_annotate_removes_nothing(tmp_path):
    # the mapping only confirms MGA_1, so the overlapping same-strand MGA_2 is not removed (E. coli aceK case)
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 100, 100 + len(GENE), 1), ("MGA_2", 110, 150, 1)])
    comp.apply([[hit]])
    assert sorted(_cds_features(rec)) == ["MGA_1", "MGA_2"] and comp.counts["removed"] == 0

def test_no_start_with_same_gene_only_annotates(tmp_path):
    gene = "CCC" + GENE[3:]
    comp, rec, hit = _setup(tmp_path, gene, cds=[("MGA_1", 130, 100 + len(gene), 1)])
    comp.apply([[hit]])
    f = _cds_features(rec)
    assert list(f) == ["MGA_1"] and int(f["MGA_1"].location.start) == 130 and f["MGA_1"].primary_hit.id == "ARO:1"

def test_no_start_without_same_gene_is_pseudogene(tmp_path):
    comp, rec, hit = _setup(tmp_path, "CCC" + GENE[3:])
    comp.apply([[hit]])
    f = _cds_features(rec)["MAP_1"]
    assert f.qualifiers["pseudogene"] == ["unknown"] and "no start codon" in f.qualifiers["note"]

def test_invalid_reference_only_annotates(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE, validity=(False, True))
    comp.apply([[hit]])
    assert _cds_features(rec) == {}
    comp, rec, hit = _setup(tmp_path, GENE, validity=(False, True), cds=[("MGA_1", 130, 100 + len(GENE), 1)])
    comp.apply([[hit]])
    f = _cds_features(rec)["MGA_1"]
    assert int(f.location.start) == 130 and f.primary_hit.id == "ARO:1"

def test_cds_overlapping_two_regions_removed_once(tmp_path):
    gene2 = GENE + "C" * 30 + GENE
    comp, rec, hit1 = _setup(tmp_path, gene2, ref_len=len(GENE), cds=[("MGA_1", 110, 90 + len(gene2), 1)])
    hit2 = _hit(qlen=len(GENE), qend=len(GENE), sstart=101 + len(GENE) + 30, send=100 + len(gene2), slen=len(rec), btop=str(len(GENE)))
    comp.apply([[hit1], [hit2]])
    assert sorted(_cds_features(rec)) == ["MAP_1", "MAP_2"]

def test_second_instance_does_not_touch_mapped_feature(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    mapped = ExtendedFeature(location=FeatureLocation(100, 100 + len(GENE), strand=1), type="CDS", id="MAP_9", seq_id="sequence1",
                             qualifiers={"translation": ["M"]}, annotations={"partial_flag": "00", "nucl_mapping": "rescued"})
    mapped.primary_hit = NuclHit("VFG1", comp.nucl_refs["ARO:1"], "VFDB", 0.0, 1.0, 99.0, 100.0, 100.0, "")
    rec.features.append(mapped); comp.genome.set_feature_dictionary()
    comp.apply([[hit]])
    assert list(_cds_features(rec)) == ["MAP_9"] and [h.id for h in mapped.secondary_hits] == ["ARO:1"]


GENE_A = "ATG" + "GCT" * 40 + "TAA"                       # 126 bp
GENE_B = GENE_A[-35:] + "G" + "GCT" * 40 + "TAA"          # starts (CTG) 35 bp before the end of GENE_A
B0 = 100 + len(GENE_A) - 35                                 # start0 of GENE_B in the genome

def _two_genes(tmp_path, cds):
    """GENE_A at 100..226 and GENE_B overlapping it by 35 bp (ecpE/ecpD), each with its own full-length reference."""
    comp, rec, hit_a = _setup(tmp_path, GENE_A[:-35] + GENE_B, ref_len=len(GENE_A), cds=cds)
    _add_reference(comp, "ARO:2", GENE_B)
    hit_b = _hit(ref_id="ARO:2", qlen=len(GENE_B), qend=len(GENE_B), sstart=B0 + 1, send=B0 + len(GENE_B),
                 slen=len(rec), btop=str(len(GENE_B)))
    return comp, rec, hit_a, hit_b

def test_overlapping_neighbours_keep_their_cds(tmp_path):
    comp, rec, hit_a, hit_b = _two_genes(tmp_path, [("MGA_1", 100, 100 + len(GENE_A), 1),
                                                    ("MGA_2", B0, B0 + len(GENE_B), 1)])
    comp.apply([[hit_a], [hit_b]])
    f = _cds_features(rec)
    assert list(f) == ["MGA_1", "MGA_2"] and comp.counts["removed"] == 0
    assert f["MGA_1"].primary_hit.id == "ARO:1" and f["MGA_2"].primary_hit.id == "ARO:2"

def test_neighbour_of_rescued_gene_is_processed(tmp_path):
    comp, rec, hit_a, hit_b = _two_genes(tmp_path, [("MGA_2", B0, B0 + len(GENE_B), 1)])
    comp.apply([[hit_a], [hit_b]])
    f = _cds_features(rec)
    assert sorted(f) == ["MAP_1", "MGA_2"] and f["MGA_2"].primary_hit.id == "ARO:2"

def test_small_overlap_with_earlier_instance_feature_is_ignored(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    mapped = ExtendedFeature(location=FeatureLocation(100 + len(GENE) - 20, 250, strand=1), type="CDS", id="MAP_9",
                             seq_id="sequence1", qualifiers={"translation": ["M"]},
                             annotations={"partial_flag": "00", "nucl_mapping": "rescued"})
    rec.features.append(mapped); comp.genome.set_feature_dictionary()
    comp.apply([[hit]])
    f = _cds_features(rec)
    assert sorted(f) == ["MAP_10", "MAP_9"] and f["MAP_10"].primary_hit.id == "ARO:1" and not mapped.secondary_hits


def _second_instance(tmp_path, comp):
    """A VFDB instance on the genome of comp, with the same references."""
    vfdb = _component(tmp_path, comp.genome, db_name="VFDB")
    vfdb.nucl_refs, vfdb.ref_validity = comp.nucl_refs, comp.ref_validity
    return vfdb

def test_neighbours_mapped_by_two_instances(tmp_path):
    comp, rec, hit_a, hit_b = _two_genes(tmp_path, [("MGA_1", 100, 100 + len(GENE_A), 1),
                                                    ("MGA_2", B0, B0 + len(GENE_B), 1)])
    comp.apply([[hit_a]])                       # CARD: GENE_A only
    vfdb = _second_instance(tmp_path, comp)
    vfdb.apply([[hit_b]])                       # VFDB: GENE_B only
    f = _cds_features(rec)
    assert f["MGA_1"].primary_hit.id == "ARO:1" and comp.counts["removed"] == vfdb.counts["removed"] == 0
    # CARD only confirms MGA_1, so it removes nothing; VFDB then annotates MGA_2 in place.
    assert sorted(f) == ["MGA_1", "MGA_2"] and f["MGA_2"].primary_hit.id == "ARO:2"

UPSTREAM = "ATG" + "GCT" * 9                    # an in-frame start codon 30 bp upstream of GENE

def _longer_reference(tmp_path, cds):
    """Genome 100 C + UPSTREAM + GENE + 100 C; hit_short maps GENE (ARO:1), hit_long maps UPSTREAM + GENE (ARO:2)."""
    comp, rec, hit_short = _setup(tmp_path, UPSTREAM + GENE, cds=cds)
    hit_short = _hit(qlen=len(GENE), qend=len(GENE), sstart=131, send=130 + len(GENE), slen=len(rec), btop=str(len(GENE)))
    _add_reference(comp, "ARO:2", UPSTREAM + GENE)
    n = len(UPSTREAM + GENE)
    hit_long = _hit(ref_id="ARO:2", qlen=n, qend=n, sstart=101, send=100 + n, slen=len(rec), btop=str(n))
    return comp, rec, hit_short, hit_long

def test_later_instance_does_not_relocate_annotated_cds(tmp_path):
    end = 130 + len(GENE)
    comp, rec, hit_short, hit_long = _longer_reference(tmp_path, [("MGA_1", 130, end, 1)])
    comp.apply([[hit_short]])                   # CARD annotates MGA_1 in place
    vfdb = _second_instance(tmp_path, comp)
    vfdb.apply([[hit_long]])                    # same locus, different start: skipped
    f = _cds_features(rec)
    assert list(f) == ["MGA_1"] and (int(f["MGA_1"].location.start), int(f["MGA_1"].location.end)) == (130, end)
    assert f["MGA_1"].primary_hit.id == "ARO:1" and not f["MGA_1"].secondary_hits

def test_same_locus_different_location_of_rescued_feature_is_skipped(tmp_path):
    comp, rec, hit_short, hit_long = _longer_reference(tmp_path, [])
    mapped = ExtendedFeature(location=FeatureLocation(130, 130 + len(GENE), strand=1), type="CDS", id="MAP_9",
                             seq_id="sequence1", qualifiers={"translation": ["M"]},
                             annotations={"partial_flag": "00", "nucl_mapping": "rescued"})
    rec.features.append(mapped); comp.genome.set_feature_dictionary()
    comp.apply([[hit_long]])
    assert list(_cds_features(rec)) == ["MAP_9"] and not mapped.secondary_hits
    assert int(mapped.location.start) == 130



def test_write_hits_json(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    comp.apply([[hit]]); comp.write_hits()
    data = json.load(open(tmp_path / "nucl_hits.json"))
    assert data["MAP_1"]["model"]["note"].startswith("rescued; ") and data["MAP_1"]["id"] == "ARO:1"

def test_summarize_amr_reads_nucl_mapping(tmp_path):
    work = tmp_path / "NuclMapping"; work.mkdir()
    comp, rec, hit = _setup(work, GENE)
    comp.apply([[hit]]); comp.write_hits()
    rec.features[0].qualifiers["locus_tag"] = ["LOCUS_010"]
    summarize_amr(comp.genome, str(tmp_path), {}, str(tmp_path / "amr_summary.tsv"))
    rows = [l.split("\t") for l in open(tmp_path / "amr_summary.tsv") if not l.startswith("#")]
    assert rows[0][0] == "LOCUS_010" and rows[0][-1].startswith("rescued; ")

def test_summarize_pseudo_accepts_nucl_hit(tmp_path):
    gene = GENE[:10] + GENE[11:]
    comp, rec, hit = _setup(tmp_path, gene, ref_len=len(GENE), btop="10A-" + str(len(GENE) - 11))
    comp.apply([[hit]])
    rec.features[0].qualifiers["locus_tag"] = ["LOCUS_010"]
    summarize_pseudo(comp.genome, str(tmp_path / "pseudogene_summary.tsv"))
    rows = [l.split("\t") for l in open(tmp_path / "pseudogene_summary.tsv")]
    assert rows[1][1] == "LOCUS_010" and rows[1][3] == "frameshift" and rows[1][6] == "ARO:1"

def test_prepare_queries_skips_cds_without_translation(tmp_path):
    from dfc.components.baseComponent import BaseAnnotationComponent
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 0, 30, 1)])
    rec.features.append(new_pseudogene(MappedRegion("sequence1", 100, 100 + len(GENE), 1, "frameshift", "", [], set(), set()), "MAP_1", INFERENCE))
    comp.genome.set_feature_dictionary()
    base = BaseAnnotationComponent.__new__(BaseAnnotationComponent)
    base.genome, base.workDir, base.CPU, base.skipAnnotatedFeatures, base.logger = comp.genome, str(tmp_path), 1, False, comp.logger
    base.query_files, base.query_sequences = {}, {}
    base.prepareQueries(split_query=False)
    assert open(base.query_files[0]).read().count(">") == 1

def test_fix_partial_cds_keeps_complete_nucl_hit_cds(tmp_path):
    from tests.test_overlap_rules import _Genome, _Config, _feature
    from dfc.utils.feature_util import FeatureUtil
    comp, rec, hit = _setup(tmp_path, GENE)
    nucl_hit = NuclHit("ARO:1", comp.nucl_refs["ARO:1"], "CARD", 0.0, 500.0, 99.0, 100.0, 100.0, "")
    complete = _feature("CDS", 0, 66, "MAP_1", qualifiers={"product": ["x"], "translation": ["M"]})
    complete.annotations["partial_flag"] = "00"; complete.primary_hit = nucl_hit
    partial = _feature("CDS", 0, 66, "MGA_1", qualifiers={"product": ["x"], "translation": ["M"]})
    partial.annotations["partial_flag"] = "10"; partial.primary_hit = nucl_hit   # shorter than 500 bp: not fixable
    fu = FeatureUtil(_Genome([complete, partial]), _Config())
    fu.fix_partial_CDS(complete); fu.fix_partial_CDS(partial)
    assert complete.type == "CDS" and partial.type == "misc_feature"


def test_summarize_pseudo_keeps_empty_protein_description(tmp_path):
    from dfc.models.hit import ProteinHit
    comp, rec, _ = _setup(tmp_path, GENE)
    f = new_pseudogene(MappedRegion("sequence1", 100, 100 + len(GENE), 1, "frameshift", "", [], set(), {111}), "PG_1", INFERENCE)
    f.secondary_hits.append(PseudoGene("P1", set(), {111}))
    f.primary_hit = ProteinHit("P1", "", "", "", "", "", "DB", 0.0, 1.0, 99.0, 100.0, 100.0, "")
    f.seq_id = "sequence1"
    rec.features.append(f); comp.genome.set_feature_dictionary()
    summarize_pseudo(comp.genome, str(tmp_path / "pseudogene_summary.tsv"))
    rows = [l.rstrip("\n").split("\t") for l in open(tmp_path / "pseudogene_summary.tsv")]
    assert rows[1][3] == "frameshift" and rows[1][6] == "P1" and rows[1][7] == ""

def test_summarize_pseudo_lists_stop_lost_and_no_start(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE[:-3] + "GCT")        # C flank: no stop within 100 bp
    comp.apply([[hit]])
    comp2, rec2, hit2 = _setup(tmp_path, "CCC" + GENE[3:])
    comp2.apply([[hit2]])
    for c in (comp, comp2):
        summarize_pseudo(c.genome, str(tmp_path / f"{id(c)}.tsv"))
    r1 = [l.split("\t") for l in open(tmp_path / f"{id(comp)}.tsv")][1]
    r2 = [l.split("\t") for l in open(tmp_path / f"{id(comp2)}.tsv")][1]
    assert r1[3] == "stop_codon_lost" and r1[6] == "ARO:1" and r1[7] == "OXA beta-lactamase"
    assert r2[3] == "no_start_codon"

def test_summarize_amr_dedups_feature_in_two_json_files(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    comp.apply([[hit]])
    h = comp.hits_json["MAP_1"]
    for d, acc in (("NuclMapping", "ARO:1"), ("NuclMapping_2", "ARO:2")):
        (tmp_path / d).mkdir()
        json.dump({"MAP_1": dict(h, id=acc)}, open(tmp_path / d / "nucl_hits.json", "w"))
    summarize_amr(comp.genome, str(tmp_path), {}, str(tmp_path / "amr_summary.tsv"))
    rows = [l for l in open(tmp_path / "amr_summary.tsv") if not l.startswith("#")]
    assert len(rows) == 2

def test_pseudogene_note_kept_on_other_features():
    f = ExtendedFeature(location=FeatureLocation(100, 166, strand=1), type="CDS", id="MGA_1", seq_id="sequence1")
    f.secondary_hits.append(PseudoGene("P1", set(), {111}))     # from PseudoGeneDetection
    f.assign_hit()
    assert any(n.startswith("frameshifted") for n in f.qualifiers["note"])

def test_invalid_reference_attaches_to_in_frame_overlapping_cds(tmp_path):
    right = 100 + len(GENE)
    # MGA_1 shares the start and frame but not the stop; MGA_2 overlaps more but in another frame
    comp, rec, hit = _setup(tmp_path, GENE, validity=(False, True), cds=[("MGA_1", 100, right - 30, 1),
                                                                        ("MGA_2", 101, right + 1, 1)])
    comp.apply([[hit]])
    f = _cds_features(rec)
    assert list(f) == ["MGA_1", "MGA_2"] and comp.counts["removed"] == 0
    assert f["MGA_1"].primary_hit.id == "ARO:1" and f["MGA_2"].primary_hit is None
    assert (int(f["MGA_1"].location.start), int(f["MGA_1"].location.end)) == (100, right - 30)

def test_distant_defect_attaches_to_in_frame_overlapping_cds(tmp_path):
    gene = GENE[:10] + GENE[11:]
    comp, rec, hit = _setup(tmp_path, gene, ref_len=len(GENE), btop="10A-" + str(len(GENE) - 11),
                            cds=[("MGA_1", 100, 139, 1), ("MGA_2", 99, 100 + len(gene), -1)])
    hit.pident = 93.6
    comp.apply([[hit]])
    f = _cds_features(rec)
    assert sorted(f) == ["MGA_1", "MGA_2"] and comp.counts["removed"] == 0
    assert f["MGA_1"].primary_hit.id == "ARO:1" and f["MGA_2"].primary_hit is None

def test_pseudogene_called_at_cutoff_identity(tmp_path):
    gene = GENE[:10] + GENE[11:]
    comp, rec, hit = _setup(tmp_path, gene, ref_len=len(GENE), btop="10A-" + str(len(GENE) - 11),
                            cds=[("MGA_1", 100, 139, 1)])
    hit.pident = 95.0                                       # >= pseudogene_pident_cutoff
    comp.apply([[hit]])
    f = _cds_features(rec)
    assert list(f) == ["MAP_1"] and f["MAP_1"].qualifiers["pseudogene"] == ["unknown"]

def test_invalid_reference_ranks_after_valid_one(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    hit.pident = 95.0
    _add_reference(comp, "ARO:2", GENE, validity=(False, True))
    invalid = _hit(ref_id="ARO:2", pident=100.0, bitscore=600.0, qlen=len(GENE), qend=len(GENE), sstart=101,
                   send=100 + len(GENE), slen=len(rec), btop=str(len(GENE)))
    comp.apply([[invalid, hit]])
    f = _cds_features(rec)
    assert list(f) == ["MAP_1"] and f["MAP_1"].primary_hit.id == "ARO:1"
    assert [h.id for h in f["MAP_1"].secondary_hits] == ["ARO:2"]

def test_secondary_hit_keeps_hits_json_of_created_feature(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    _add_reference(comp, "ARO:2", GENE)
    other = _hit(ref_id="ARO:2", qlen=len(GENE), qend=len(GENE), sstart=101, send=100 + len(GENE), slen=len(rec),
                 btop=str(len(GENE)))
    comp.apply([[hit], [other]])                # a later group at the same location of the same instance
    assert list(_cds_features(rec)) == ["MAP_1"] and comp.hits_json["MAP_1"]["id"] == "ARO:1"
    assert comp.hits_json["MAP_1"]["model"]["note"].startswith("rescued; ")

def test_missing_reference_files_abort_with_hint(tmp_path, caplog):
    comp = _component(tmp_path, FakeGenome([]), database=str(tmp_path / "CARD"))
    try:
        comp.load_references()
    except SystemExit as e:
        assert e.code == 1
    else:
        raise AssertionError("no exit")
    assert "reference_util_for_nucl.py --card --vfdb" in caplog.text


# Truncated CDS (yddK case): a 1-bp insertion leaves an ORF of 85% of the reference from its start codon
LONG = "ATG" + "GCT" * 16 + "AAA" * 3 + "TAA"
LONG_INS = LONG[:51] + "T" + LONG[51:]


def test_late_frameshift_corrects_the_truncated_cds(tmp_path):
    comp, rec, hit = _setup(tmp_path, LONG_INS, ref_len=len(LONG), btop="51-T12", cds=[("MGA_1", 130, 154, 1)])
    comp.apply([[hit]])
    f = _cds_features(rec)
    assert list(f) == ["MGA_1"] and (int(f["MGA_1"].location.start), int(f["MGA_1"].location.end)) == (100, 154)
    assert "85% of the reference length; frameshift at around 152" in f["MGA_1"].qualifiers["note"]
    assert "pseudogene" not in f["MGA_1"].qualifiers
    assert comp.hits_json["MGA_1"]["model"]["note"].startswith("start corrected, truncated (85% of the reference); ")


def test_late_frameshift_without_cds_rescues_the_truncated_cds(tmp_path):
    comp, rec, hit = _setup(tmp_path, LONG_INS, ref_len=len(LONG), btop="51-T12")
    comp.apply([[hit]])
    f = _cds_features(rec)["MAP_1"]
    assert (int(f.location.start), int(f.location.end)) == (100, 154) and "pseudogene" not in f.qualifiers
    assert f.qualifiers["translation"] == ["M" + "A" * 16]


def test_truncated_fraction_option(tmp_path):
    comp, rec, hit = _setup(tmp_path, LONG_INS, ref_len=len(LONG), btop="51-T12")
    comp.truncated_fraction = 0.9            # 85% is not enough: pseudogene as before
    comp.apply([[hit]])
    assert _cds_features(rec)["MAP_1"].qualifiers["pseudogene"] == ["unknown"]


def test_custom_plasmid_db_reference_is_loaded(tmp_path, capsys):
    from dfc.models.nucref import PLASMID_DB
    ref = PLASMID_DB("P1", "repA", "replication protein", "", "", GENE, "", "AB000001", "pTest")
    (tmp_path / "PDB.nucl.ref").write_text(ref.to_tabular() + "\n")
    (tmp_path / "PDB.nucl.fasta").write_text(ref.to_nucl_fasta())
    comp = _component(tmp_path, FakeGenome([]), database=str(tmp_path / "PDB"), db_name="PLASMID_DB")
    comp.load_references()
    assert list(comp.nucl_refs) == ["P1"] and comp.ref_validity["P1"] == (True, True)
    assert comp.nucl_refs["P1"].to_dict()["accession"] == "PLASMID_DB:P1" and capsys.readouterr().out == ""


def test_truncated_cds_is_reported_when_only_annotated(tmp_path):
    comp, rec, hit = _setup(tmp_path, LONG_INS, ref_len=len(LONG), btop="51-T12", cds=[("MGA_1", 100, 154, 1)])
    comp.apply([[hit]])
    assert comp.hits_json["MGA_1"]["model"]["note"].startswith("truncated (85% of the reference); ")


def test_truncated_cds_needs_pseudogene_identity(tmp_path):
    # 94% < pseudogene_pident_cutoff: no truncated CDS; the hit is only attached to the in-frame predicted CDS
    comp, rec, hit = _setup(tmp_path, LONG_INS, ref_len=len(LONG), btop="51-T12", cds=[("MGA_1", 130, 154, 1)])
    hit.pident = 94.0
    comp.apply([[hit]])
    f = _cds_features(rec)
    assert list(f) == ["MGA_1"] and int(f["MGA_1"].location.start) == 130 and f["MGA_1"].primary_hit.id == "ARO:1"
    assert comp.counts["start corrected"] == comp.counts["rescued"] == comp.counts["removed"] == 0


def test_intact_reference_ranks_before_truncated_one(tmp_path):
    # ARO:1 (LONG, 99.9%) is truncated by the insertion; ARO:2 is the intact 54-bp ORF (99.0%) at the same locus
    comp, rec, hit_long = _setup(tmp_path, LONG_INS, ref_len=len(LONG), btop="51-T12")
    hit_long.pident = 99.9
    orf = LONG_INS[:54]
    _add_reference(comp, "ARO:2", orf)
    hit_orf = _hit(ref_id="ARO:2", qlen=54, qend=54, sstart=101, send=154, slen=len(rec), btop="54", pident=99.0)
    comp.apply([[hit_long, hit_orf]])
    f = _cds_features(rec)["MAP_1"]
    assert f.primary_hit.id == "ARO:2" and not [n for n in f.qualifiers["note"] if "of the reference length" in n]
