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
    comp.result_file = str(tmp_path / "mapping.out")
    comp.references, comp.ref_validity, comp.hits_json = {}, {}, {}
    for k, v in options.items(): setattr(comp, k, v)
    return comp

def test_registered():
    assert COMPONENTS["NuclMapping"] is NuclMapping

def test_select_hits_keeps_full_length_hits_on_known_contigs(tmp_path, caplog):
    genome = FakeGenome([SeqRecord(Seq("A" * 5000), id="sequence1")])
    comp = _component(tmp_path, genome)
    rows = ["ARO:1\tsequence1\t99.0\t300\t3\t0\t1\t300\t1001\t1300\t0.0\t500\t300\t5000\t300",
            "ARO:2\tsequence1\t99.0\t200\t2\t0\t50\t250\t2001\t2200\t0.0\t300\t300\t5000\t200",   # partial
            "ARO:3\tunknown\t99.0\t300\t3\t0\t1\t300\t1001\t1300\t0.0\t500\t300\t5000\t300"]      # unknown contig
    (tmp_path / "mapping.out").write_text("\n".join(rows) + "\n")
    groups = comp.select_hits()
    assert [[h.ref_id for h in g] for g in groups] == [["ARO:1"]]
    assert "unknown" in caplog.text

def test_blast_command(tmp_path):
    comp = _component(tmp_path, FakeGenome([]), database="/db/card/CARD", evalue_cutoff=1e-10, aligner_options="", CPU=4)
    cmd = comp.blast_command("/w/genome")
    assert cmd[:3] == ["blastn", "-task", "dc-megablast"] and "/db/card/CARD.nucl.fasta" in cmd
    assert cmd[cmd.index("-perc_identity") + 1] == "90" and cmd[cmd.index("-num_threads") + 1] == "4"


def _setup(tmp_path, gene, cds=(), validity=(True, True), btop=None, ref_len=None):
    """One contig: 100 C + gene + 100 C, with predicted CDSs given as (id, start0, end, strand)."""
    rec = SeqRecord(Seq("C" * 100 + gene + "C" * 100), id="sequence1")
    rec.features = [ExtendedFeature(location=FeatureLocation(s, e, strand=st), type="CDS", id=i, seq_id="sequence1",
                                    qualifiers={"translation": ["M"], "transl_table": ["11"], "codon_start": [1]},
                                    annotations={"partial_flag": "00"}) for i, s, e, st in cds]
    comp = _component(tmp_path, FakeGenome([rec]))
    comp.references = {"ARO:1": CARD("ARO:1", "OXA-48", "OXA beta-lactamase", "P1", "", "N1", gene, "E. coli", "")}
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
    assert comp.genome.features["MAP_1"] is f

def test_annotate_keeps_identical_cds(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 100, 100 + len(GENE), 1)])
    comp.apply([[hit]])
    assert list(_cds_features(rec)) == ["MGA_1"]
    f = rec.features[0]
    assert f.primary_hit.id == "ARO:1" and "nucl_mapping" not in f.annotations
    assert not comp.hits_json["MGA_1"]["model"]["note"].startswith(("rescued", "start corrected", "pseudogene"))

def test_start_correction_relocates_same_gene(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 130, 100 + len(GENE), 1)])
    comp.apply([[hit]])
    f = _cds_features(rec)["MGA_1"]
    assert int(f.location.start) == 100 and f.qualifiers["translation"] == ["M" + "A" * 20]
    assert "start codon moved 30 bp upstream by mapping of CARD:OXA-48" in f.qualifiers["note"]
    assert f.annotations["nucl_mapping"] == "start corrected" and f.primary_hit.id == "ARO:1"
    assert comp.hits_json["MGA_1"]["model"]["note"].startswith("start corrected; ")

def test_frameshift_makes_pseudogene_and_removes_fragments(tmp_path):
    gene = GENE[:10] + GENE[11:]   # 1-bp deletion at position 111 of the genome
    comp, rec, hit = _setup(tmp_path, gene, ref_len=len(GENE), btop="10A-" + str(len(GENE) - 11),
                            cds=[("MGA_1", 100, 140, 1), ("MGA_2", 140, 100 + len(gene), 1)])
    comp.apply([[hit]])
    assert list(_cds_features(rec)) == ["MAP_1"]
    f = rec.features[0]
    assert f.qualifiers["pseudogene"] == ["unknown"] and "translation" not in f.qualifiers
    assert "frameshift at around 111" in f.qualifiers["note"]
    assert any(isinstance(h, PseudoGene) for h in f.secondary_hits) and f.primary_hit.id == "ARO:1"
    assert comp.hits_json["MAP_1"]["model"]["note"].startswith("pseudogene: frameshift; ")

def test_overlapping_cds_other_strand_removed_neighbour_kept(tmp_path):
    right = 100 + len(GENE)
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 110, right - 10, -1), ("MGA_2", right - 20, right + 100, 1)])
    comp.apply([[hit]])
    assert sorted(_cds_features(rec)) == ["MAP_1", "MGA_2"]

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
    comp, rec, hit1 = _setup(tmp_path, gene2, ref_len=len(GENE), cds=[("MGA_1", 100, 100 + len(gene2), -1)])
    hit2 = _hit(qlen=len(GENE), qend=len(GENE), sstart=101 + len(GENE) + 30, send=100 + len(gene2), slen=len(rec), btop=str(len(GENE)))
    comp.apply([[hit1], [hit2]])
    assert sorted(_cds_features(rec)) == ["MAP_1", "MAP_2"]

def test_second_instance_does_not_touch_mapped_feature(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    mapped = ExtendedFeature(location=FeatureLocation(100, 100 + len(GENE), strand=1), type="CDS", id="MAP_9", seq_id="sequence1",
                             qualifiers={"translation": ["M"]}, annotations={"partial_flag": "00", "nucl_mapping": "rescued"})
    mapped.primary_hit = NuclHit("VFG1", comp.references["ARO:1"], "VFDB", 0.0, 1.0, 99.0, 100.0, 100.0, "")
    rec.features.append(mapped); comp.genome.set_feature_dictionary()
    comp.apply([[hit]])
    assert list(_cds_features(rec)) == ["MAP_9"] and [h.id for h in mapped.secondary_hits] == ["ARO:1"]



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
    rec.features.append(new_pseudogene(MappedRegion("sequence1", 100, 100 + len(GENE), 1, "frameshift", "", [], set(), set()), "MAP_1", 11))
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
    nucl_hit = NuclHit("ARO:1", comp.references["ARO:1"], "CARD", 0.0, 500.0, 99.0, 100.0, 100.0, "")
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
    f = new_pseudogene(MappedRegion("sequence1", 100, 100 + len(GENE), 1, "frameshift", "", [], set(), {111}), "PG_1", 11)
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
