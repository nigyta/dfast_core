# coding: UTF8
"""Overlap handling through FeatureUtil: specific rules first, resolve_overlap() as the fallback."""
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, ExactPosition
from Bio.SeqRecord import SeqRecord
from dfc.genome import Genome
from dfc.models.bio_feature import ExtendedFeature
from dfc.models.hit import ProteinHit
from dfc.utils.feature_util import FeatureUtil

ALIGNED_42 = "aligned only 42 percent of the 23S ribosomal RNA"


class _Genome:
    sort_features = Genome.sort_features
    set_feature_dictionary = Genome.set_feature_dictionary

    def __init__(self, features, length=6000):
        record = SeqRecord(Seq("A" * length), id="seq1")
        for f in features:
            f.seq_id = "seq1"
        record.features = features
        self.seq_records = {"seq1": record}
        self.set_feature_dictionary()


class _Config:
    FEATURE_ADJUSTMENT = {}


def _feature(type_, start, end, fid, strand=1, qualifiers=None):
    return ExtendedFeature(location=FeatureLocation(ExactPosition(start), ExactPosition(end), strand),
                           type=type_, id=fid, qualifiers=qualifiers or {})


def _rrna(start, end, fid="Barrnap_1", strand=1, note=None):
    qualifiers = {"product": ["5S ribosomal RNA"]}
    if note:
        qualifiers["note"] = [note]
    return _feature("rRNA", start, end, fid, strand, qualifiers)


def _trna(start, end, fid="Aragorn_1", strand=1):
    return _feature("tRNA", start, end, fid, strand, {"product": ["tRNA-Gly"]})


def _cds(start, end, fid="MGA_1", strand=1, product=None):
    cds = _feature("CDS", start, end, fid, strand, {"product": ["hypothetical protein"]})
    if product:
        cds.primary_hit = ProteinHit("WP_000001.1", product, gene="", ec_number="", source_db="RefSeq", organism="",
                                     db_name="test", e_value=1e-50, score=500, identity=90.0, q_cov=90.0, s_cov=90.0, flag="")
    return cds


def _gap(start, end):
    return _feature("assembly_gap", start, end, "GAP_1")


def _run(features, length=6000):
    genome = _Genome(features, length)
    fu = FeatureUtil(genome, _Config())
    fu.execute()
    fu.execute_after_annotation()
    return [f.id for f in genome.seq_records["seq1"].features], genome


# ---- rRNA vs assembly_gap (adjust_rrna_features) ----

def test_rrna_overlapping_gap_by_more_than_10_percent_is_split_not_removed():
    # Codex review case: rRNA [0:1500] with gap [400:600] used to be removed whole by resolve_overlap
    ids, genome = _run([_rrna(0, 1500), _gap(400, 600)])
    features = genome.seq_records["seq1"].features
    assert [f.type for f in features] == ["misc_feature", "assembly_gap", "misc_feature"]
    assert [(int(f.location.start), int(f.location.end)) for f in features if f.type == "misc_feature"] == [(0, 400), (600, 1500)]


def test_rrna_inside_gap_is_removed_and_logged(caplog):
    caplog.set_level("DEBUG")
    ids, _ = _run([_rrna(450, 550), _gap(400, 600)])
    assert ids == ["GAP_1"]
    assert "Removed 1 features by specific overlap rules" in caplog.text
    assert "Removed feature by overlap rule" in caplog.text  # per-feature line at DEBUG


# ---- rRNA vs rRNA (ANN5310) ----

def test_shorter_of_overlapping_rrnas_is_removed():
    ids, _ = _run([_rrna(1000, 1120, "Barrnap_1"), _rrna(1100, 2600, "Barrnap_2")])
    assert ids == ["Barrnap_2"]


def test_rrnas_on_opposite_strands_are_kept():
    ids, _ = _run([_rrna(1000, 1120, "Barrnap_1"), _rrna(1100, 2600, "Barrnap_2", strand=-1)])
    assert ids == ["Barrnap_1", "Barrnap_2"]


# ---- rRNA / tRNA vs CDS (ANN5310 / ANN5320), decided by CDS products ----

def test_hypothetical_cds_overlapping_rrna_by_one_base_is_removed():
    # SAMD01959111 case: 5S rRNA 219..329 vs hypothetical CDS 324..782 (6 bp overlap)
    ids, _ = _run([_rrna(218, 329), _cds(323, 782)])
    assert ids == ["Barrnap_1"]


def test_rrna_overlapping_functional_cds_is_removed():
    ids, _ = _run([_rrna(218, 329), _cds(323, 782, product="DNA polymerase III subunit alpha")])
    assert ids == ["MGA_1"]


def test_trna_in_hypothetical_cds_removes_cds():
    ids, _ = _run([_cds(4004, 5432), _trna(5041, 5129)])
    assert ids == ["Aragorn_1"]


def test_trna_in_functional_cds_is_removed():
    # SAMD01959133 case: tRNA-Pro inside an MFS transporter
    ids, _ = _run([_cds(4004, 5432, product="MFS transporter"), _trna(5041, 5129)])
    assert ids == ["MGA_1"]


def test_hypothetical_product_names_count_as_hypothetical():
    ids, _ = _run([_cds(4004, 5432, product="Conserved hypothetical protein"), _trna(5041, 5129)])
    assert ids == ["Aragorn_1"]


def test_trna_partially_overlapping_cds_is_not_a_rule_case():
    # not contained: ANN5320 does not apply; the overlap (9 of 1000 bp) is below the fallback threshold
    ids, _ = _run([_cds(5120, 6000, product="MFS transporter"), _trna(5041, 5129)])
    assert ids == ["Aragorn_1", "MGA_1"]


def test_opposite_strand_containment_is_left_to_fallback(caplog):
    # validator allows opposite strands; the fallback still removes a CDS masked by >10%
    ids, _ = _run([_cds(5000, 5600, strand=-1), _trna(5041, 5129)])
    assert ids == ["Aragorn_1"]
    assert any(r.levelname == "WARNING" and "not handled by specific overlap rules" in r.message for r in caplog.records)


def test_rna_overlapping_hypothetical_and_functional_cds_keeps_both_cds():
    ids, _ = _run([_cds(100, 400, "MGA_1"), _rrna(300, 1500), _cds(1400, 2000, "MGA_2", product="ribosomal protein L2")])
    assert ids == ["MGA_1", "MGA_2"]


def test_putative_rrna_misc_feature_still_masks_cds_in_fallback():
    # partial rRNA -> misc_feature is not checked by ANN5310 but still masks CDS as rRNA in resolve_overlap
    ids, _ = _run([_rrna(1000, 2000, note=ALIGNED_42), _cds(1500, 2400)])
    assert ids == ["Barrnap_1"]


# ---- eukaryotic rRNA (possible contamination) vs CDS ----

def _euk18s(start, end, fid="Barrnap_1", strand=1):
    f = _feature("rRNA", start, end, fid, strand, {"product": ["18S ribosomal RNA"]})
    f.annotations["eukaryotic_rrna"] = True
    return f


def test_eukaryotic_rrna_becomes_misc_and_removes_hypothetical_cds_on_both_strands():
    ids, genome = _run([_euk18s(100, 1970), _cds(200, 434, "MGA_1"), _cds(1000, 1234, "MGA_2", strand=-1),
                        _cds(1950, 2500, "MGA_3", product="DNA polymerase"), _cds(3000, 3300, "MGA_4")])
    assert ids == ["Barrnap_1", "MGA_3", "MGA_4"]  # a functional CDS is kept (overlap < 10%, else resolve_overlap removes it)
    misc = genome.seq_records["seq1"].features[0]
    assert misc.type == "misc_feature"
    assert misc.qualifiers["note"] == ["eukaryotic 18S ribosomal RNA-like sequence, possible contamination"]
