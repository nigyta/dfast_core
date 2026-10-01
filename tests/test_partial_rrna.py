# coding: UTF8
from Bio.SeqFeature import FeatureLocation, ExactPosition, BeforePosition, AfterPosition
from dfc.models.bio_feature import ExtendedFeature
from dfc.tools.barrnap import Barrnap
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dfc.genome import Genome
from dfc.utils.feature_util import FeatureUtil, adjust_rrna

ALIGNED_42 = "aligned only 42 percent of the 23S ribosomal RNA"


def _rrna(start, end, strand=1, note=None, location=None):
    qualifiers = {"product": ["23S ribosomal RNA"], "inference": ["COORDINATES:profile:Barrnap:0.8"]}
    if note:
        qualifiers["note"] = [note]
    return ExtendedFeature(location=location or FeatureLocation(ExactPosition(start), ExactPosition(end), strand),
                           type="rRNA", id="Barrnap_1", qualifiers=qualifiers)


def _gap(start, end):
    return ExtendedFeature(location=FeatureLocation(start, end, 1), type="assembly_gap", id="GAP_1")


def _coords(f):
    return int(f.location.start) + 1, int(f.location.end)


def test_partial_rrna_overlapping_gap_is_trimmed_and_converted():
    # rRNA 1574..2951 / assembly_gap 2948..3046 (1-based) -> misc_feature 1574..>2947
    f = _rrna(1573, 2951, note=ALIGNED_42)
    [p], _ = adjust_rrna(f, 4221, [_gap(2947, 3046)])
    assert p is f
    assert _coords(p) == (1574, 2947)
    assert isinstance(p.location.end, AfterPosition) and isinstance(p.location.start, ExactPosition)
    assert p.type == "misc_feature"
    assert "product" not in p.qualifiers
    assert p.qualifiers["note"] == ["putative rRNA, " + ALIGNED_42, "truncated at an assembly gap"]
    assert p.qualifiers["inference"] == ["COORDINATES:profile:Barrnap:0.8"]


def test_partial_rrna_not_truncated_is_converted_only():
    f = _rrna(1000, 2000, strand=-1, note=ALIGNED_42)
    [p], _ = adjust_rrna(f, 5000, [])
    assert _coords(p) == (1001, 2000) and p.location.strand == -1
    assert p.type == "misc_feature"
    assert p.qualifiers["note"] == ["putative rRNA, " + ALIGNED_42]


def test_complete_rrna_trimmed_at_gap_stays_rrna_with_partial_location():
    f = _rrna(3040, 4195, strand=-1)
    [p], _ = adjust_rrna(f, 5000, [_gap(2947, 3046)])
    assert _coords(p) == (3047, 4195)
    assert isinstance(p.location.start, BeforePosition) and isinstance(p.location.end, ExactPosition)
    assert p.type == "rRNA" and p.qualifiers["product"] == ["23S ribosomal RNA"]
    assert p.qualifiers["note"] == ["truncated at an assembly gap"]


def test_complete_rrna_near_contig_end_is_not_extended():
    f = _rrna(3, 1405)  # complete (>= 80%): starts at 1-based position 4
    [p], _ = adjust_rrna(f, 4221, [])
    assert _coords(p) == (4, 1405) and isinstance(p.location.start, ExactPosition)
    assert p.type == "rRNA" and "note" not in p.qualifiers


def test_partial_rrna_near_contig_end_is_extended_and_converted():
    f = _rrna(3, 1405, note=ALIGNED_42)
    [p], _ = adjust_rrna(f, 4221, [])
    assert _coords(p) == (1, 1405) and isinstance(p.location.start, BeforePosition)
    assert p.type == "misc_feature"
    assert p.qualifiers["note"] == ["putative rRNA, " + ALIGNED_42, "truncated at the contig end"]


def test_contig_end_margin_is_less_than_10_on_both_sides():
    # 9 bases before the end -> extended, 10 bases -> not extended
    [p], _ = adjust_rrna(_rrna(500, 4221 - 9, note=ALIGNED_42), 4221, [])
    assert int(p.location.end) == 4221 and isinstance(p.location.end, AfterPosition)
    [p], _ = adjust_rrna(_rrna(500, 4221 - 10, note=ALIGNED_42), 4221, [])
    assert int(p.location.end) == 4211 and isinstance(p.location.end, ExactPosition)
    [p], _ = adjust_rrna(_rrna(9, 1000, note=ALIGNED_42), 4221, [])
    assert int(p.location.start) == 0 and isinstance(p.location.start, BeforePosition)
    [p], _ = adjust_rrna(_rrna(10, 1000, note=ALIGNED_42), 4221, [])
    assert int(p.location.start) == 10 and isinstance(p.location.start, ExactPosition)


def test_end_exceeding_contig_is_clamped():
    f = _rrna(2500, 4300)
    [p], _ = adjust_rrna(f, 4221, [])
    assert int(p.location.end) == 4221 and isinstance(p.location.end, AfterPosition)
    assert p.type == "rRNA" and p.qualifiers["note"] == ["truncated at the contig end"]


def test_barrnap_contig_end_marker_kept_for_other_tools():
    loc = FeatureLocation(BeforePosition(0), ExactPosition(1405), 1)
    [p], _ = adjust_rrna(_rrna(0, 0, location=loc, note=ALIGNED_42), 4221, [])
    assert isinstance(p.location.start, BeforePosition) and p.type == "misc_feature"
    assert p.qualifiers["note"] == ["putative rRNA, " + ALIGNED_42, "truncated at the contig end"]


def test_gap_inside_splits_into_misc_features():
    # complete 23S spanning a scaffold gap: 1574..4497 with gap 2948..3045
    f = _rrna(1573, 4497)
    (left, right), _ = adjust_rrna(f, 4600, [_gap(2947, 3045)])
    assert left is f and right.id == "Barrnap_1_2"
    assert _coords(left) == (1574, 2947) and isinstance(left.location.end, AfterPosition)
    assert _coords(right) == (3046, 4497) and isinstance(right.location.start, BeforePosition)
    for p in (left, right):
        assert p.type == "misc_feature" and "product" not in p.qualifiers
        assert p.qualifiers["note"] == ["putative rRNA overlapping an assembly gap, 23S ribosomal RNA",
                                        "truncated at an assembly gap"]


def test_partial_rrna_with_gap_inside_splits_too():
    f = _rrna(1573, 4497, note=ALIGNED_42)
    (left, right), _ = adjust_rrna(f, 4600, [_gap(2947, 3045)])
    assert left.qualifiers["note"][0] == "putative rRNA overlapping an assembly gap, " + ALIGNED_42
    assert right.type == "misc_feature"


def test_complete_rrna_untouched():
    f = _rrna(100, 1600)
    [p], _ = adjust_rrna(f, 5000, [_gap(2000, 2100)])
    assert _coords(p) == (101, 1600) and isinstance(p.location.start, ExactPosition)
    assert p.type == "rRNA" and "note" not in p.qualifiers


def test_barrnap_reject_default():
    Barrnap.version = "0.8"  # skip the external version check
    assert Barrnap(options={}).cmd_options == "--reject 0.5"
    assert Barrnap(options={"cmd_options": "--lencutoff 0.6"}).cmd_options == "--reject 0.5 --lencutoff 0.6"
    assert Barrnap(options={"cmd_options": "--reject 0.25"}).cmd_options == "--reject 0.25"


def test_rrna_entirely_inside_gap_is_removed():
    kept, removed = adjust_rrna(_rrna(450, 550), 2000, [_gap(400, 600)])
    assert kept == [] and len(removed) == 1 and "inside an assembly gap" in removed[0]


def test_short_piece_after_trimming_is_removed():
    # left piece 381..400 (20 bp) is dropped, right piece is kept
    kept, removed = adjust_rrna(_rrna(380, 1500), 2000, [_gap(400, 600)])
    assert [_coords(p) for p in kept] == [(601, 1500)]
    assert len(removed) == 1 and "381..400" in removed[0] and "shorter than 30 bp" in removed[0]
    # an untrimmed feature is never removed for being short
    kept, removed = adjust_rrna(_rrna(100, 120), 2000, [])
    assert len(kept) == 1 and removed == []


# ---- integration through FeatureUtil.execute(): specific rules run before resolve_overlap() ----

class _Genome:
    sort_features = Genome.sort_features
    set_feature_dictionary = Genome.set_feature_dictionary

    def __init__(self, features, length=4000):
        record = SeqRecord(Seq("A" * length), id="seq1")
        for f in features:
            f.seq_id = "seq1"
        record.features = features
        self.seq_records = {"seq1": record}
        self.set_feature_dictionary()


class _Config:
    FEATURE_ADJUSTMENT = {}


def _cds(start, end, fid):
    return ExtendedFeature(location=FeatureLocation(start, end, 1), type="CDS", id=fid)


def _run(features, length=4000):
    genome = _Genome(features, length)
    FeatureUtil(genome, _Config()).execute()
    return genome.seq_records["seq1"].features


def test_execute_splits_rrna_overlapping_gap_by_more_than_10_percent():
    # Codex review case: rRNA [0:1500] with gap [400:600] was removed whole by resolve_overlap
    features = _run([_rrna(0, 1500), _gap(400, 600)])
    misc = [f for f in features if f.type == "misc_feature"]
    assert [(int(f.location.start), int(f.location.end)) for f in misc] == [(0, 400), (600, 1500)]
    assert [f.type for f in features] == ["misc_feature", "assembly_gap", "misc_feature"]  # sorted


def test_execute_removes_rrna_inside_gap_and_logs_summary(caplog):
    caplog.set_level("DEBUG")
    features = _run([_rrna(450, 550), _gap(400, 600)])
    assert [f.type for f in features] == ["assembly_gap"]
    assert "Removed 1 features by specific overlap rules" in caplog.text
    assert "Removed feature by overlap rule" in caplog.text  # per-feature line at DEBUG


def test_execute_putative_rrna_still_masks_cds():
    # CDS overlapping a partial rRNA (now misc_feature) by >10% is still removed by the fallback
    features = _run([_rrna(1000, 2000, note=ALIGNED_42), _cds(1500, 2400, "CDS_1")])
    assert [f.type for f in features] == ["misc_feature"]


def test_execute_fallback_removal_is_a_warning(caplog):
    features = _run([_rrna(1000, 2000), _cds(1500, 2400, "CDS_1")])
    assert [f.type for f in features] == ["rRNA"]
    assert any(r.levelname == "WARNING" and "not handled by specific overlap rules" in r.message
               for r in caplog.records)
