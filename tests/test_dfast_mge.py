# coding: UTF8
import pytest
from dfc.tools.dfast_mge import classify_mge


@pytest.mark.parametrize("type_code,evidence,expected", [
    # (feature_key, mobile_element_type, is_putative)
    ("is", 1, ("mobile_element", "insertion sequence", False)),
    ("mite", 1, ("mobile_element", "MITE", False)),
    ("tn", 1, ("mobile_element", "transposon", False)),
    ("cn", 1, ("mobile_element", "transposon", False)),       # DB ヒット
    ("cn", 2, ("misc_feature", None, True)),                  # putative
    ("integron", 1, ("mobile_element", "integron", False)),
    ("retrotransposon", 1, ("mobile_element", "retrotransposon", False)),
    ("ice", 1, ("misc_feature", None, False)),
    ("ime", 1, ("misc_feature", None, False)),
    ("cime", 1, ("misc_feature", None, False)),
    ("mic", 1, ("misc_feature", None, False)),
    ("unknown_code", 1, ("misc_feature", None, False)),       # 未知コードは misc_feature
])
def test_classify_mge(type_code, evidence, expected):
    assert classify_mge(type_code, evidence) == expected


from dfc.tools.dfast_mge import build_qualifiers


def _entry(**kw):
    base = dict(name="Tn1999", type="cn", evidence=1, identity=1.0,
                coverage=1.0, template={"accession": "AY236073"})
    base.update(kw)
    return base


def test_qualifiers_mobile_element():
    q = build_qualifiers(_entry())
    assert q["mobile_element_type"] == ["transposon:Tn1999"]
    assert any("MobileElementFinder: Tn1999" in n for n in q["note"])
    assert any("identity:100.0%" in n and "coverage:100.0%" in n for n in q["note"])
    assert any("AY236073" in n for n in q["note"])


def test_qualifiers_misc_putative():
    q = build_qualifiers(_entry(evidence=2))
    assert "mobile_element_type" not in q
    assert any("putative composite transposon" in n for n in q["note"])


def test_qualifiers_misc_ice():
    q = build_qualifiers(_entry(name="SXT", type="ice", evidence=1,
                                template={"accession": "AY055428"}))
    assert "mobile_element_type" not in q
    assert any("integrative conjugative element (ICE)" in n for n in q["note"])


def test_qualifiers_omit_empty_accession():
    # accession が空のときは "similar to  (MGEdb)" の空欄を出さない
    q = build_qualifiers(_entry(template={"accession": ""}))
    assert all("similar to" not in n for n in q["note"])
    assert all("(MGEdb)" not in n for n in q["note"])


from Bio.SeqFeature import BeforePosition, AfterPosition, ExactPosition
from dfc.tools.dfast_mge import entry_to_feature


def _full_entry(**kw):
    base = dict(name="Tn1999", type="cn", evidence=1, identity=1.0, coverage=1.0,
                template={"accession": "AY236073"},
                contig="JN626286.1 Klebsiella pneumoniae plasmid pOXA-48",
                start=2292, end=8429, strand=-1, trunc_5p=0, trunc_3p=0)
    base.update(kw)
    return base


def test_entry_to_feature_basic():
    seq_id, feat = entry_to_feature(_full_entry(), index=1)
    assert seq_id == "JN626286.1"                      # defline 先頭トークン
    assert feat.type == "mobile_element"
    assert int(feat.location.start) == 2291            # 1-based 2292 -> 0-based 2291
    assert int(feat.location.end) == 8429
    assert feat.location.strand == -1
    assert feat.id == "MGE_1"
    assert feat.qualifiers["mobile_element_type"] == ["transposon:Tn1999"]


def test_entry_to_feature_putative_is_misc():
    seq_id, feat = entry_to_feature(_full_entry(evidence=2), index=2)
    assert feat.type == "misc_feature"
    assert feat.id == "MGE_2"


def test_entry_to_feature_complete_not_partial():
    # trunc_5p=1 は 1-based の alignment 開始位置 = 5' 端まで完全。partial にしない。
    seq_id, feat = entry_to_feature(_full_entry(strand=1, trunc_5p=1, trunc_3p=0), index=1)
    assert isinstance(feat.location.start, ExactPosition)
    assert isinstance(feat.location.end, ExactPosition)


def test_entry_to_feature_truncation_partial():
    # trunc_5p>1 で 5' truncated。strand=+1 なら左端(start)が BeforePosition、右端は Exact。
    seq_id, feat = entry_to_feature(_full_entry(strand=1, trunc_5p=50, trunc_3p=0), index=1)
    assert isinstance(feat.location.start, BeforePosition)
    assert isinstance(feat.location.end, ExactPosition)


def test_entry_to_feature_truncation_complement():
    # strand=-1 で trunc_5p>1 なら 5' は右端 -> AfterPosition、左端は Exact。
    seq_id, feat = entry_to_feature(_full_entry(strand=-1, trunc_5p=50, trunc_3p=0), index=1)
    assert isinstance(feat.location.end, AfterPosition)
    assert isinstance(feat.location.start, ExactPosition)


# ---- 座標系判定 (_zero_based_left) : MEF composite の off-by-one 吸収 ----
from dfc.tools.dfast_mge import _zero_based_left, _location


@pytest.mark.parametrize("start,end,asl,expected", [
    # IS など 1-based: end-start+1 == allele_seq_length -> left = start-1
    (6350896, 6353278, 2383, 6350895),
    # DBヒット composite (cn evidence=1) も 1-based (pOXA48 Tn1999 実データ)
    (2292, 8429, 6138, 2291),
    # putative composite (cn evidence=2) は 0-based: end-start == asl -> left = start (無補正)
    (2979776, 2986173, 6397, 2979776),
    # contig 先頭にかかる putative composite: start=0 -> left=0 (<0 にしない)
    (0, 1811, 1811, 0),
    # 将来 MEF が composite を 1-based に修正した場合: end-start+1==asl -> -1 が復活し二重補正なし
    (100, 200, 101, 99),
    # allele_seq_length 無し -> 文書化された 1-based を仮定
    (500, 600, None, 499),
    # どちらにも一致しない不整合 -> 1-based 仮定
    (500, 600, 55, 499),
])
def test_zero_based_left(start, end, asl, expected):
    assert _zero_based_left(start, end, asl) == expected


def test_putative_composite_zero_based_not_shifted():
    # putative composite の 0-based start はシフトされず IS と同じ基準に揃う
    e = _full_entry(evidence=2, start=2979776, end=2986173, allele_seq_length=6397,
                    trunc_5p=1, trunc_3p=0)
    _seq_id, feat = entry_to_feature(e, index=1)
    assert feat.type == "misc_feature"
    assert int(feat.location.start) == 2979776   # start-1 されない
    assert int(feat.location.end) == 2986173


def test_boundary_composite_not_negative():
    # contig 先頭の putative composite: <0 でなく <1 (BeforePosition(0)) になり None 化しない
    e = _full_entry(evidence=2, strand=1, start=0, end=1811, allele_seq_length=1811,
                    trunc_5p=50, trunc_3p=0, contig="sequence129")
    _seq_id, feat = entry_to_feature(e, index=33)
    assert feat.location is not None
    assert int(feat.location.start) == 0
    assert isinstance(feat.location.start, BeforePosition)


def test_clamp_guards_negative_left_without_length():
    # allele_seq_length 欠落 + start=0 では left=-1 になるが 0 にクランプされる (防御)
    loc = _location(0, 500, 1, trunc_5p=1, trunc_3p=0, allele_seq_length=None)
    assert int(loc.start) == 0
    assert isinstance(loc.start, BeforePosition)


# ---- putative composite の内部マーカーと DDBJ ann からの除外 ----

def test_entry_to_feature_marks_putative_composite():
    _sid, feat = entry_to_feature(
        _full_entry(evidence=2, start=2979776, end=2986173, allele_seq_length=6397), index=1)
    assert feat.type == "misc_feature"
    assert feat.annotations.get("mge_putative_composite") is True


def test_entry_to_feature_no_marker_for_mobile_element():
    # cn evidence=1 (DBヒット) は mobile_element でマーカー無し
    _sid, feat = entry_to_feature(
        _full_entry(evidence=1, start=2292, end=8429, allele_seq_length=6138), index=1)
    assert feat.type == "mobile_element"
    assert "mge_putative_composite" not in feat.annotations


def _mini_genome_with_mge():
    from collections import OrderedDict
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import FeatureLocation
    from dfc.models.bio_feature import ExtendedFeature

    rec = SeqRecord(Seq("A" * 3000), id="seq1", name="seq1")
    rec.annotations = {}
    src = ExtendedFeature(location=FeatureLocation(0, 3000, strand=1),
                          type="source", id="seq1", seq_id="seq1")
    src.qualifiers = {"mol_type": ["genomic DNA"], "organism": [""]}
    me = ExtendedFeature(location=FeatureLocation(100, 900, strand=1),
                         type="mobile_element", id="MGE_1", seq_id="seq1")
    me.qualifiers = {"mobile_element_type": ["insertion sequence:ISxx"],
                     "note": ["MobileElementFinder: ISxx"]}
    pc = ExtendedFeature(location=FeatureLocation(100, 2000, strand=1),
                         type="misc_feature", id="MGE_2", seq_id="seq1")
    pc.qualifiers = {"note": ["MobileElementFinder: cn_1900_ISxx",
                              "putative composite transposon (predicted by MobileElementFinder)"]}
    pc.annotations["mge_putative_composite"] = True
    rec.features = [src, me, pc]

    class FakeGenome:
        pass
    g = FakeGenome()
    g.seq_records = OrderedDict({"seq1": rec})
    g.features = OrderedDict({"seq1": src, "MGE_1": me, "MGE_2": pc})
    g.complete = False
    g.project_type = ""
    return g


def test_ann_excludes_putative_composite_by_default(tmp_path):
    from dfc.utils.ddbj_submission import create_ddbj_submission_file
    g = _mini_genome_with_mge()
    ann = str(tmp_path / "out.ann"); fa = str(tmp_path / "out.fasta")
    create_ddbj_submission_file(g, {}, ann, fa, verbosity=1)  # include_putative_composite 既定=False
    text = open(ann).read()
    assert "putative composite transposon" not in text   # putative composite は除外
    assert "insertion sequence:ISxx" in text             # mobile_element(IS) は残る


def test_ann_includes_putative_composite_when_enabled(tmp_path):
    from dfc.utils.ddbj_submission import create_ddbj_submission_file
    g = _mini_genome_with_mge()
    ann = str(tmp_path / "out.ann"); fa = str(tmp_path / "out.fasta")
    create_ddbj_submission_file(g, {}, ann, fa, verbosity=1, include_putative_composite=True)
    text = open(ann).read()
    assert "putative composite transposon" in text


import json
import os
from dfc.tools.dfast_mge import parse_mge_results

_FIXTURE = os.path.join(os.path.dirname(__file__), "data", "mge_pOXA48.json")


def test_parse_mge_results_fixture():
    data = json.load(open(_FIXTURE))
    result = parse_mge_results(data)
    # seq_id 別 dict
    assert "JN626286.1" in result
    feats = result["JN626286.1"]
    assert len(feats) == 1
    feat = feats[0]
    assert feat.type == "mobile_element"
    assert feat.qualifiers["mobile_element_type"] == ["transposon:Tn1999"]
    assert int(feat.location.start) == 2291
    assert int(feat.location.end) == 8429
    assert feat.location.strand == -1


def test_parse_mge_results_empty():
    assert parse_mge_results({"result": []}) == {}


from dfc.tools.dfast_mge import MobileElementFinder


def test_getcommand_includes_cmd_options(tmp_path):
    tool = MobileElementFinder(options={"cmd_options": "--min-coverage 0.95"},
                               workDir=str(tmp_path))
    cmd = tool.getCommand()
    assert cmd[0] == "mefinder"
    assert "find" in cmd
    assert "-c" in cmd
    assert "--json" in cmd
    # cmd_options が分割されて含まれる
    assert "--min-coverage" in cmd and "0.95" in cmd
    # temp-dir が contig_annotation/mge_finder 配下
    ti = cmd.index("--temp-dir")
    assert cmd[ti + 1].endswith(os.path.join("contig_annotation", "mge_finder"))


def test_getcommand_includes_db_path(tmp_path):
    tool = MobileElementFinder(options={"db_path": "/db/mefinder_db"},
                               workDir=str(tmp_path))
    cmd = tool.getCommand()
    assert "--db-path" in cmd
    assert cmd[cmd.index("--db-path") + 1] == "/db/mefinder_db"


def test_getcommand_no_db_path_when_unset(tmp_path):
    tool = MobileElementFinder(options={}, workDir=str(tmp_path))
    cmd = tool.getCommand()
    assert "--db-path" not in cmd


def test_getfeatures_reads_json(tmp_path):
    # 出力先に フィクスチャ JSON を置いて getFeatures が読めること
    tool = MobileElementFinder(options={}, workDir=str(tmp_path))
    os.makedirs(tool.output_directory, exist_ok=True)
    import shutil
    shutil.copy(_FIXTURE, tool.result_file)
    feats = tool.getFeatures()
    assert "JN626286.1" in feats
    assert feats["JN626286.1"][0].type == "mobile_element"


def test_getresult_source_notes_empty_report_populated(tmp_path):
    # MEF は座標付き feature を getFeatures() で返すため source_notes は空。
    # report には amr_summary.tsv の ## 行用の MGE 要約を入れる（PlasmidFinder 同形式）。
    tool = MobileElementFinder(options={}, workDir=str(tmp_path))
    os.makedirs(tool.output_directory, exist_ok=True)
    import shutil
    shutil.copy(_FIXTURE, tool.result_file)
    source_notes, report = tool.getResult()
    assert source_notes == {}
    assert "JN626286.1" in report
    assert any("MobileElementFinder: Tn1999" in line for line in report["JN626286.1"])


def test_base_getfeatures_default(tmp_path):
    from dfc.tools.base_tools import ContigAnnotationTool

    class Dummy(ContigAnnotationTool):
        NAME = "Dummy"
        def getCommand(self):
            return ["true"]
        def getResult(self):
            return {}, {}

    tool = Dummy(options={}, workDir=str(tmp_path))
    assert tool.getFeatures() == {}


def test_add_contig_features_appends_and_sorts():
    from collections import OrderedDict
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from dfc.models.bio_feature import ExtendedFeature
    from Bio.SeqFeature import FeatureLocation

    # genome.add_contig_features を単体で叩くため、最小の擬似 genome を作る
    import dfc.genome as genome_mod

    class FakeGenome:
        def __init__(self):
            rec = SeqRecord(Seq("A" * 10000), id="JN626286.1")
            self.seq_records = OrderedDict({"JN626286.1": rec})
            self.features = OrderedDict()
        # 実装をバインド
        add_contig_features = genome_mod.Genome.add_contig_features
        sort_features = genome_mod.Genome.sort_features
        set_feature_dictionary = genome_mod.Genome.set_feature_dictionary

    g = FakeGenome()
    feat = ExtendedFeature(location=FeatureLocation(2291, 8429, strand=-1),
                           type="mobile_element", id="MGE_1", seq_id="JN626286.1")
    feat.qualifiers = {"mobile_element_type": ["transposon:Tn1999"]}
    g.add_contig_features({"JN626286.1": [feat]})

    rec = g.seq_records["JN626286.1"]
    assert any(f.type == "mobile_element" for f in rec.features)
    assert "MGE_1" in g.features
