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
