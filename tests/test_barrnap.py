# coding: UTF8
"""rRNA prediction with nhmmer and Barrnap's rules, rRNA model selection, and building models from Rfam."""
import gzip
import shutil

import pytest

from dfc.tools.barrnap import Barrnap, apply_barrnap_rules, compete, read_nhmmer_tblout
from dfc.utils.config_util import set_rrna_kingdoms, set_rrna_model
from dfc.utils.reffile_util import build_rrna_hmm

# nhmmer --tblout: target, acc, query, acc, hmmfrom, hmm to, alifrom, alito, envfrom, envto, sq len, strand, E-value, ...
TBLOUT = """\
# target name  accession  query name  accession  hmmfrom  hmm to  alifrom  alito ...
contig2 - 23S_rRNA RF02541 1 2906 3049 102 3049 102 3100 - 0 1500.0 10.0 -
contig1 - 16S_rRNA RF00177 1 1533 2000 2800 2000 2800 5000 + 1.6e-68 300.0 1.0 -
contig1 - 5S_rRNA RF00001 1 119 101 195 101 195 5000 + 2e-09 50.0 0.1 -
contig1 - 5S_rRNA RF00001 1 119 300 358 300 358 5000 + 1e-07 30.0 0.1 -
contig1 - 16S_rRNA RF00177 1 1533 10 1500 10 1500 5000 + 0 1400.0 1.0 -
"""


def test_parse_applies_barrnap_rules(tmp_path):
    tbl = tmp_path / "out.txt"
    tbl.write_text(TBLOUT)
    hits = apply_barrnap_rules(read_nhmmer_tblout(str(tbl), "bac"))
    # sorted by sequence and left end; the 59 bp 5S is kept (not below int(0.5 * 119) = 59) as a partial hit
    assert [(h[0], h[1], h[2], h[3], h[4]) for h in hits] == [
        ("contig1", 10, 1500, "+", "16S_rRNA"), ("contig1", 101, 195, "+", "5S_rRNA"),
        ("contig1", 300, 358, "+", "5S_rRNA"), ("contig1", 2000, 2800, "+", "16S_rRNA"),
        ("contig2", 102, 3049, "-", "23S_rRNA")]
    full16s, s5_95bp, s5_59bp, partial16s, s23 = hits
    assert full16s[5] == "16S ribosomal RNA" and full16s[6] is None
    assert s5_95bp[6] is None  # 95 bp is not below int(0.8 * 119) = 95, as in Barrnap
    assert s5_59bp[6] == "aligned only 49 percent of the 5S ribosomal RNA"
    assert partial16s[6] == "aligned only 50 percent of the 16S ribosomal RNA" and s23[7] == "0"
    # a higher reject threshold drops the short hits
    assert [h[2] - h[1] + 1 for h in apply_barrnap_rules(read_nhmmer_tblout(str(tbl), "bac"), reject=0.6)] == [1491, 95, 2948]


def _hit(kingdom, gene, left, right, score, strand="+"):
    return {"seq_id": "c1", "gene": gene, "left": left, "right": right, "strand": strand, "evalue": "0", "score": score, "kingdom": kingdom}


def test_compete_keeps_best_kingdom_per_subunit():
    hits = [_hit("bac", "16S_rRNA", 100, 1500, 800), _hit("arc", "16S_rRNA", 99, 1500, 1300),
            _hit("euk", "18S_rRNA", 90, 1600, 400),
            _hit("bac", "5S_rRNA", 3000, 3110, 80), _hit("euk", "5S_rRNA", 3000, 3110, 80),  # same model: tie goes to bac
            _hit("bac", "23S_rRNA", 1600, 4500, 2000), _hit("arc", "23S_rRNA", 1600, 4500, 2500, strand="-"),  # other strand
            _hit("bac", "16S_rRNA", 5000, 5600, 300), _hit("bac", "16S_rRNA", 5500, 6000, 200)]  # same kingdom: both kept
    kept = sorted((h["kingdom"], h["gene"], h["left"], h["strand"]) for h in compete(hits))
    assert kept == [("arc", "16S_rRNA", 99, "+"), ("arc", "23S_rRNA", 1600, "-"), ("bac", "16S_rRNA", 5000, "+"),
                    ("bac", "16S_rRNA", 5500, "+"), ("bac", "23S_rRNA", 1600, "+"), ("bac", "5S_rRNA", 3000, "+")]


@pytest.fixture
def rrna_db(tmp_path, monkeypatch):
    monkeypatch.setattr(Barrnap, "version", "3.4")  # skip the external version check
    for kingdom in ("bac", "arc", "euk"):
        (tmp_path / "barrnap_0.9_{}.hmm".format(kingdom)).write_text("")
        (tmp_path / "Rfam_15.1_{}.hmm".format(kingdom)).write_text("")
    return str(tmp_path)


def test_model_selection_and_options(rrna_db):
    tool = Barrnap(options={"db_dir": rrna_db, "cmd_options": "--reject 0.25 --evalue 1e-5"})
    assert tool.kingdoms == ["bac", "euk"]  # euk is always searched as a contamination check
    assert tool.hmm_files["bac"].endswith("barrnap_0.9_bac.hmm")
    assert tool.model_label == "BarrnapHmm:0.9"
    assert (tool.params.reject, tool.params.lencutoff, tool.params.evalue) == (0.25, 0.8, 1e-5)
    cmd = tool.getCommand()
    assert cmd[:5] == ["nhmmer", "--cpu", "1", "-E", "1e-05"] and cmd[-2] == tool.hmm_files["bac"]
    tool = Barrnap(options={"db_dir": rrna_db, "model": "rfam", "kingdoms": ["bac", "arc"]})
    assert tool.kingdoms == ["bac", "arc", "euk"]
    assert tool.hmm_files["euk"].endswith("Rfam_15.1_euk.hmm") and tool.model_label == "Rfam:15.1"
    assert tool.getCommand("arc")[-2].endswith("Rfam_15.1_arc.hmm") and tool.tblout("arc").endswith("Barrnap.arc.txt")
    assert tool.params.reject == 0.5
    assert Barrnap(options={"db_dir": rrna_db, "kingdoms": ["arc"], "check_eukaryotic": False}).kingdoms == ["arc"]


def test_missing_model_or_unknown_kingdom_aborts(tmp_path, rrna_db, monkeypatch):
    with pytest.raises(SystemExit):
        Barrnap(options={"db_dir": str(tmp_path / "empty")})
    with pytest.raises(SystemExit):
        Barrnap(options={"db_dir": rrna_db, "kingdoms": ["bac", "mito"]})
    with pytest.raises(SystemExit):
        Barrnap(options={"db_dir": rrna_db, "kingdoms": ["bac", "euk"]})  # euk is not a choice; it is always searched


def test_set_rrna_model():
    class Config:
        STRUCTURAL_ANNOTATION = [{"tool_name": "Barrnap", "options": {"model": "barrnap"}}, {"tool_name": "Aragorn"}]
    set_rrna_model(Config, "rfam")
    set_rrna_kingdoms(Config, ["bac", "arc"])
    assert Config.STRUCTURAL_ANNOTATION[0]["options"] == {"model": "rfam", "kingdoms": ["bac", "arc"]}
    assert "options" not in Config.STRUCTURAL_ANNOTATION[1]


SEED = """\
# STOCKHOLM 1.0
#=GF AC   RF00001
#=GF ID   5S_rRNA
s1  ACGUACGUACGUAGCUAGCUAGCUAGCUAGCAUCGAUCGA
s2  ACGUACGUACGUAGCUAGCUAGCUAGCUAGCAUCGAUCGA
s3  ACGUACGAACGUAGCUAGCUAGCUAGCUAGCAUCGAUCGA
//
# STOCKHOLM 1.0
#=GF AC   RF00005
#=GF ID   tRNA
s1  GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUC
s2  GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUC
//
"""


@pytest.mark.skipif(shutil.which("hmmbuild") is None, reason="hmmbuild (HMMER) is not on PATH")
def test_build_rrna_hmm(tmp_path):
    seed, out = tmp_path / "Rfam.seed.gz", tmp_path / "bac.hmm"
    with gzip.open(seed, "wt") as f:
        f.write(SEED)
    build_rrna_hmm(str(seed), str(out), [("5S_rRNA", "RF00001")])
    text = out.read_text()
    assert text.count("//") == 1 and "NAME  5S_rRNA" in text and "ALPH  RNA" in text
    with pytest.raises(ValueError):
        build_rrna_hmm(str(seed), str(out), [("16S_rRNA", "RF00177")])
