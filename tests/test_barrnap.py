# coding: UTF8
"""rRNA prediction with nhmmer and Barrnap's rules, rRNA model selection, and building models from Rfam."""
import gzip
import shutil

import pytest

from dfc.tools.barrnap import Barrnap, parse_nhmmer_tblout
from dfc.utils.config_util import set_rrna_model
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
    hits = parse_nhmmer_tblout(str(tbl))
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
    assert [h[2] - h[1] + 1 for h in parse_nhmmer_tblout(str(tbl), reject=0.6)] == [1491, 95, 2948]


@pytest.fixture
def rrna_db(tmp_path, monkeypatch):
    monkeypatch.setattr(Barrnap, "version", "3.4")  # skip the external version check
    (tmp_path / "barrnap_0.9_bac.hmm").write_text("")
    (tmp_path / "Rfam_15.1_bac.hmm").write_text("")
    return str(tmp_path)


def test_model_selection_and_options(rrna_db):
    tool = Barrnap(options={"db_dir": rrna_db, "cmd_options": "--reject 0.25 --evalue 1e-5"})
    assert tool.hmm_file.endswith("barrnap_0.9_bac.hmm") and tool.model_label == "BarrnapHmm:0.9"
    assert (tool.params.reject, tool.params.lencutoff, tool.params.evalue) == (0.25, 0.8, 1e-5)
    cmd = tool.getCommand()
    assert cmd[:5] == ["nhmmer", "--cpu", "1", "-E", "1e-05"] and cmd[-2] == tool.hmm_file
    tool = Barrnap(options={"db_dir": rrna_db, "model": "rfam"})
    assert tool.hmm_file.endswith("Rfam_15.1_bac.hmm") and tool.model_label == "Rfam:15.1"
    assert tool.params.reject == 0.5


def test_missing_model_aborts(tmp_path, monkeypatch):
    monkeypatch.setattr(Barrnap, "version", "3.4")
    with pytest.raises(SystemExit):
        Barrnap(options={"db_dir": str(tmp_path)})


def test_set_rrna_model():
    class Config:
        STRUCTURAL_ANNOTATION = [{"tool_name": "Barrnap", "options": {"model": "barrnap"}}, {"tool_name": "Aragorn"}]
    set_rrna_model(Config, "rfam")
    assert Config.STRUCTURAL_ANNOTATION[0]["options"]["model"] == "rfam" and "options" not in Config.STRUCTURAL_ANNOTATION[1]


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
