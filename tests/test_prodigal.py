# coding: UTF8
"""Prodigal wrapper: partial-gene parsing, broken-build detection and meta mode for short inputs.
The fixtures are cut from real Prodigal 2.6.3 output (bioconda for prodigal_edges.gff, the Ubuntu apt build for
prodigal_broken.gff); the partial=11 line is hand-written."""
import os

import pytest
from Bio.SeqFeature import BeforePosition, AfterPosition

from dfc.tools.prodigal import Prodigal

DATA = os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture(autouse=True)
def _no_version_check(monkeypatch):
    monkeypatch.setattr(Prodigal, "version", "2.6.3")


def _write_fasta(path, lengths):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for name, n in lengths.items():
            f.write(">{0}\n{1}\n".format(name, "A" * n))


def _tool(tmp_path, lengths, transl_table=11, gff=None):
    _write_fasta(os.path.join(str(tmp_path), "input", "genome.fna"), lengths)  # Genome writes it before tools are created
    tool = Prodigal({"transl_table": transl_table}, workDir=str(tmp_path))
    if gff:
        tool.outputFile = os.path.join(DATA, gff)
    return tool


def _features(tmp_path):
    return _tool(tmp_path, {"edge01": 30000, "edge11": 30000, "both": 300}, gff="prodigal_edges.gff").getFeatures()


def test_parse_edge_gene_plus(tmp_path):
    f = _features(tmp_path)["edge01"][0]
    assert str(f.location) == "[<0:2463](+)"
    assert f.annotations["partial_flag"] == "10"
    assert f.qualifiers["codon_start"] == [1]


def test_parse_edge_gene_minus(tmp_path):
    f = _features(tmp_path)["edge11"][0]
    assert str(f.location) == "[29214:>30000](-)"
    assert f.annotations["partial_flag"] == "01"


def test_parse_both_ends_partial(tmp_path):
    f = _features(tmp_path)["both"][0]
    assert f.location.start == BeforePosition(0) and isinstance(f.location.start, BeforePosition)
    assert f.location.end == AfterPosition(300) and isinstance(f.location.end, AfterPosition)
    assert f.qualifiers["codon_start"] == [2]


def test_broken_gff_exits(tmp_path, caplog):
    tool = _tool(tmp_path, {"a": 30000}, gff="prodigal_broken.gff")
    with pytest.raises(SystemExit):
        tool.getFeatures()
    assert "partial" in caplog.text


def test_short_input_uses_meta(tmp_path):
    assert "-p meta" in " ".join(_tool(tmp_path, {"a": 4000, "b": 6000}).getCommand())


def test_short_input_non11_exits(tmp_path, caplog):
    with pytest.raises(SystemExit):
        _tool(tmp_path, {"a": 10000}, transl_table=25)
    assert "20,000" in caplog.text and "genetic code" in caplog.text


def test_long_input_single_mode(tmp_path):
    cmd = " ".join(_tool(tmp_path, {"a": 25000}, transl_table=25).getCommand())
    assert "-p meta" not in cmd and "-g 25" in cmd
