# coding: UTF8
"""Diced wrapper. tests/data/crispr_region.fna is NC_000913.3:2876201..2879963 (E. coli K-12). The CRISPR at
2877701..2878463 (13 repeats) was found identically by CRT and Diced with CRT parameters, so it is 0-based
[1500:2263] in the region."""
import os
from types import SimpleNamespace

import diced
import pytest
from Bio.SeqFeature import BeforePosition, AfterPosition

from dfc.structuralAnnotation import TOOLS
from dfc.tools.diced import Diced

DATA = os.path.join(os.path.dirname(__file__), "data", "crispr_region.fna")
REPEAT = "cggtttatccccgctggcgcggggaactc"
CRT_PARAMS = {"min_repeat_length": 19, "max_repeat_length": 38, "min_spacer_length": 19, "max_spacer_length": 48}


def _seq():
    return "".join(l.strip() for l in open(DATA) if not l.startswith(">"))


def _features(tmp_path, seq, name="region"):
    os.makedirs(os.path.join(str(tmp_path), "input"), exist_ok=True)
    with open(os.path.join(str(tmp_path), "input", "genome.fna"), "w") as f:
        f.write(">{0}\n{1}\n".format(name, seq))
    tool = Diced({}, workDir=str(tmp_path))
    tool.run()
    return tool.getFeatures()


def test_diced_finds_crispr(tmp_path):
    (f,) = _features(tmp_path, _seq())["region"]
    assert (int(f.location.start), int(f.location.end), f.location.strand) == (1500, 2263, 1)
    assert f.type == "CRISPR" and f.id == "Diced_1"
    assert f.qualifiers["rpt_unit_seq"] == [REPEAT]
    assert f.qualifiers["rpt_family"] == ["CRISPR"] and f.qualifiers["rpt_type"] == ["direct"]
    assert f.qualifiers["inference"] == ["COORDINATES:alignment:Diced:" + diced.__version__]


def test_diced_lowercase_input(tmp_path):
    (f,) = _features(tmp_path, _seq().lower())["region"]
    assert (int(f.location.start), int(f.location.end)) == (1500, 2263)
    assert f.qualifiers["rpt_unit_seq"] == [REPEAT]


def test_diced_crispr_at_contig_end(tmp_path):
    (f,) = _features(tmp_path, _seq()[:2263])["region"]
    assert f.location.end == AfterPosition(2263)
    assert not isinstance(f.location.start, BeforePosition)


def test_diced_crispr_at_contig_start(tmp_path, monkeypatch):
    # Diced never reports a CRISPR starting at 0 for this region (it drops the first repeat), so fake the hit.
    fake = SimpleNamespace(start=0, end=763, repeats=["ACGT" * 6] * 3)
    monkeypatch.setattr(diced, "scan", lambda seq, **kw: [fake])
    (f,) = _features(tmp_path, _seq()[1500:])["region"]
    assert f.location.start == BeforePosition(0) and f.location.end == 763


def test_diced_crt_params_passed(tmp_path, monkeypatch):
    seen = {}
    monkeypatch.setattr(diced, "scan", lambda seq, **kw: seen.update(kw) or [])
    assert _features(tmp_path, _seq()) == {}
    assert seen == CRT_PARAMS


def test_registered():
    assert "Diced" in TOOLS and "CRT" not in TOOLS
