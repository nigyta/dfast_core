# coding: UTF8
"""LAST for PseudoGeneDetection: genetic code via lastal -G, minimum version, and MAF score lines."""
import pytest

from dfc.components.PseudoGeneDetection import PseudoGeneDetection
from dfc.tools.last import Lastal


@pytest.fixture
def lastal_1654(monkeypatch):
    monkeypatch.setattr(Lastal, "version", "1654")  # skip the external version check


@pytest.mark.parametrize("transl_table", [11, 4, 25])
def test_lastal_passes_genetic_code(lastal_1654, transl_table):
    cmd = Lastal(options={"transl_table": transl_table}).get_command("query.fna", "reference", "out.maf")
    assert cmd[:3] == ["lastal", "-G", str(transl_table)]
    assert "lastal4" not in cmd and "-F 15" in cmd


def test_old_lastal_is_rejected(monkeypatch):
    monkeypatch.setattr(Lastal, "version", "959")  # bundled until DFAST 1.4.3; no NCBI genetic codes for -G
    with pytest.raises(SystemExit):
        Lastal(options={})


@pytest.mark.parametrize("line", ["a score=479 EG2=1.8e-59 E=3e-70",   # LAST < 1648
                                  "a score=479 EG=1.5e-63 E=2.8e-70"])  # LAST >= 1648
def test_maf_score_line(line):
    m = PseudoGeneDetection.PAT_SCORE.search(line)
    assert m and m.group(1) == "479" and float(m.group(3)) < 1e-60
