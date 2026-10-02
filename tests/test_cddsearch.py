# coding: UTF8
"""CDDsearch with rpsbproc 0.5: version check and parsing of its output."""
import os
import stat
from logging import getLogger

import pytest
from Bio.SeqFeature import FeatureLocation

from dfc.components.CDDsearch import CDDsearch
from dfc.models.bio_feature import ExtendedFeature
from dfc.tools.rpsbproc import Rpsbproc

NEW_VERSION_OUTPUT = "rpsbproc: 0.5.1\n Package: blast 2.17.0, build Jan 22 2026 10:11:34\n"
OLD_VERSION_OUTPUT = ("Invalid switch 'v', ignored..\nInvalid switch 'e', ignored..\n"
                      "Post-RPSBLAST Processing Utility v0.11\n")


def _fake_rpsbproc(tmp_path, output, monkeypatch):
    exe = tmp_path / "rpsbproc"
    exe.write_text("#!/bin/sh\ncat <<'EOF'\n{}EOF\n".format(output))
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", "{}:{}".format(tmp_path, os.environ["PATH"]))
    monkeypatch.setattr(Rpsbproc, "version", None)


def test_rpsbproc_05_version_is_detected(tmp_path, monkeypatch):
    _fake_rpsbproc(tmp_path, NEW_VERSION_OUTPUT, monkeypatch)
    assert Rpsbproc(options={"rpsbproc_data": "x"}).version == "0.5.1"


def test_old_rpsbproc_is_rejected(tmp_path, monkeypatch):
    _fake_rpsbproc(tmp_path, OLD_VERSION_OUTPUT, monkeypatch)
    with pytest.raises(SystemExit):
        Rpsbproc(options={"rpsbproc_data": "x"})


# rpsbproc 0.5.1 output (-q): a blank line after the header, and a superfamily hit with PSSM-ID 0.
RPSBPROC_OUTPUT = """#Post-RPSBLAST Processing Utility

#Input data file:\talignment0.asn
#DATA
#SESSION\t<session-ordinal>\t<program>\t<database>\t<score-matrix>\t<evalue-threshold>
DATA
SESSION\t1\tblastp\t2.17.0+\tCog\tBLOSUM62\t1e-06
QUERY\tQuery_1\tPeptide\t698\tMGA_1
DOMAINS
1\tQuery_1\tSpecific\t441819\t2\t698\t0\t694.195\tCOG2217\tZntA\t-\t-
1\tQuery_1\tSuperfamily\t0\t29\t247\t1.80725e-35\t128.973\tcl46809\tMaf_flag10_N\tN\t-
ENDDOMAINS
ENDQUERY
ENDSESSION
ENDDATA
"""
CDDID = "441819\tCOG2217\tZntA\tCation-transporting P-type ATPase [Inorganic ion transport and metabolism]. \t717\n"


def test_parse_rpsbproc_05_output(tmp_path):
    (tmp_path / "rpsbproc0.out").write_text(RPSBPROC_OUTPUT)
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "cddid.tbl").write_text(CDDID)

    feature = ExtendedFeature(location=FeatureLocation(0, 2094, 1), type="CDS", id="MGA_1")
    search = CDDsearch.__new__(CDDsearch)  # skip tool version checks in __init__
    search.logger = getLogger(__name__)
    search.workDir = str(tmp_path)
    search.query_files = {0: "query0.fasta"}
    search.rpsbproc = type("FakeRpsbproc", (), {"rpsbproc_data": str(data_dir)})()
    search.genome = type("FakeGenome", (), {"features": {"MGA_1": feature}})()

    search.parse_result()

    assert len(feature.secondary_hits) == 1  # the PSSM-ID 0 superfamily line is skipped
    hit = feature.secondary_hits[0]
    assert (hit.accession, hit.result_type, hit.category, hit.description) == \
        ("COG2217", "COG", "P", "Cation-transporting P-type ATPase")
