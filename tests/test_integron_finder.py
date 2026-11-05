# coding: UTF8
"""IntegronFinder result parsing (dfc.tools.integron_finder). The fixture is cut from a real run on A. baumannii AYE."""
import os

import pytest

from dfc.tools.integron_finder import read_integrons_file, group_integrons, integron_to_features, integron_summary

DATA = os.path.join(os.path.dirname(__file__), "data", "integron_finder_AYE.integrons")


def _integrons():
    return {i["id"]: i for i in group_integrons(read_integrons_file(DATA))}


def test_group_integrons():
    integrons = _integrons()
    assert set(integrons) == {"integron_01", "integron_02", "integron_03"}
    complete = integrons["integron_01"]
    assert (complete["seq_id"], complete["type"], complete["start"], complete["end"]) == \
        ("sequence1", "complete", 3620981, 3625513)
    assert len(complete["attC"]) == 4 and complete["intI"]
    assert integrons["integron_02"]["type"] == "CALIN" and not integrons["integron_02"]["intI"]
    assert integrons["integron_03"]["type"] == "In0" and integrons["integron_03"]["attC"] == []


def test_integron_to_features():
    features = integron_to_features(_integrons()["integron_01"])
    integron, attcs = features[0], features[1:]
    assert integron.type == "mobile_element"
    assert integron.qualifiers["mobile_element_type"] == ["integron"]
    assert integron.qualifiers["note"] == ["IntegronFinder: complete integron; attC sites: 4; integrase: intI"]
    assert (int(integron.location.start), int(integron.location.end)) == (3620980, 3625513)
    assert [f.type for f in attcs] == ["misc_recomb"] * 4
    assert (int(attcs[0].location.start), int(attcs[0].location.end), attcs[0].location.strand) == (3621430, 3621490, -1)
    assert len({f.id for f in features}) == 5


def test_integron_summary_labels():
    integrons = _integrons()
    assert integron_summary(integrons["integron_02"]).startswith("CALIN (cluster of attC sites lacking integron integrase)")
    assert integron_summary(integrons["integron_03"]) == \
        "In0 (integron integrase without attC sites); attC sites: 0; integrase: intI"


def test_no_integron(tmp_path):
    empty = tmp_path / "genome.integrons"
    empty.write_text("# integron_finder 2.0.6 \n# cmd: integron_finder genome.fna\n# No Integron found\n")
    assert group_integrons(read_integrons_file(str(empty))) == []


def test_missing_executable_exits_with_hint(tmp_path, monkeypatch, caplog):
    from dfc.tools.integron_finder import IntegronFinder
    monkeypatch.setenv("PATH", str(tmp_path))  # no integron_finder
    monkeypatch.setattr(IntegronFinder, "version", None)
    with pytest.raises(SystemExit):
        IntegronFinder(workDir=str(tmp_path))
    assert 'conda install -c bioconda "integron_finder>=2.0.6"' in caplog.text
