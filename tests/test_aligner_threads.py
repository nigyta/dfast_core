# coding: UTF8
"""DIAMOND and GHOSTX run as one process with CPU threads; single-threaded aligners get one query file per CPU."""
from logging import getLogger

from Bio.SeqFeature import FeatureLocation

from dfc.components.baseComponent import BaseAnnotationComponent
from dfc.models.bio_feature import ExtendedFeature
from dfc.tools.diamond import Diamond
from dfc.tools.blastp import Blastp
from dfc.tools.ghostx import Ghostx


def _component(tmp_path, aligner):
    features = {}
    for i in range(10):
        f = ExtendedFeature(location=FeatureLocation(0, 30, 1), type="CDS", id="MGA_{}".format(i),
                            qualifiers={"translation": ["MKKLLPTAAAGLLLLAAQPAMA"]})
        features[f.id] = f
    component = BaseAnnotationComponent.__new__(BaseAnnotationComponent)  # skip directory and option setup
    component.logger, component.workDir, component.CPU = getLogger(__name__), str(tmp_path), 4
    component.skipAnnotatedFeatures, component.query_files, component.query_sequences = False, {}, {}
    component.genome = type("FakeGenome", (), {"features": features})()
    component.aligner = aligner
    return component


def test_diamond_one_process_with_cpu_threads(tmp_path, monkeypatch):
    monkeypatch.setattr(Diamond, "version", "2.2.8")  # skip the external version check
    component = _component(tmp_path, Diamond())
    component.prepareAlignerQueries()
    assert len(component.query_files) == 1
    cmd = " ".join(component.aligner.get_command("q.faa", "db", "out.tsv"))
    assert "--very-sensitive --threads 4 " in cmd


def test_ghostx_one_process_with_cpu_threads(tmp_path, monkeypatch):
    monkeypatch.setattr(Ghostx, "version", "1.3.6")
    component = _component(tmp_path, Ghostx())
    component.prepareAlignerQueries()
    assert len(component.query_files) == 1
    assert " ".join(component.aligner.get_command("q.faa", "db", "out.tsv")).endswith("-b 1 -v 1 -a 4")


def test_blastp_one_single_threaded_process_per_cpu(tmp_path, monkeypatch):
    monkeypatch.setattr(Blastp, "version", "2.17.0")
    component = _component(tmp_path, Blastp())
    component.prepareAlignerQueries()
    assert len(component.query_files) == 4 and component.aligner.threads == 1
