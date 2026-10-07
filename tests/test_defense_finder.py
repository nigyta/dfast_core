# coding: UTF8
"""DefenseFinder component (dfc.components.DefenseFinder) without running defense-finder.
The genes table is cut from a real run on pOXA-48, with the gembase ids that the component generates."""
import json
import logging
import os
import shutil
from collections import OrderedDict

from Bio.SeqFeature import FeatureLocation
from dfc.components.DefenseFinder import DefenseFinder
from dfc.models.bio_feature import ExtendedFeature
from dfc.models.hit import DefenseHit
from dfc.utils.config_util import load_config, enable_amr, enable_mge

APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(os.path.dirname(__file__), "data", "defense_finder_genes.tsv")


class _Record:
    def __init__(self, features):
        self.features = features


class _Genome:
    def __init__(self, records):
        self.seq_records = OrderedDict(records)
        self.features = {f.id: f for rec in self.seq_records.values() for f in rec.features}


def _cds(fid, start, seq_id):
    return ExtendedFeature(location=FeatureLocation(start, start + 300, strand=1), type="CDS", id=fid,
                           seq_id=seq_id, qualifiers={"translation": ["MAAA"]})


def _component(tmp_path, genome):
    # bypass __init__, which checks the defense-finder executable and the models directory
    component = DefenseFinder.__new__(DefenseFinder)
    component.genome, component.workDir = genome, str(tmp_path)
    component.query_file = str(tmp_path / "query.faa")
    component.query_ids = {}
    component.logger = logging.getLogger(__name__)
    return component


def test_prepare_query_gembase_ids(tmp_path):
    genome = _Genome([
        ("seq1", _Record([_cds("MGA_1", 0, "seq1"), ExtendedFeature(type="tRNA", id="t1"), _cds("MGA_2", 500, "seq1")])),
        ("seq2", _Record([_cds("MGA_3", 0, "seq2")])),
    ])
    component = _component(tmp_path, genome)
    component.prepare_query()
    assert component.query_ids == {"R0001_00001": "MGA_1", "R0001_00002": "MGA_2", "R0002_00001": "MGA_3"}
    assert open(component.query_file).read().startswith(">R0001_00001\nMAAA\n")


def test_set_results(tmp_path):
    genome = _Genome([("seq1", _Record([_cds("MGA_10", 0, "seq1"), _cds("MGA_11", 400, "seq1"), _cds("MGA_47", 900, "seq1")]))])
    component = _component(tmp_path, genome)
    component.prepare_query()
    shutil.copy(DATA, str(tmp_path / "query_defense_finder_genes.tsv"))
    component.set_results()

    hit = genome.features["MGA_10"].secondary_hits[0]
    assert str(hit) == "DefenseFinder: MazE of MazEF defense system"
    assert str(genome.features["MGA_47"].secondary_hits[0]) == \
        "DefenseFinder: klca of ardb_klca_acric11 (Anti_RM) anti-defense system"
    systems = json.load(open(str(tmp_path / "systems.json")))
    assert sorted(len(s["features"]) for s in systems.values()) == [1, 2]


def test_defense_hit_note_at_any_verbosity():
    feature = ExtendedFeature(type="CDS", id="MGA_1", qualifiers={})
    DefenseHit("RM_Type_I", "N6_RM_Type_I", "Defense", "RM__HsdR", "sys1").assign_as_note(feature, verbosity=1)
    assert feature.qualifiers["note"] == ["DefenseFinder: HsdR of N6_RM_Type_I (RM_Type_I) defense system"]


def test_enable_amr_and_mge():
    config = load_config(APP_ROOT, os.path.join(APP_ROOT, "dfc", "default_config.py"))
    enabled = lambda key, name, settings: [s["enabled"] for s in settings if s.get(key) == name]
    assert enabled("component_name", "DefenseFinder", config.FUNCTIONAL_ANNOTATION) == [False]
    assert enabled("tool_name", "IntegronFinder", config.CONTIG_ANNOTATION) == [False]
    enable_amr(config)
    enable_mge(config)
    assert enabled("component_name", "DefenseFinder", config.FUNCTIONAL_ANNOTATION) == [True]
    assert enabled("tool_name", "IntegronFinder", config.CONTIG_ANNOTATION) == [True]
