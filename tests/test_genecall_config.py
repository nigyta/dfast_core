import os
import pytest
from dfc.utils.config_util import (load_config, enable_mga, enable_metagenome_mode, set_genetic_code,
                                   check_genetic_code)
from dfc.utils import fix_origin
from tests.test_overlap_rules import _run, _cds, _gap

APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _config(path="dfc/default_config.py"):
    return load_config(APP_ROOT, os.path.join(APP_ROOT, path))


def _cds_tools(config):
    # example/dfast_web_config.py has no "target" keys, so select by tool name
    return {s["tool_name"]: s for s in config.STRUCTURAL_ANNOTATION if s["tool_name"] in ("MGA", "Prodigal")}


@pytest.mark.parametrize("path", ["dfc/default_config.py", "example/dfc_config/dfast_web_config.py",
                                  "example/dfast_web_config.py", "example/test_config.py"])
def test_default_is_prodigal(path):
    tools = _cds_tools(_config(path))
    assert tools["Prodigal"]["enabled"] and not tools["MGA"]["enabled"]


def test_use_mga():
    config = _config()
    enable_mga(config)
    tools = _cds_tools(config)
    assert tools["MGA"]["enabled"] and not tools["Prodigal"]["enabled"]


def test_gcode25_sets_all_tables(caplog):
    config = _config()
    set_genetic_code(config, 25)
    structural = {s["tool_name"]: s for s in config.STRUCTURAL_ANNOTATION if "transl_table" in s["options"]}
    assert structural["Prodigal"]["options"]["transl_table"] == 25
    assert structural["Aragorn"]["options"]["transl_table"] == 25
    for name in ("PseudoGeneDetection", "NuclMapping"):
        assert all(s["options"]["transl_table"] == 25 for s in config.FUNCTIONAL_ANNOTATION
                   if s["component_name"] == name)
    assert not [r for r in caplog.records if r.levelname == "WARNING"]


def test_unsupported_gcode_warns(caplog):
    set_genetic_code(_config(), 2)
    assert any("Unsupported genetic code" in r.message for r in caplog.records)


def test_mga_with_gcode25_exits():
    config = _config()
    enable_mga(config)
    set_genetic_code(config, 25)
    with pytest.raises(SystemExit):
        check_genetic_code(config)


def test_metagenome_with_gcode4_exits():
    config = _config()
    set_genetic_code(config, 4)
    enable_metagenome_mode(config)
    with pytest.raises(SystemExit):
        check_genetic_code(config)


def test_genetic_code_guard_passes_for_valid_combinations():
    check_genetic_code(_config())
    config = _config()
    set_genetic_code(config, 4)
    check_genetic_code(config)
    config = _config()
    enable_mga(config)
    enable_metagenome_mode(config)
    check_genetic_code(config)


def test_prodigal_cds_across_gap_removed():
    # a Prodigal CDS spanning an assembly_gap is handled exactly like an MGA one
    results = {}
    for name in ("MGA_1", "Prodigal_1"):
        ids, _ = _run([_cds(100, 1000, name), _gap(400, 600)])
        results[name] = [i for i in ids if i != "GAP_1"]
    assert results["MGA_1"] == [] and results["Prodigal_1"] == []


def test_fix_origin_config_uses_prodigal():
    assert fix_origin.Config.STRUCTURAL_ANNOTATION[0]["tool_name"] == "Prodigal"
