import os
import pytest
from dfc.utils.config_util import (load_config, enable_mga, enable_prodigal, enable_metagenome_mode, set_genetic_code,
                                   check_genetic_code, disable_cds_prediction)
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


def test_use_prodigal_reselects_prodigal_in_mga_config():
    # a 1.5-era custom config with MGA enabled: --use_prodigal (enable_prodigal) must still switch to Prodigal
    config = _config()
    enable_mga(config)
    enable_prodigal(config)
    tools = _cds_tools(config)
    assert tools["Prodigal"]["enabled"] and not tools["MGA"]["enabled"]


def test_config_without_target_keys():
    config = _config("example/test_config.py")
    enable_mga(config)
    tools = _cds_tools(config)
    assert tools["MGA"]["enabled"] and not tools["Prodigal"]["enabled"]
    set_genetic_code(config, 25)
    with pytest.raises(SystemExit):
        check_genetic_code(config)


@pytest.mark.parametrize("path", ["dfc/default_config.py", "example/dfc_config/dfast_web_config.py",
                                  "example/dfast_web_config.py", "example/test_config.py"])
def test_no_rnammer(path):
    from dfc.structuralAnnotation import TOOLS
    assert "RNAmmer" not in TOOLS
    assert all(s.get("tool_name") != "RNAmmer" for s in _config(path).STRUCTURAL_ANNOTATION)


def test_dfast_help_has_no_use_rnammer():
    import subprocess
    import sys
    out = subprocess.run([sys.executable, os.path.join(APP_ROOT, "dfast"), "-h"], capture_output=True, text=True).stdout
    assert "--use_rnammer" not in out
    assert "--use_mga" in out


def test_metagenome_with_target_less_config():
    config = _config("example/test_config.py")
    enable_metagenome_mode(config)
    assert _cds_tools(config)["Prodigal"]["options"]["cmd_options"] == "-p meta"


def test_disable_cds_with_target_less_config():
    config = _config("example/test_config.py")
    disable_cds_prediction(config)
    assert not any(s["enabled"] for s in _cds_tools(config).values())


def test_fix_origin_passes_genetic_code(monkeypatch, tmp_path):
    seen = {}

    class Fake:
        def __init__(self, *args):
            seen["config"] = args[-1]  # Genome(config) / StructuralAnnotation(genome, config)
            self.seq_records = {}

        def execute(self):
            pass

    for name in ("Genome", "StructuralAnnotation", "FunctionalAnnotation"):
        monkeypatch.setattr(fix_origin, name, Fake)
    fix_origin.fix_origin("x.fna", output_file=str(tmp_path / "out.fna"), transl_table=4)
    assert seen["config"].STRUCTURAL_ANNOTATION[0]["options"]["transl_table"] == 4
