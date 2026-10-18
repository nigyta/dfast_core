import os
from dfc.utils.config_util import load_config, enable_amr, enable_amr_map, set_genetic_code
APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _enabled(config, name):
    return [s["enabled"] for s in config.FUNCTIONAL_ANNOTATION if s["component_name"] == name]

def _config(path="dfc/default_config.py"):
    return load_config(APP_ROOT, os.path.join(APP_ROOT, path))

def test_entries_disabled_by_default():
    for path in ["dfc/default_config.py", "example/dfc_config/dfast_web_config.py"]:
        assert _enabled(_config(path), "NuclMapping") == [False, False]

def test_amr_alone_leaves_mapping_disabled():
    config = _config(); enable_amr(config)
    assert _enabled(config, "NuclMapping") == [False, False]

def test_amr_map_enables_mapping_and_plasmidfinder():
    config = _config(); enable_amr_map(config)
    assert _enabled(config, "NuclMapping") == [True, True]
    assert all(not e for e in _enabled(config, "NuclSearch"))
    assert [t["enabled"] for t in config.CONTIG_ANNOTATION if t["tool_name"] == "Plasmidfinder"] == [True]

def test_amr_map_with_amr_disables_nuclsearch(caplog):
    config = _config(); enable_amr(config); enable_amr_map(config)
    assert all(not e for e in _enabled(config, "NuclSearch")) and "NuclSearch" in caplog.text

def test_genetic_code_propagates():
    config = _config(); set_genetic_code(config, 4)
    assert [s["options"]["transl_table"] for s in config.FUNCTIONAL_ANNOTATION if s["component_name"] == "NuclMapping"] == [4, 4]
