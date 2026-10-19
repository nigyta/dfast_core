import logging
import os
from dfc.utils.config_util import load_config, enable_amr, enable_amr_map, set_genetic_code
APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _enabled(config, name, db_names=("CARD", "VFDB")):
    return [s["enabled"] for s in config.FUNCTIONAL_ANNOTATION if s["component_name"] == name
            and s["options"].get("db_name") in db_names]

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
    assert all(not e for e in _enabled(config, "NuclSearch"))
    assert any(r.levelno == logging.WARNING and "NuclSearch" in r.getMessage() for r in caplog.records)

def test_amr_map_alone_does_not_warn(caplog):
    config = _config(); enable_amr_map(config)
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING and "NuclSearch" in r.getMessage()]

def test_genetic_code_propagates():
    config = _config(); set_genetic_code(config, 4)
    assert [s["options"]["transl_table"] for s in config.FUNCTIONAL_ANNOTATION if s["component_name"] == "NuclMapping"] == [4, 4, 4]

def test_mapping_option_defaults():
    for path in ["dfc/default_config.py", "example/dfc_config/dfast_web_config.py"]:
        options = [s["options"] for s in _config(path).FUNCTIONAL_ANNOTATION if s["component_name"] == "NuclMapping"]
        assert [(o["blast_task"], o["pseudogene_pident_cutoff"], o["overlap_fraction"], o["truncated_fraction"])
                for o in options] == [("dc-megablast", 95, 0.2, 0.7)] * 3

def test_amr_map_keeps_other_nuclsearch_entries(caplog):
    config = _config()
    custom = {"component_name": "NuclSearch", "enabled": True,   # a user's own NuclSearch entry
              "options": {"database": "/my/db/MY_DB", "db_name": "MY_DB"}}
    config.FUNCTIONAL_ANNOTATION.append(custom)
    enable_amr(config); enable_amr_map(config)
    assert custom["enabled"]
    assert _enabled(config, "NuclSearch") == [False, False]
    assert "MY_DB" not in caplog.text


def test_custom_nucleotide_reference_uses_nuclmapping():
    for path in ["dfc/default_config.py", "example/dfc_config/dfast_web_config.py"]:
        config = _config(path)
        assert _enabled(config, "NuclSearch", ("PLASMID_DB",)) == []
        assert _enabled(config, "NuclMapping", ("PLASMID_DB",)) == [False]
        enable_amr_map(config)                  # --amr_map maps CARD/VFDB only
        assert _enabled(config, "NuclMapping", ("PLASMID_DB",)) == [False]
