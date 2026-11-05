import logging
import os
from dfc.utils.config_util import load_config, enable_amr, set_genetic_code
APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _enabled(config, name, db_names=("CARD", "VFDB")):
    return [s["enabled"] for s in config.FUNCTIONAL_ANNOTATION if s["component_name"] == name
            and s["options"].get("db_name") in db_names]

def _config(path="dfc/default_config.py"):
    return load_config(APP_ROOT, os.path.join(APP_ROOT, path))

def test_entries_disabled_by_default():
    for path in ["dfc/default_config.py", "example/dfc_config/dfast_web_config.py"]:
        assert _enabled(_config(path), "NuclMapping") == [False, False]

def test_amr_enables_mapping_and_plasmidfinder():
    for path in ["dfc/default_config.py", "example/dfc_config/dfast_web_config.py"]:
        config = _config(path); enable_amr(config)
        assert _enabled(config, "NuclMapping") == [True, True]
        assert _enabled(config, "NuclSearch") == [False, False]       # NuclSearch is used only when set in the config
        assert [t["enabled"] for t in config.CONTIG_ANNOTATION if t["tool_name"] == "Plasmidfinder"] == [True]


def test_nuclsearch_enabled_in_config_replaces_mapping_of_that_db(caplog):
    caplog.set_level(logging.INFO)
    config = _config()
    card = next(s for s in config.FUNCTIONAL_ANNOTATION
                if s["component_name"] == "NuclSearch" and s["options"]["db_name"] == "CARD")
    card["enabled"] = True                      # a user config: CARD with NuclSearch
    enable_amr(config)
    assert _enabled(config, "NuclSearch") == [True, False]
    assert _enabled(config, "NuclMapping") == [False, True]   # VFDB is still mapped
    assert "CARD" in caplog.text and "NuclSearch" in caplog.text


def test_genetic_code_propagates():
    config = _config(); set_genetic_code(config, 4)
    assert [s["options"]["transl_table"] for s in config.FUNCTIONAL_ANNOTATION if s["component_name"] == "NuclMapping"] == [4, 4, 4]

def test_mapping_option_defaults():
    for path in ["dfc/default_config.py", "example/dfc_config/dfast_web_config.py"]:
        options = [s["options"] for s in _config(path).FUNCTIONAL_ANNOTATION if s["component_name"] == "NuclMapping"]
        assert [(o["blast_task"], o["pseudogene_pident_cutoff"], o["overlap_fraction"], o["truncated_fraction"])
                for o in options] == [("dc-megablast", 95, 0.2, 0.7)] * 3

def test_amr_keeps_other_nuclsearch_entries():
    config = _config()
    custom = {"component_name": "NuclSearch", "enabled": True,   # a user's own NuclSearch entry
              "options": {"database": "/my/db/MY_DB", "db_name": "MY_DB"}}
    config.FUNCTIONAL_ANNOTATION.append(custom)
    enable_amr(config)
    assert custom["enabled"] and _enabled(config, "NuclMapping") == [True, True]


def test_custom_nucleotide_reference_uses_nuclmapping():
    for path in ["dfc/default_config.py", "example/dfc_config/dfast_web_config.py"]:
        config = _config(path)
        assert _enabled(config, "NuclSearch", ("PLASMID_DB",)) == []
        assert _enabled(config, "NuclMapping", ("PLASMID_DB",)) == [False]
        enable_amr(config)                      # --amr maps CARD/VFDB only
        assert _enabled(config, "NuclMapping", ("PLASMID_DB",)) == [False]
