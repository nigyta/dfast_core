import logging
import os
import subprocess
import sys

import pytest

from dfc.utils import config_util
from dfc.utils.config_util import load_config, enable_modules, MODULES

APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DFAST = os.path.join(APP_ROOT, "dfast")
GENOME = os.path.join(APP_ROOT, "example", "sample.genome.fna")


def _entries(config):
    return list(config.FUNCTIONAL_ANNOTATION) + list(config.CONTIG_ANNOTATION)


def _entry(config, name, db_name=None):
    for s in _entries(config):
        if name in (s.get("component_name"), s.get("tool_name")) and (db_name is None or s["options"].get("db_name") == db_name):
            return s


def _enabled(config, name, db_name=None):
    return _entry(config, name, db_name)["enabled"]


def _config(tmp_path, available=True):
    """Default config whose data paths point into tmp_path (created when `available`)."""
    config = load_config(APP_ROOT, os.path.join(APP_ROOT, "dfc", "default_config.py"))
    for name, db_name, key in [("NuclMapping", "CARD", "database"), ("NuclMapping", "VFDB", "database"),
                               ("Plasmidfinder", None, "db_path"), ("MobileElementFinder", None, "db_path"),
                               ("DefenseFinder", None, "models_dir")]:
        entry = _entry(config, name, db_name)
        path = tmp_path / (db_name or name)
        entry["options"][key] = str(path)
        if available:
            if key == "database":
                for suffix in (".nucl.ref", ".nucl.fasta"):
                    open(str(path) + suffix, "w").close()
            else:
                path.mkdir()
    return config


@pytest.fixture
def tools(monkeypatch):
    monkeypatch.setattr(config_util.shutil, "which", lambda cmd: "/usr/bin/" + cmd)
    monkeypatch.setattr(config_util, "_plasmidfinder_available", lambda: True)


@pytest.fixture
def no_tools(monkeypatch):
    monkeypatch.setattr(config_util.shutil, "which", lambda cmd: None)
    monkeypatch.setattr(config_util, "_plasmidfinder_available", lambda: False)


def test_explicit_modules_enable_entries(tmp_path, tools):
    config = _config(tmp_path)
    assert enable_modules(config, ["integron", "amr"]) == ["amr", "integron"]
    assert _enabled(config, "NuclMapping", "CARD") and _enabled(config, "IntegronFinder")
    for name, db_name in [("NuclMapping", "VFDB"), ("Plasmidfinder", None), ("MobileElementFinder", None), ("DefenseFinder", None)]:
        assert not _enabled(config, name, db_name)


def test_all_enables_six(tmp_path, tools):
    config = _config(tmp_path)
    assert enable_modules(config, ["all"]) == list(MODULES)
    for name, db_name in [("NuclMapping", "CARD"), ("NuclMapping", "VFDB"), ("Plasmidfinder", None),
                          ("MobileElementFinder", None), ("IntegronFinder", None), ("DefenseFinder", None)]:
        assert _enabled(config, name, db_name)


def test_all_skips_missing(tmp_path, no_tools, caplog):
    config = _config(tmp_path, available=False)
    with caplog.at_level(logging.INFO):
        assert enable_modules(config, ["all"]) == []
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == len(MODULES)
    assert not any(s["enabled"] for s in _entries(config) if s.get("tool_name") in ("Plasmidfinder", "MobileElementFinder", "IntegronFinder"))


def test_explicit_missing_exits(tmp_path, tools, caplog):
    config = _config(tmp_path)
    config_util_models = _entry(config, "DefenseFinder")["options"]["models_dir"]
    os.rmdir(config_util_models)
    with caplog.at_level(logging.INFO), pytest.raises(SystemExit) as e:
        enable_modules(config, ["defense"])
    assert e.value.code == 1
    assert "dfast_file_downloader.py --modules defense" in caplog.text


def test_explicit_with_all_is_strict(tmp_path, tools):
    config = _config(tmp_path)
    os.rmdir(_entry(config, "DefenseFinder")["options"]["models_dir"])
    with pytest.raises(SystemExit):
        enable_modules(config, ["all", "defense"])


def test_missing_executable_message(tmp_path, no_tools, caplog):
    config = _config(tmp_path)
    with caplog.at_level(logging.INFO), pytest.raises(SystemExit):
        enable_modules(config, ["integron"])
    assert "conda install -c bioconda integron_finder" in caplog.text


def test_amr_respects_nuclsearch_in_config(tmp_path, tools, caplog):
    config = _config(tmp_path, available=False)    # no CARD files: the NuclSearch entry is used instead
    _entry(config, "NuclSearch", "CARD")["enabled"] = True
    with caplog.at_level(logging.INFO):
        assert enable_modules(config, ["amr"]) == ["amr"]
    assert not _enabled(config, "NuclMapping", "CARD")
    assert _enabled(config, "NuclSearch", "CARD")
    assert "NuclSearch" in caplog.text


@pytest.mark.parametrize("flag,modules", [("--amr", "--modules amr vf plasmid defense"), ("--mge", "--modules mge integron")])
def test_removed_flags(flag, modules):
    p = subprocess.run([sys.executable, DFAST, "-g", GENOME, flag, "--show_config"], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True)
    assert p.returncode == 1
    assert "{} was removed in Ver 1.6. Use {} (or --modules all).".format(flag, modules) in p.stdout


def test_modules_choices():
    p = subprocess.run([sys.executable, DFAST, "-g", GENOME, "--modules", "foo", "--show_config"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert p.returncode == 2
