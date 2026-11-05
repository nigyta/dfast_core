# coding: UTF8
"""StructuralAnnotation config handling: disabled entries and removed tools."""
import pytest

from dfc.structuralAnnotation import StructuralAnnotation


class _Config:
    CPU = 1

    def __init__(self, tmp_path, entries):
        self.WORK_DIR = str(tmp_path)
        self.STRUCTURAL_ANNOTATION = entries


def _entry(name, enabled):
    return {"tool_name": name, "enabled": enabled, "options": {}}


def test_disabled_removed_tools_are_ignored(tmp_path):
    # a Ver 1.5 default_config.py still lists CRT and RNAmmer (disabled)
    sa = StructuralAnnotation(None, _Config(tmp_path, [_entry("CRT", False), _entry("RNAmmer", False)]))
    assert sa.tools == []


@pytest.mark.parametrize("name", ["CRT", "RNAmmer"])
def test_enabled_removed_tool_exits(tmp_path, caplog, name):
    with pytest.raises(SystemExit):
        StructuralAnnotation(None, _Config(tmp_path, [_entry(name, True)]))
    assert "removed in Ver 1.6" in caplog.text and name in caplog.text


def test_enabled_unknown_tool_exits(tmp_path, caplog):
    with pytest.raises(SystemExit):
        StructuralAnnotation(None, _Config(tmp_path, [_entry("Nonexistent", True)]))
    assert "not registered" in caplog.text
