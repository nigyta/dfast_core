# coding: UTF8
"""PlasmidFinder 3.x from PyPI has no command: DFAST runs it as a module unless plasmidfinder.py is on PATH."""
import os
import subprocess
import sys

APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = "import dfc.tools.dfast_plasmidfinder as m; print(m.PLASMIDFINDER)"


def _command(path):
    env = dict(os.environ, PATH=path)
    return subprocess.run([sys.executable, "-c", CODE], cwd=APP_ROOT, env=env,
                          stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout.strip()


def test_module_when_no_command(tmp_path):
    assert _command(str(tmp_path)) == str([sys.executable, "-m", "plasmidfinder"])


def test_command_on_path(tmp_path):
    script = tmp_path / "plasmidfinder.py"
    script.write_text("#!/bin/sh\n")
    script.chmod(0o755)
    assert _command(str(tmp_path)) == "['plasmidfinder.py']"


def _fake(tmp_path, body):
    script = tmp_path / "fake_plasmidfinder"
    script.write_text("#!/bin/sh\n" + body + "\n")
    script.chmod(0o755)
    return [str(script)]


def test_version_3_is_accepted(tmp_path, monkeypatch):
    import pytest
    from dfc.tools.dfast_plasmidfinder import Plasmidfinder
    monkeypatch.setattr(Plasmidfinder, "version", None)
    monkeypatch.setattr(Plasmidfinder, "VERSION_CHECK_CMD", _fake(tmp_path, "echo 3.0.3"))
    Plasmidfinder(workDir=str(tmp_path))
    assert Plasmidfinder.version == "3.0.3"
    # PlasmidFinder 2.x has no -v option (usage error), and a missing module fails: both abort
    for body in ("echo 'error: the following arguments are required: -i/--infile' >&2; exit 2", "echo 2.1.6"):
        monkeypatch.setattr(Plasmidfinder, "version", None)
        monkeypatch.setattr(Plasmidfinder, "VERSION_CHECK_CMD", _fake(tmp_path, body))
        with pytest.raises(SystemExit):
            Plasmidfinder(workDir=str(tmp_path))
