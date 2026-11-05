import os
import subprocess
import sys

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "dfast_file_downloader.py")


def run(*args):
    return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)


def test_downloader_modules_integron(tmp_path):
    r = run("--modules", "integron", "--dbroot", str(tmp_path))
    assert r.returncode == 0
    assert "IntegronFinder needs no reference data" in r.stdout + r.stderr
    assert list(tmp_path.iterdir()) == []


def test_downloader_modules_help_lists_names():
    r = run("-h")
    assert r.returncode == 0
    assert "--modules" in r.stdout
    assert "integron" in r.stdout


def test_downloader_modules_bad_name():
    assert run("--modules", "foo").returncode == 2
