# coding: UTF8
"""Tests for `dfast_file_downloader.py --mefinder` (retrieve_mefinder_reference).

The downloader is a script that parses argv and dispatches at import time, so the
function cannot be imported in isolation. Instead these tests run the script as a
subprocess with a fake `mefinder` executable placed first on PATH, exercising the
install-skip / --no_indexing / index-build / index-verification branches without
needing a real MobileElementFinder install or network access.
"""
import os
import sys
import stat
import subprocess

import pytest

APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOWNLOADER = os.path.join(APP_ROOT, "scripts", "dfast_file_downloader.py")

# Fake `mefinder`: always reports version 1.1.2. For `index --db-path DIR`,
# "good" mode creates the expected BLAST index files; "broken" mode creates none.
_FAKE_GOOD = r'''#!/bin/bash
if [ "$1" = "--version" ]; then echo "1.1.2"; exit 0; fi
if [ "$1" = "index" ]; then
  db=""
  while [ $# -gt 0 ]; do if [ "$1" = "--db-path" ]; then db="$2"; fi; shift; done
  mkdir -p "$db"
  for ext in nin nhr nsq ndb njs not ntf nto; do : > "$db/mge_records.$ext"; done
  echo "Created databases in: $db"; exit 0
fi
exit 0
'''

_FAKE_BROKEN = r'''#!/bin/bash
if [ "$1" = "--version" ]; then echo "1.1.2"; exit 0; fi
if [ "$1" = "index" ]; then echo "noop (no index created)"; exit 0; fi
exit 0
'''


def _make_fake_mefinder(bindir, script):
    os.makedirs(bindir, exist_ok=True)
    path = os.path.join(bindir, "mefinder")
    with open(path, "w") as f:
        f.write(script)
    os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def _run_downloader(bindir, dbroot, extra_args=()):
    env = dict(os.environ)
    env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
    cmd = [sys.executable, DOWNLOADER, "--mefinder", "-d", str(dbroot)] + list(extra_args)
    return subprocess.run(cmd, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, universal_newlines=True)


def test_mefinder_no_indexing_skips_index_build(tmp_path):
    # With --no_indexing, the index is not built (returns before makedirs).
    bindir = str(tmp_path / "bin")
    _make_fake_mefinder(bindir, _FAKE_GOOD)
    dbroot = tmp_path / "db"
    proc = _run_downloader(bindir, dbroot, extra_args=["--no_indexing"])
    assert proc.returncode == 0, proc.stderr
    assert not (dbroot / "mefinder_db" / "mge_records.nsq").exists()


def test_mefinder_builds_index(tmp_path):
    # Without --no_indexing, the BLAST index is built into DB_ROOT/mefinder_db.
    bindir = str(tmp_path / "bin")
    _make_fake_mefinder(bindir, _FAKE_GOOD)
    dbroot = tmp_path / "db"
    proc = _run_downloader(bindir, dbroot)
    assert proc.returncode == 0, proc.stderr
    for ext in ("nin", "nhr", "nsq"):
        assert (dbroot / "mefinder_db" / ("mge_records." + ext)).exists()


def test_mefinder_incomplete_index_fails(tmp_path):
    # If `mefinder index` produces no index files, verification fails (exit != 0).
    bindir = str(tmp_path / "bin")
    _make_fake_mefinder(bindir, _FAKE_BROKEN)
    dbroot = tmp_path / "db"
    proc = _run_downloader(bindir, dbroot)
    assert proc.returncode != 0
    assert "incomplete" in (proc.stdout + proc.stderr).lower()
