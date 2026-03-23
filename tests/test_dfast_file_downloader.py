#! /usr/bin/env python
# coding: UTF8

import os
import sys
import tempfile
import shutil
import unittest
from unittest.mock import patch, MagicMock

# The downloader script runs argparse at module level. We must
# set sys.argv before importing so argparse doesn't consume test args.
# With no flags, argparse defaults are all None/False — the "all is None"
# check passes over --plasmidfinder (False != None), so no download runs.
_original_argv = sys.argv
sys.argv = [sys.argv[0]]

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scripts import dfast_file_downloader

sys.argv = _original_argv


class TestIndexPlasmidfinerWithBlast(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        for name in ["enterobacteriales.fsa", "Inc18.fsa", "Rep1.fsa"]:
            with open(os.path.join(self.tmpdir, name), "w") as f:
                f.write(">dummy\nACGT\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    @patch("subprocess.Popen")
    def test_indexes_all_fsa_files(self, mock_popen):
        proc = MagicMock()
        proc.communicate.return_value = (b"", b"")
        proc.returncode = 0
        mock_popen.return_value = proc

        dfast_file_downloader._index_plasmidfinder_with_blast(self.tmpdir)

        self.assertEqual(mock_popen.call_count, 3)
        called_fsa = set()
        for c in mock_popen.call_args_list:
            cmd = c[0][0]
            self.assertEqual(cmd[0], "makeblastdb")
            self.assertIn("-dbtype", cmd)
            self.assertIn("nucl", cmd)
            self.assertIn("-hash_index", cmd)
            self.assertIn("-parse_seqids", cmd)
            in_index = cmd.index("-in")
            called_fsa.add(os.path.basename(cmd[in_index + 1]))
        self.assertEqual(called_fsa, {"enterobacteriales.fsa", "Inc18.fsa", "Rep1.fsa"})

    @patch("subprocess.Popen")
    def test_output_db_name_matches_fsa_stem(self, mock_popen):
        proc = MagicMock()
        proc.communicate.return_value = (b"", b"")
        proc.returncode = 0
        mock_popen.return_value = proc

        dfast_file_downloader._index_plasmidfinder_with_blast(self.tmpdir)

        for c in mock_popen.call_args_list:
            cmd = c[0][0]
            in_index = cmd.index("-in")
            out_index = cmd.index("-out")
            fsa_path = cmd[in_index + 1]
            db_path = cmd[out_index + 1]
            self.assertEqual(db_path, fsa_path.replace(".fsa", ""))

    @patch("shutil.which", return_value=None)
    def test_exits_if_makeblastdb_not_found(self, mock_which):
        with self.assertRaises(SystemExit):
            dfast_file_downloader._index_plasmidfinder_with_blast(self.tmpdir)

    def test_exits_if_no_fsa_files(self):
        empty_dir = tempfile.mkdtemp()
        try:
            with self.assertRaises(SystemExit):
                dfast_file_downloader._index_plasmidfinder_with_blast(empty_dir)
        finally:
            shutil.rmtree(empty_dir, ignore_errors=True)

    @patch("subprocess.Popen")
    def test_exits_on_makeblastdb_failure(self, mock_popen):
        proc = MagicMock()
        proc.communicate.return_value = (b"", b"error occurred")
        proc.returncode = 1
        mock_popen.return_value = proc

        with self.assertRaises(SystemExit):
            dfast_file_downloader._index_plasmidfinder_with_blast(self.tmpdir)


class TestRetrievePlasmidfinder(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.db_dir = os.path.join(self.tmpdir, "plasmidfinder_db")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    @patch("scripts.dfast_file_downloader._index_plasmidfinder_with_blast")
    @patch("subprocess.Popen")
    @patch("shutil.which", return_value=None)
    def test_falls_back_to_blast_when_no_kma(self, mock_which, mock_popen, mock_blast_index):
        proc = MagicMock()
        proc.communicate.return_value = (b"", b"")
        proc.returncode = 0
        mock_popen.return_value = proc

        os.makedirs(self.db_dir, exist_ok=True)

        dfast_file_downloader.retrieve_plasmidfinder_reference(self.tmpdir)

        mock_blast_index.assert_called_once_with(self.db_dir)

    @patch("subprocess.Popen")
    @patch("shutil.which", return_value="/usr/bin/kma")
    def test_uses_kma_when_available(self, mock_which, mock_popen):
        proc = MagicMock()
        proc.communicate.return_value = (b"", b"")
        proc.returncode = 0
        mock_popen.return_value = proc

        os.makedirs(self.db_dir, exist_ok=True)

        dfast_file_downloader.retrieve_plasmidfinder_reference(self.tmpdir)

        # Two Popen calls: git clone + INSTALL.py
        self.assertEqual(mock_popen.call_count, 2)
        install_cmd = mock_popen.call_args_list[1][0][0]
        self.assertIn("INSTALL.py", install_cmd)

    @patch("scripts.dfast_file_downloader._index_plasmidfinder_with_blast")
    @patch("subprocess.Popen")
    @patch("shutil.which", return_value="/usr/bin/kma")
    def test_falls_back_to_blast_on_kma_failure(self, mock_which, mock_popen, mock_blast_index):
        clone_proc = MagicMock()
        clone_proc.communicate.return_value = (b"", b"")
        clone_proc.returncode = 0

        install_proc = MagicMock()
        install_proc.communicate.return_value = (b"", b"kma error")
        install_proc.returncode = 1

        mock_popen.side_effect = [clone_proc, install_proc]

        os.makedirs(self.db_dir, exist_ok=True)

        dfast_file_downloader.retrieve_plasmidfinder_reference(self.tmpdir)

        mock_blast_index.assert_called_once_with(self.db_dir)

    @patch("subprocess.Popen")
    def test_exits_on_clone_failure(self, mock_popen):
        proc = MagicMock()
        proc.communicate.return_value = (b"", b"clone failed")
        proc.returncode = 1
        mock_popen.return_value = proc

        with self.assertRaises(SystemExit):
            dfast_file_downloader.retrieve_plasmidfinder_reference(self.tmpdir)

    @patch("subprocess.Popen")
    @patch("shutil.which", return_value=None)
    def test_removes_existing_db_dir(self, mock_which, mock_popen):
        os.makedirs(self.db_dir, exist_ok=True)
        old_file = os.path.join(self.db_dir, "old_file.txt")
        with open(old_file, "w") as f:
            f.write("old")

        proc = MagicMock()
        proc.communicate.return_value = (b"", b"")
        proc.returncode = 0
        mock_popen.return_value = proc

        try:
            dfast_file_downloader.retrieve_plasmidfinder_reference(self.tmpdir)
        except (SystemExit, FileNotFoundError):
            pass

        self.assertFalse(os.path.exists(old_file))


class TestDownloadPlasmifinderExtraDatabases(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create a minimal config file
        with open(os.path.join(self.tmpdir, "config"), "w") as f:
            f.write("#db_prefix\tname\tdescription\n")
            f.write("enterobacteriales\tEnterobacteriales\tEnterobacteriales\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    @patch("scripts.dfast_file_downloader.request.urlretrieve")
    def test_downloads_and_registers_extra_databases(self, mock_urlretrieve):
        def fake_download(url, output_file):
            with open(output_file, "w") as f:
                f.write(">seq1\nACGT\n")

        mock_urlretrieve.side_effect = fake_download

        dfast_file_downloader._download_plasmidfinder_extra_databases(self.tmpdir)

        # Check files were "downloaded"
        self.assertTrue(os.path.exists(os.path.join(self.tmpdir, "repP_database_v2.fsa")))
        self.assertTrue(os.path.exists(os.path.join(self.tmpdir, "AcinetobacterPlasmidTyping_v3.fsa")))

        # Check config was updated
        with open(os.path.join(self.tmpdir, "config")) as f:
            config_content = f.read()
        self.assertIn("repP_database_v2", config_content)
        self.assertIn("AcinetobacterPlasmidTyping_v3", config_content)
        self.assertIn("repP", config_content)
        self.assertIn("Acinetobacter", config_content)

    @patch("scripts.dfast_file_downloader.request.urlretrieve")
    def test_downloads_correct_urls(self, mock_urlretrieve):
        def fake_download(url, output_file):
            with open(output_file, "w") as f:
                f.write(">seq1\nACGT\n")

        mock_urlretrieve.side_effect = fake_download

        dfast_file_downloader._download_plasmidfinder_extra_databases(self.tmpdir)

        called_urls = [call[0][0] for call in mock_urlretrieve.call_args_list]
        self.assertIn("https://ndownloader.figshare.com/files/58979998", called_urls)
        self.assertIn("https://ndownloader.figshare.com/files/58980199", called_urls)

    @patch("scripts.dfast_file_downloader.request.urlretrieve")
    def test_skips_on_download_failure(self, mock_urlretrieve):
        mock_urlretrieve.side_effect = Exception("Network error")

        # Should not raise — failures are logged and skipped
        dfast_file_downloader._download_plasmidfinder_extra_databases(self.tmpdir)

        # Config should have no extra entries
        with open(os.path.join(self.tmpdir, "config")) as f:
            lines = [l for l in f.readlines() if not l.startswith("#")]
        self.assertEqual(len(lines), 1)  # only original enterobacteriales

    @patch("scripts.dfast_file_downloader.request.urlretrieve")
    def test_skips_empty_download(self, mock_urlretrieve):
        def fake_empty_download(url, output_file):
            with open(output_file, "w") as f:
                pass  # empty file

        mock_urlretrieve.side_effect = fake_empty_download

        dfast_file_downloader._download_plasmidfinder_extra_databases(self.tmpdir)

        with open(os.path.join(self.tmpdir, "config")) as f:
            lines = [l for l in f.readlines() if not l.startswith("#")]
        self.assertEqual(len(lines), 1)


if __name__ == "__main__":
    unittest.main()
