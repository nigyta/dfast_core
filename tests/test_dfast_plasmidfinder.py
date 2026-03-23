#! /usr/bin/env python
# coding: UTF8

import json
import os
import tempfile
import unittest
from unittest.mock import patch, MagicMock

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dfc.tools.dfast_plasmidfinder import Plasmidfinder


# Sample PlasmidFinder 3.x BeOne JSON output (single hit)
SAMPLE_BEONE_SINGLE = {
    "type": "software_result",
    "databases": {
        "PlasmidFinder-2.2.0": {
            "type": "database",
            "database_name": "PlasmidFinder",
            "database_version": "2.2.0",
        }
    },
    "seq_regions": {
        "sequence1:55802..56462:IncL_1__JN626286:100.000000": {
            "type": "seq_region",
            "ref_database": "PlasmidFinder-2.2.0:enterobacteriales",
            "name": "IncL",
            "identity": 100.0,
            "alignment_length": 661,
            "ref_seq_lenght": 661,
            "coverage": 100.0,
            "ref_id": "IncL_1__JN626286",
            "ref_acc": "JN626286",
            "query_id": "sequence1",
            "query_start_pos": 55802,
            "query_end_pos": 56462,
            "note": "",
        }
    },
    "seq_variations": {},
    "phenotypes": {},
    "software_name": "PlasmidFinder",
    "software_version": "3.0.2",
}

# Sample with multiple hits on different contigs
SAMPLE_BEONE_MULTI = {
    "type": "software_result",
    "databases": {},
    "seq_regions": {
        "sequence1:100..500:IncFIA_1__AP001918:99.5": {
            "type": "seq_region",
            "ref_database": "PlasmidFinder-2.2.0:enterobacteriales",
            "name": "IncFIA",
            "identity": 99.5,
            "alignment_length": 388,
            "ref_seq_lenght": 400,
            "coverage": 97.0,
            "ref_acc": "AP001918",
            "query_id": "sequence1",
            "query_start_pos": 100,
            "query_end_pos": 500,
            "note": "",
        },
        "sequence2:2000..3000:IncN_1__AY046276:95.0": {
            "type": "seq_region",
            "ref_database": "PlasmidFinder-2.2.0:enterobacteriales",
            "name": "IncN",
            "identity": 95.0,
            "alignment_length": 950,
            "ref_seq_lenght": 1000,
            "coverage": 95.0,
            "ref_acc": "AY046276",
            "query_id": "sequence2",
            "query_start_pos": 2000,
            "query_end_pos": 3000,
            "note": "",
        },
    },
    "software_name": "PlasmidFinder",
    "software_version": "3.0.2",
}

# Sample with no hits
SAMPLE_BEONE_EMPTY = {
    "type": "software_result",
    "databases": {},
    "seq_regions": {},
    "software_name": "PlasmidFinder",
    "software_version": "3.0.2",
}

# Sample with corrected typo key (ref_seq_length instead of ref_seq_lenght)
SAMPLE_BEONE_CORRECTED_KEY = {
    "type": "software_result",
    "databases": {},
    "seq_regions": {
        "contig1:10..200:IncX_1__AB123456:98.0": {
            "type": "seq_region",
            "ref_database": "PlasmidFinder-2.2.0:enterobacteriales",
            "name": "IncX",
            "identity": 98.0,
            "alignment_length": 190,
            "ref_seq_length": 200,
            "coverage": 95.0,
            "ref_acc": "AB123456",
            "query_id": "contig1",
            "query_start_pos": 10,
            "query_end_pos": 200,
            "note": "",
        },
    },
    "software_name": "PlasmidFinder",
    "software_version": "3.1.0",
}


def _create_plasmidfinder(tmpdir, options=None):
    """Helper to create a Plasmidfinder instance with mocked version check."""
    if options is None:
        options = {"db_path": "/mock/db/plasmidfinder_db"}
    Plasmidfinder.version = "3.0.2"
    pf = Plasmidfinder(options=options, workDir=tmpdir)
    return pf


class TestPlasmidfinder(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Reset class-level version before each test
        Plasmidfinder.version = None

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # --- Initialization tests ---

    def test_init_creates_output_directory(self):
        pf = _create_plasmidfinder(self.tmpdir)
        expected_dir = os.path.join(self.tmpdir, "contig_annotation", "plasmidfinder")
        self.assertTrue(os.path.isdir(expected_dir))

    def test_init_sets_result_file(self):
        pf = _create_plasmidfinder(self.tmpdir)
        expected = os.path.join(self.tmpdir, "contig_annotation", "plasmidfinder", "results.json")
        self.assertEqual(pf.result_file, expected)

    def test_init_default_options(self):
        pf = _create_plasmidfinder(self.tmpdir, options={})
        self.assertEqual(pf.db_path, "")
        self.assertEqual(pf.method_path, "")
        self.assertEqual(pf.cmd_options, "")

    def test_init_with_method_path(self):
        pf = _create_plasmidfinder(self.tmpdir, options={
            "db_path": "/mock/db",
            "method_path": "/usr/bin/blastn",
        })
        self.assertEqual(pf.method_path, "/usr/bin/blastn")

    # --- getCommand tests ---

    def test_get_command_basic(self):
        pf = _create_plasmidfinder(self.tmpdir, options={
            "db_path": "/mock/db/plasmidfinder_db",
        })
        cmd = pf.getCommand()
        self.assertEqual(cmd[0], "plasmidfinder.py")
        self.assertIn("--infile", cmd)
        self.assertIn("--outputPath", cmd)
        self.assertIn("--tmp_dir", cmd)
        self.assertIn("--databasePath", cmd)
        self.assertIn("-j", cmd)
        # No -mp when method_path is empty
        self.assertNotIn("-mp", cmd)

    def test_get_command_with_method_path(self):
        pf = _create_plasmidfinder(self.tmpdir, options={
            "db_path": "/mock/db/plasmidfinder_db",
            "method_path": "/usr/local/bin/blastn",
        })
        cmd = pf.getCommand()
        mp_index = cmd.index("-mp")
        self.assertEqual(cmd[mp_index + 1], "/usr/local/bin/blastn")

    def test_get_command_infile_points_to_genome_fasta(self):
        pf = _create_plasmidfinder(self.tmpdir, options={
            "db_path": "/mock/db",
        })
        cmd = pf.getCommand()
        infile_index = cmd.index("--infile")
        expected_fasta = os.path.join(self.tmpdir, "input", "genome.fna")
        self.assertEqual(cmd[infile_index + 1], expected_fasta)

    def test_get_command_result_file_via_j_flag(self):
        pf = _create_plasmidfinder(self.tmpdir, options={
            "db_path": "/mock/db",
        })
        cmd = pf.getCommand()
        j_index = cmd.index("-j")
        self.assertEqual(cmd[j_index + 1], pf.result_file)

    # --- getResult tests ---

    def _write_result_file(self, pf, data):
        """Write JSON data to the result file path."""
        with open(pf.result_file, "w") as f:
            json.dump(data, f)

    def test_get_result_single_hit(self):
        pf = _create_plasmidfinder(self.tmpdir)
        self._write_result_file(pf, SAMPLE_BEONE_SINGLE)

        notes, report = pf.getResult()

        self.assertIn("sequence1", notes)
        note_list = notes["sequence1"]
        self.assertEqual(len(note_list), 3)
        self.assertEqual(note_list[0], "PlasmidFinder: IncL")
        self.assertIn("enterobacteriales", note_list[1])
        self.assertIn("IncL", note_list[1])
        self.assertIn("JN626286", note_list[2])
        self.assertIn("55802..56462", note_list[2])
        self.assertIn("100.0%", note_list[2])
        self.assertIn("661/661", note_list[2])

        self.assertIn("sequence1", report)
        self.assertEqual(len(report["sequence1"]), 1)

    def test_get_result_multiple_hits(self):
        pf = _create_plasmidfinder(self.tmpdir)
        self._write_result_file(pf, SAMPLE_BEONE_MULTI)

        notes, report = pf.getResult()

        self.assertIn("sequence1", notes)
        self.assertIn("sequence2", notes)
        # Each hit produces 3 note lines
        self.assertEqual(len(notes["sequence1"]), 3)
        self.assertEqual(len(notes["sequence2"]), 3)
        self.assertEqual(notes["sequence1"][0], "PlasmidFinder: IncFIA")
        self.assertEqual(notes["sequence2"][0], "PlasmidFinder: IncN")
        self.assertIn("99.5%", notes["sequence1"][2])
        self.assertIn("388/400", notes["sequence1"][2])
        self.assertIn("95.0%", notes["sequence2"][2])
        self.assertIn("950/1000", notes["sequence2"][2])

        self.assertEqual(len(report["sequence1"]), 1)
        self.assertEqual(len(report["sequence2"]), 1)

    def test_get_result_empty(self):
        pf = _create_plasmidfinder(self.tmpdir)
        self._write_result_file(pf, SAMPLE_BEONE_EMPTY)

        notes, report = pf.getResult()

        self.assertEqual(notes, {})
        self.assertEqual(report, {})

    def test_get_result_corrected_key_spelling(self):
        """PlasmidFinder 3.x has a typo 'ref_seq_lenght'; if a future version
        fixes it to 'ref_seq_length', the parser should still work."""
        pf = _create_plasmidfinder(self.tmpdir)
        self._write_result_file(pf, SAMPLE_BEONE_CORRECTED_KEY)

        notes, report = pf.getResult()

        self.assertIn("contig1", notes)
        self.assertIn("190/200", notes["contig1"][2])

    def test_get_result_organism_extraction(self):
        pf = _create_plasmidfinder(self.tmpdir)
        self._write_result_file(pf, SAMPLE_BEONE_SINGLE)

        notes, _ = pf.getResult()

        # The description should contain the organism extracted from ref_database
        description = notes["sequence1"][1]
        self.assertIn("enterobacteriales", description)

    def test_get_result_report_format(self):
        pf = _create_plasmidfinder(self.tmpdir)
        self._write_result_file(pf, SAMPLE_BEONE_SINGLE)

        _, report = pf.getResult()

        report_str = report["sequence1"][0]
        self.assertTrue(report_str.startswith("PlasmidFinder: IncL"))
        self.assertIn("Description:", report_str)
        self.assertIn("HitInfo:", report_str)

    # --- setVersion tests ---

    @patch.object(Plasmidfinder, "executeCommand")
    def test_set_version(self, mock_exec):
        mock_exec.return_value = (b"3.0.2", b"")
        Plasmidfinder.version = None
        pf = Plasmidfinder(options={"db_path": "/mock/db"}, workDir=self.tmpdir)
        self.assertEqual(Plasmidfinder.version, "3.0.2")

    @patch.object(Plasmidfinder, "executeCommand")
    def test_set_version_with_extra_output(self, mock_exec):
        mock_exec.return_value = (b"PlasmidFinder version 3.1.0\n", b"")
        Plasmidfinder.version = None
        pf = Plasmidfinder(options={}, workDir=self.tmpdir)
        self.assertEqual(Plasmidfinder.version, "3.1.0")

    # --- Class attributes ---

    def test_class_attributes(self):
        self.assertEqual(Plasmidfinder.NAME, "Plasmidfinder")
        self.assertEqual(Plasmidfinder.TYPE, "source")
        self.assertEqual(Plasmidfinder.VERSION_CHECK_CMD, ["plasmidfinder.py", "-v"])


if __name__ == "__main__":
    unittest.main()
