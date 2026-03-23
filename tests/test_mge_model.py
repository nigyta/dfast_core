#! /usr/bin/env python
# coding: UTF8

import os
import sys
import tempfile
import shutil
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dfc.models.mge import MGE


SAMPLE_FASTA = """\
>Tn1|1|GQ160960.2
GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGCAAC
>ISAba1|3|CP000521.1
ATGCGTTATATTCGCCTGTGTATTATCTCCCTGTTAGCCGCC
>ICEPaePA14-1|1|NC_008463.1
ATGATAGGTTTGATTGTTGCGAGGTCAAAGAATAATGTTATAG
"""

SAMPLE_TSV = (
    "Tn1|1|GQ160960.2\tTn1\t1\tGQ160960.2\tGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGCAAC\t\n"
    "ISAba1|3|CP000521.1\tISAba1\t3\tCP000521.1\tATGCGTTATATTCGCCTGTGTATTATCTCCCTGTTAGCCGCC\t\n"
)


class TestMGEModel(unittest.TestCase):

    def test_init(self):
        mge = MGE("Tn1|1|GQ160960.2", "Tn1", "1", "GQ160960.2", "ACGT")
        self.assertEqual(mge.mge_id, "Tn1|1|GQ160960.2")
        self.assertEqual(mge.name, "Tn1")
        self.assertEqual(mge.allele, "1")
        self.assertEqual(mge.accession, "GQ160960.2")
        self.assertEqual(mge.nucl_seq, "ACGT")

    def test_str(self):
        mge = MGE("Tn1|1|GQ160960.2", "Tn1", "1", "GQ160960.2", "ACGT")
        self.assertIn("Tn1", str(mge))

    def test_info(self):
        mge = MGE("Tn1|1|GQ160960.2", "Tn1", "1", "GQ160960.2", "ACGT")
        info = mge.info()
        self.assertIn("Tn1", info)
        self.assertIn("allele 1", info)
        self.assertIn("GQ160960.2", info)

    def test_to_tabular(self):
        mge = MGE("Tn1|1|GQ160960.2", "Tn1", "1", "GQ160960.2", "ACGT", "test note")
        tab = mge.to_tabular()
        fields = tab.split("\t")
        self.assertEqual(fields[0], "Tn1|1|GQ160960.2")
        self.assertEqual(fields[1], "Tn1")
        self.assertEqual(fields[2], "1")
        self.assertEqual(fields[3], "GQ160960.2")
        self.assertEqual(fields[4], "ACGT")
        self.assertEqual(fields[5], "test note")

    def test_to_nucl_fasta(self):
        mge = MGE("Tn1|1|GQ160960.2", "Tn1", "1", "GQ160960.2", "ACGT")
        fasta = mge.to_nucl_fasta()
        self.assertEqual(fasta, ">Tn1|1|GQ160960.2\nACGT\n")

    def test_to_prot_fasta_empty(self):
        mge = MGE("Tn1|1|GQ160960.2", "Tn1", "1", "GQ160960.2", "ACGT")
        self.assertEqual(mge.to_prot_fasta(), "")

    def test_parse_line(self):
        line = "Tn1|1|GQ160960.2\tTn1\t1\tGQ160960.2\tACGT\tsome note\n"
        mge_id, mge = MGE.parse_line(line)
        self.assertEqual(mge_id, "Tn1|1|GQ160960.2")
        self.assertEqual(mge.name, "Tn1")
        self.assertEqual(mge.allele, "1")
        self.assertEqual(mge.accession, "GQ160960.2")
        self.assertEqual(mge.nucl_seq, "ACGT")
        self.assertEqual(mge.note, "some note")

    def test_parse_line_no_note(self):
        line = "ISAba1|3|CP000521.1\tISAba1\t3\tCP000521.1\tATGC\n"
        mge_id, mge = MGE.parse_line(line)
        self.assertEqual(mge.note, "")

    def test_parse_mge_fasta(self):
        tmpdir = tempfile.mkdtemp()
        try:
            fasta_file = os.path.join(tmpdir, "test.fna")
            with open(fasta_file, "w") as f:
                f.write(SAMPLE_FASTA)

            result = MGE.parse_mge_fasta(fasta_file)
            self.assertEqual(len(result), 3)
            self.assertIn("Tn1|1|GQ160960.2", result)
            self.assertIn("ISAba1|3|CP000521.1", result)
            self.assertIn("ICEPaePA14-1|1|NC_008463.1", result)

            tn1 = result["Tn1|1|GQ160960.2"]
            self.assertEqual(tn1.name, "Tn1")
            self.assertEqual(tn1.allele, "1")
            self.assertEqual(tn1.accession, "GQ160960.2")
        finally:
            shutil.rmtree(tmpdir)

    def test_to_dict(self):
        mge = MGE("Tn1|1|GQ160960.2", "Tn1", "1", "GQ160960.2", "ACGT")
        d = mge.to_dict()
        self.assertEqual(d["gene"], "Tn1")
        self.assertIn("mobile element", d["product"])
        self.assertIn("MGE:", d["accession"])
        self.assertIn("Tn1", d["note"])

    def test_set_info_to_feature(self):
        mge = MGE("ISAba1|3|CP000521.1", "ISAba1", "3", "CP000521.1", "ACGT")

        class MockFeature:
            qualifiers = {}

        feature = MockFeature()
        mge.set_info_to_feature(feature)
        self.assertIn("mobile element ISAba1", feature.qualifiers["product"])
        notes = feature.qualifiers["note"]
        self.assertTrue(any("MGEdb" in n for n in notes))

    def test_load_from_tsv_file(self):
        tmpdir = tempfile.mkdtemp()
        try:
            tsv_file = os.path.join(tmpdir, "test.nucl.ref")
            with open(tsv_file, "w") as f:
                f.write(SAMPLE_TSV)

            result = MGE.load_from_tsv_file(tsv_file)
            self.assertEqual(len(result), 2)
            self.assertIn("Tn1|1|GQ160960.2", result)
            self.assertIn("ISAba1|3|CP000521.1", result)
        finally:
            shutil.rmtree(tmpdir)


if __name__ == "__main__":
    unittest.main()
