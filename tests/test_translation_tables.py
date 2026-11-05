from logging import getLogger
from types import SimpleNamespace

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from dfc.models.bio_feature import ExtendedFeature
from dfc.structuralAnnotation import StructuralAnnotation, start_codons_of


class FakeCDSTool:
    TYPE = "CDS"
    NAME = "Fake"

    def __init__(self, transl_table, features):
        self.transl_table = transl_table
        self._features = features

    def getFeatures(self):
        return {"seq1": self._features}


def translate_cds(nuc, transl_table):
    """Run setFeatures with a fake CDS tool and return the translation of the single CDS."""
    feature = ExtendedFeature(FeatureLocation(0, len(nuc), 1), type="CDS", id="cds1")
    feature.qualifiers["transl_table"] = [str(transl_table)]
    tool = FakeCDSTool(transl_table, [feature])
    record = SimpleNamespace(seq=Seq(nuc), features=[])
    genome = SimpleNamespace(seq_records={"seq1": record}, sort_features=lambda: None,
                             set_feature_dictionary=lambda: None)
    sa = StructuralAnnotation.__new__(StructuralAnnotation)
    sa.genome, sa.tools, sa.logger = genome, [tool], getLogger("test")
    sa.setFeatures()
    return feature.qualifiers["translation"][0]


def test_start_codons_table25():
    assert start_codons_of(25) == {"ATG", "GTG", "TTG"}


def test_start_codons_other_tables():
    assert "TTA" in start_codons_of(4)
    assert "ATT" in start_codons_of(11)
    assert start_codons_of("11") == start_codons_of(11)


def test_ctg_start_not_m_in_table25():
    nuc = "CTGAAATAA"
    assert not translate_cds(nuc, 25).startswith("M")
    assert translate_cds(nuc, 11).startswith("M")


def test_tga_is_gly_in_table25():
    assert translate_cds("ATGTGAAAATAA", 25) == "MGK"
