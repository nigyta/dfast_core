# coding: UTF8
"""HMMscan with the NCBI HMM collection: TIGR subset, naming by HMM attributes, and --hmm_db."""
from logging import getLogger

import pytest
from Bio.SeqFeature import FeatureLocation

from dfc.components.HMMscan import HMMscan, read_hmm_attributes, select_naming_hit
from dfc.models.bio_feature import ExtendedFeature
from dfc.models.hit import HmmHit
from dfc.utils.config_util import set_hmm_db
from dfc.utils.reffile_util import extract_hmm_models

HMM_LIB = """HMMER3/f [3.1b2 | February 2015]
NAME  rpmI_bact
ACC   TIGR00001.1
DESC  50S ribosomal protein L35
//
HMMER3/f [3.1b2 | February 2015]
NAME  vanH-D
ACC   NF000004.1
DESC  NCBIFAM: D-lactate dehydrogenase VanH-D
//
"""

# hmm_PGAP.tsv: 24 columns. 1 accession, 7 family_type, 9 for_naming, 11 product_name, 12 gene_symbol, 14 ec_number
def _row(acc, family_type, for_naming, product, gene="", ec=""):
    fields = [""] * 24
    fields[0], fields[6], fields[8], fields[10], fields[11], fields[13] = acc, family_type, for_naming, product, gene, ec
    return "\t".join(fields) + "\n"


ATTRIBUTES_TSV = ("#ncbi_accession\tsource_identifier\tlabel" + "\t" * 21 + "\n"
                  + _row("TIGR00001.1", "equivalog", "Y", "50S ribosomal protein L35", "rpmI")
                  + _row("NF000004.1", "subfamily", "Y", "D-lactate dehydrogenase VanH-D", "vanH-D", "1.1.1.28,1.1.1.-")
                  + _row("TIGR00002.1", "hypoth_equivalog", "Y", "YbaB/EbfC family nucleoid-associated protein")
                  + _row("TIGR00003.1", "equivalog", "Y", "hypothetical protein")
                  + _row("TIGR00004.1", "domain", "N", "DUF123 domain-containing protein"))


def test_extract_tigr_models(tmp_path):
    lib, out = tmp_path / "all.LIB", tmp_path / "tigr.LIB"
    lib.write_text(HMM_LIB)
    assert extract_hmm_models(str(lib), str(out), "TIGR") == 1
    assert "ACC   TIGR00001.1" in out.read_text() and "NF000004.1" not in out.read_text()


@pytest.fixture
def attributes(tmp_path):
    tsv = tmp_path / "hmm_PGAP.tsv"
    tsv.write_text(ATTRIBUTES_TSV)
    return read_hmm_attributes(str(tsv))


def _hit(acc, score, attributes):
    return HmmHit(acc, "name", "desc", 1e-50, score, 0.1, "TIGR", attributes.get(acc))


def test_read_attributes(attributes):
    assert attributes["NF000004.1"]["ec_number"] == "1.1.1.28,1.1.1.-"
    assert attributes["TIGR00001.1"]["family_type"] == "equivalog"


def test_naming_prefers_specific_family_type_then_score(attributes):
    hits = [_hit("NF000004.1", 900, attributes), _hit("TIGR00001.1", 100, attributes)]
    assert select_naming_hit(hits).accession == "TIGR00001.1"  # equivalog beats subfamily despite lower score


def test_naming_skips_hypothetical_and_not_for_naming(attributes):
    hits = [_hit("TIGR00003.1", 500, attributes), _hit("TIGR00004.1", 400, attributes), _hit("TIGR99999.1", 300, attributes)]
    assert select_naming_hit(hits) is None
    # hypoth_equivalog with a real family name can still name the protein
    assert select_naming_hit([_hit("TIGR00002.1", 50, attributes)]).accession == "TIGR00002.1"


def test_assign_sets_product_gene_ec_and_inference(attributes):
    feature = ExtendedFeature(location=FeatureLocation(0, 900, 1), type="CDS", id="MGA_1",
                              qualifiers={"product": ["hypothetical protein"]})
    _hit("NF000004.1", 900, attributes).assign(feature)
    q = feature.qualifiers
    assert q["product"] == ["D-lactate dehydrogenase VanH-D"] and q["gene"] == ["vanH-D"]
    assert q["EC_number"] == ["1.1.1.28", "1.1.1.-"]
    assert q["inference"] == ["protein motif:HMM:NF000004.1"]
    assert q["note"][0].startswith("TIGR:NF000004.1; desc")


def test_hit_without_attributes_is_only_a_note():
    feature = ExtendedFeature(location=FeatureLocation(0, 900, 1), type="CDS", id="MGA_1",
                              qualifiers={"product": ["hypothetical protein"]})
    HmmHit("TIGR00001", "rpmI_bact", "desc", 1e-50, 100, 0.1, "TIGR").assign(feature)
    assert feature.qualifiers["product"] == ["hypothetical protein"] and "inference" not in feature.qualifiers


def test_set_results_names_unannotated_cds_only(tmp_path, attributes):
    # hmmscan --tblout lines: target name, accession, query name, accession, E-value, score, bias, ... description
    tbl = ("rpmI_bact TIGR00001.1 MGA_1 - 1e-50 100.0 0.1 1e-50 100.0 0.1 1.0 1 0 0 1 1 1 1 50S ribosomal protein L35\n"
           "vanH-D NF000004.1 MGA_1 - 1e-60 900.0 0.1 1e-60 900.0 0.1 1.0 1 0 0 1 1 1 1 NCBIFAM: VanH-D\n"
           "rpmI_bact TIGR00001.1 MGA_2 - 1e-50 100.0 0.1 1e-50 100.0 0.1 1.0 1 0 0 1 1 1 1 50S ribosomal protein L35\n")
    (tmp_path / "result0.out").write_text(tbl)
    named, annotated = (ExtendedFeature(location=FeatureLocation(0, 900, 1), type="CDS", id=i) for i in ("MGA_1", "MGA_2"))
    annotated.primary_hit = "a database hit"
    search = HMMscan.__new__(HMMscan)  # skip tool and database checks in __init__
    search.logger, search.workDir, search.query_files, search.db_name = getLogger(__name__), str(tmp_path), {0: "q"}, "TIGR"
    search.attributes = attributes
    search.genome = type("FakeGenome", (), {"features": {"MGA_1": named, "MGA_2": annotated}})()

    search.set_results()

    assert named.primary_hit.accession == "TIGR00001.1"  # naming hit (equivalog)
    assert [h.accession for h in named.secondary_hits] == ["NF000004.1"]  # best-scoring hit kept as a note
    assert annotated.primary_hit == "a database hit"  # names from database hits are not replaced


class _Config:
    def __init__(self):
        self.FUNCTIONAL_ANNOTATION = [
            {"component_name": "HMMscan", "enabled": True, "options": {"db_name": "TIGR"}},
            {"component_name": "HMMscan", "enabled": False, "options": {"db_name": "NCBIfam"}},
            {"component_name": "HMMscan", "enabled": False, "options": {"db_name": "Pfam"}},
        ]


def test_set_hmm_db():
    config = _Config()
    set_hmm_db(config, "ncbifam")
    assert [s["enabled"] for s in config.FUNCTIONAL_ANNOTATION] == [False, True, False]
    set_hmm_db(config, "tigr")
    assert [s["enabled"] for s in config.FUNCTIONAL_ANNOTATION] == [True, False, False]
