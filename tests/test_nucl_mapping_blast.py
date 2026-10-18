import json
import os
import random
import shutil
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from dfc.components.NuclMapping import NuclMapping
from dfc.models.bio_feature import ExtendedFeature
from dfc.models.card import CARD
from tests.test_nucl_mapping import FakeGenome

pytestmark = pytest.mark.skipif(shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
                                reason="BLAST+ is not on PATH")

CODONS = ["GCT", "GGT", "CTG", "AAA", "GAA", "CGT", "TTC"]
RNG = random.Random(1)
GENE = "ATG" + "".join(RNG.choice(CODONS) for _ in range(198)) + "TAA"   # 600 bp, no internal stop
FLANK = "CG" * 150   # 300 bp without a stop codon in any frame
L = len(FLANK)


def _translation(seq):
    return str(Seq(seq[:len(seq) // 3 * 3]).translate(table=11)).rstrip("*")


def _cds(cds_id, seq_id, start, end, translation, partial_flag="00"):
    return ExtendedFeature(location=FeatureLocation(start, end, strand=1), type="CDS", id=cds_id, seq_id=seq_id,
                           qualifiers={"translation": [translation], "transl_table": ["11"], "codon_start": [1]},
                           annotations={"partial_flag": partial_flag})


def _contigs():
    """(gene variant, predicted CDSs as (id, start0, end0, translation, partial_flag)) per contig, a to g."""
    b = GENE[:147] + GENE[148:]                                   # 1-bp deletion in codon 50
    c = GENE[:297] + "TAA" + GENE[300:]                           # codon 100 is a stop
    f = GENE[:597] + "GCT" + "GCT" * 9 + "TAA"                    # stop lost, next stop 30 bp later
    g = "GTG" + GENE[3:]
    e = GENE[:400]                                                # truncated at the contig end
    return {
        "a": (GENE, [("a", L + 30, L + 600, _translation(GENE[30:]), "00")]),
        "b": (b, [("b1", L, L + 150, "M", "00"), ("b2", L + 153, L + 599, "M", "00")]),
        "c": (c, [("c", L, L + 300, "M", "00")]),
        "d": (GENE, []),
        "e": (e, [("e", L, L + 400, _translation(e), "01")]),
        "f": (f, [("f", L, L + len(f), _translation(f), "00")]),
        "g": (g, [("g", L, L + 600, "M" + _translation(g)[1:], "00")]),
    }


@pytest.fixture(scope="module")
def mapped(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("nucl_mapping")
    ref = CARD("ARO:1", "OXA-48", "OXA beta-lactamase", "P1", "", "N1", GENE, "E. coli", "")
    (tmp_path / "REF.nucl.ref").write_text(ref.to_tabular() + "\n")
    (tmp_path / "REF.nucl.fasta").write_text(ref.to_nucl_fasta())
    (tmp_path / "REF.version").write_text("test\n")
    records, ids = [], {}
    for i, (case, (gene, cds)) in enumerate(_contigs().items(), 1):
        seq_id = ids[case] = f"sequence{i}"
        seq = FLANK + gene + (FLANK if case != "e" else "")
        rec = SeqRecord(Seq(seq), id=seq_id)
        rec.features = [_cds(cid, seq_id, s, e, tr, pf) for cid, s, e, tr, pf in cds]
        records.append(rec)
    os.makedirs(tmp_path / "input")
    (tmp_path / "input" / "genome.fna").write_text("".join(f">{r.id}\n{r.seq}\n" for r in records))
    genome = FakeGenome(records)
    comp = NuclMapping(genome, {"database": str(tmp_path / "REF"), "db_name": "CARD"}, str(tmp_path), 1)
    comp.run()
    hits = json.load(open(os.path.join(comp.workDir, "nucl_hits.json")))
    return genome, ids, hits


def _cdss(mapped, case):
    genome, ids, _ = mapped
    return [f for f in genome.seq_records[ids[case]].features if f.type == "CDS"]


def _span(f):
    return int(f.location.start), int(f.location.end)


def test_a_start_corrected(mapped):
    (f,) = _cdss(mapped, "a")
    assert _span(f) == (L, L + 600)
    assert "start codon moved 30 bp upstream by mapping of CARD:OXA-48" in f.qualifiers["note"]
    assert f.qualifiers["translation"] == [_translation(GENE)]


def test_b_frameshift_pseudogene(mapped):
    (f,) = _cdss(mapped, "b")
    assert f.qualifiers["pseudogene"] == ["unknown"]
    assert any("frameshift" in n for n in f.qualifiers["note"])
    assert f.id.startswith("MAP_")   # both fragments are gone


def test_c_internal_stop_pseudogene(mapped):
    (f,) = _cdss(mapped, "c")
    assert f.qualifiers["pseudogene"] == ["unknown"]
    assert any("internal stop codon" in n for n in f.qualifiers["note"])


def test_d_rescued(mapped):
    (f,) = _cdss(mapped, "d")
    assert f.id.startswith("MAP_") and _span(f) == (L, L + 600)
    assert f.primary_hit.id == "ARO:1" and f.qualifiers["translation"] == [_translation(GENE)]


def test_e_partial_cds_untouched(mapped):
    (f,) = _cdss(mapped, "e")
    assert f.id == "e" and _span(f) == (L, L + 400) and f.annotations["partial_flag"] == "01"
    assert f.primary_hit is None and "e" not in mapped[2]


def test_f_stop_extended(mapped):
    (f,) = _cdss(mapped, "f")
    assert _span(f) == (L, L + 630) and "pseudogene" not in f.qualifiers
    assert any("extended by 30 bp" in n for n in f.qualifiers["note"])


def test_g_gtg_start(mapped):
    (f,) = _cdss(mapped, "g")
    assert _span(f) == (L, L + 600) and f.qualifiers["translation"][0].startswith("M")
    assert "pseudogene" not in f.qualifiers


def test_hits_json(mapped):
    ids = {f.id for case in "abcdfg" for f in _cdss(mapped, case)}
    assert set(mapped[2]) == ids and len(ids) == 6
