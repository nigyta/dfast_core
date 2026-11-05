# coding: UTF8
"""accessory_summary.tsv (dfc.utils.summarize_accessory) and the rows of the IntegronFinder tool."""
import json
import os
import shutil
from collections import OrderedDict

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

from dfc.models.bio_feature import ExtendedFeature
from dfc.tools.integron_finder import IntegronFinder
from dfc.utils.summarize_accessory import summarize_accessory, HEADER

INTEGRONS = os.path.join(os.path.dirname(__file__), "data", "integron_finder_AYE.integrons")


class _Genome:
    def __init__(self, records):
        self.seq_records = OrderedDict((r.id, r) for r in records)
        self.features = {f.id: f for r in records for f in r.features}


def _cds(fid, start, end, strand, seq_id, locus_tag):
    return ExtendedFeature(location=FeatureLocation(start, end, strand=strand), type="CDS", id=fid, seq_id=seq_id,
                           qualifiers={"locus_tag": [locus_tag]})


def _genome():
    rec1 = SeqRecord(Seq("A" * 5000), id="sequence1",
                     features=[_cds("MGA_1", 99, 399, -1, "sequence1", "LOCUS_010"),
                               _cds("MGA_2", 999, 1299, 1, "sequence1", "LOCUS_020"),
                               _cds("MGA_3", 1299, 1599, 1, "sequence1", "LOCUS_030")])
    rec2 = SeqRecord(Seq("A" * 3000), id="sequence2", features=[_cds("MGA_4", 9, 309, 1, "sequence2", "LOCUS_040")])
    return _Genome([rec1, rec2])


def _nucl_hit(db_name, gene, product, note):
    return {"db_name": db_name, "identity": 99.04, "model": {"gene": gene, "product": product, "note": note}}


def test_integron_rows(tmp_path):
    tool = IntegronFinder.__new__(IntegronFinder)  # skip the version check of integron_finder
    tool.output_directory = str(tmp_path)
    os.makedirs(tmp_path / "Results_Integron_Finder_genome")
    shutil.copy(INTEGRONS, tmp_path / "Results_Integron_Finder_genome" / "genome.integrons")
    rows = {r["location"]: r for r in tool.getSummaryRows()}
    assert rows["3620981..3625513"] == {"module": "integron", "sequence": "sequence1", "location": "3620981..3625513",
                                        "gene/system": "complete integron", "note": "attC sites: 4; integrase: intI"}
    assert rows["3677546..3678559"]["gene/system"].startswith("In0")


def test_summarize_accessory(tmp_path):
    genome = _genome()
    for name, hits in [("NuclMapping", {"MGA_4": _nucl_hit("CARD", "OXA-48", "OXA-48-like beta-lactamase", "rescued; similar to CARD:OXA-48")}),
                       ("NuclMapping_1", {"MGA_1": _nucl_hit("VFDB", "espL1", "effector espL1", "similar to espL1"),
                                          "MGA_9": _nucl_hit("VFDB", "removed", "removed", "removed")}),
                       ("NuclSearch", {"MGA_2": _nucl_hit("CARD", "aadA", "ANT(3'')", "similar to aadA")})]:
        os.makedirs(tmp_path / name)
        json.dump(hits, open(tmp_path / name / "nucl_hits.json", "w"))
    os.makedirs(tmp_path / "DefenseFinder")
    json.dump({"R0001_Gabija_1": {"type": "Gabija", "subtype": "Gabija", "activity": "Defense",
                                  "features": ["MGA_2", "MGA_3"], "genes": ["GajA", "GajB_2"]}},
              open(tmp_path / "DefenseFinder" / "systems.json", "w"))
    contig_rows = [{"module": "integron", "sequence": "sequence1", "location": "2001..3000",
                    "gene/system": "complete integron", "note": "attC sites: 4; integrase: intI"},
                   {"module": "plasmid", "sequence": "sequence2", "location": "101..400", "gene/system": "IncL",
                    "identity": "100.0", "note": "similar to JN626286"}]
    out = tmp_path / "accessory_summary.tsv"
    summarize_accessory(genome, str(tmp_path), contig_rows, str(out))

    lines = open(out).read().splitlines()
    assert lines[0] == "#" + "\t".join(HEADER)
    rows = [line.split("\t") for line in lines[1:]]
    # sorted by module, then sequence, then position; a CDS removed after the annotation is skipped
    assert [(r[0], r[1], r[2]) for r in rows] == [
        ("amr", "sequence1", "1000..1299"), ("amr", "sequence2", "10..309"), ("vf", "sequence1", "complement(100..399)"),
        ("plasmid", "sequence2", "101..400"), ("integron", "sequence1", "2001..3000"), ("defense", "sequence1", "1000..1599")]
    assert rows[1] == ["amr", "sequence2", "10..309", "LOCUS_040", "OXA-48", "99.0",
                       "OXA-48-like beta-lactamase; rescued; similar to CARD:OXA-48"]
    assert rows[4][3:6] == ["", "complete integron", ""]
    assert rows[5][3:] == ["LOCUS_020,LOCUS_030", "Gabija defense system", "", "genes: GajA, GajB_2"]


def test_no_rows_no_file(tmp_path):
    out = tmp_path / "accessory_summary.tsv"
    summarize_accessory(_genome(), str(tmp_path), [], str(out))
    assert not out.exists()
