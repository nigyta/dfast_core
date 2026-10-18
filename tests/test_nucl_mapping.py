import logging
from collections import OrderedDict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import dfc.genome as genome_mod
from dfc.components.NuclMapping import NuclMapping
from dfc.functionalAnnotation import COMPONENTS

class FakeGenome:
    def __init__(self, records):
        self.seq_records = OrderedDict((r.id, r) for r in records)
        self.features = OrderedDict()
        self.set_feature_dictionary()
    sort_features = genome_mod.Genome.sort_features
    set_feature_dictionary = genome_mod.Genome.set_feature_dictionary

def _component(tmp_path, genome, **options):
    comp = NuclMapping.__new__(NuclMapping)          # skip __init__ (blastn version check, reference loading)
    comp.genome, comp.workDir, comp.pipeline_dir = genome, str(tmp_path), str(tmp_path)
    comp.logger = logging.getLogger("test")
    comp.db_name, comp.pident_cutoff, comp.overlap_tolerance, comp.transl_table = "CARD", 90, 30, 11
    comp.result_file = str(tmp_path / "mapping.out")
    comp.references, comp.ref_validity, comp.hits_json = {}, {}, {}
    for k, v in options.items(): setattr(comp, k, v)
    return comp

def test_registered():
    assert COMPONENTS["NuclMapping"] is NuclMapping

def test_select_hits_keeps_full_length_hits_on_known_contigs(tmp_path, caplog):
    genome = FakeGenome([SeqRecord(Seq("A" * 5000), id="sequence1")])
    comp = _component(tmp_path, genome)
    rows = ["ARO:1\tsequence1\t99.0\t300\t3\t0\t1\t300\t1001\t1300\t0.0\t500\t300\t5000\t300",
            "ARO:2\tsequence1\t99.0\t200\t2\t0\t50\t250\t2001\t2200\t0.0\t300\t300\t5000\t200",   # partial
            "ARO:3\tunknown\t99.0\t300\t3\t0\t1\t300\t1001\t1300\t0.0\t500\t300\t5000\t300"]      # unknown contig
    (tmp_path / "mapping.out").write_text("\n".join(rows) + "\n")
    groups = comp.select_hits()
    assert [[h.ref_id for h in g] for g in groups] == [["ARO:1"]]
    assert "unknown" in caplog.text

def test_blast_command(tmp_path):
    comp = _component(tmp_path, FakeGenome([]), database="/db/card/CARD", evalue_cutoff=1e-10, aligner_options="", CPU=4)
    cmd = comp.blast_command("/w/genome")
    assert cmd[:3] == ["blastn", "-task", "dc-megablast"] and "/db/card/CARD.nucl.fasta" in cmd
    assert cmd[cmd.index("-perc_identity") + 1] == "90" and cmd[cmd.index("-num_threads") + 1] == "4"
