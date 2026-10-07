# coding: UTF8
"""Sequence names of a complete genome (dfc.genome.Genome.prepare_genome)."""
import os

from dfc.genome import Genome
from dfc.utils.config_util import load_config

APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _genome(tmp_path, n_seqs, **source_info):
    fasta = tmp_path / "genome.fna"
    fasta.write_text("".join(">contig{0} original\n{1}\n".format(i, "ACGT" * 300) for i in range(n_seqs)))
    config = load_config(APP_ROOT, os.path.join(APP_ROOT, "dfc", "default_config.py"))
    config.GENOME_FASTA, config.WORK_DIR = str(fasta), str(tmp_path / "OUT")
    config.GENOME_CONFIG["complete"] = True
    config.GENOME_SOURCE_INFORMATION.update(source_info)
    return Genome(config)


def test_complete_genome_without_seq_names_gets_default_name(tmp_path):
    # An empty name made the FASTA header '>', and nhmmer (Barrnap) failed.
    genome = _genome(tmp_path, 1)
    assert list(genome.seq_records) == ["sequence1"]
    assert genome.seq_records["sequence1"].name == "sequence1"


def test_complete_genome_without_seq_names_multiple_sequences(tmp_path):
    genome = _genome(tmp_path, 2, seq_types="c;p", seq_topologies="c;c")
    assert list(genome.seq_records) == ["sequence1", "sequence2"]


def test_complete_genome_with_seq_names(tmp_path):
    genome = _genome(tmp_path, 2, seq_names="chromosome;pA", seq_types="c;p", seq_topologies="c;c")
    assert list(genome.seq_records) == ["chromosome", "pA"]
