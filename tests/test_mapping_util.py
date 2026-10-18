from dfc.utils.mapping_util import (MappingHit, read_blast_table, is_full_length, indel_analysis,
                                    group_by_locus, is_valid_reference, resolve_region)
from Bio.Seq import Seq


def _hit(**kw):
    base = dict(ref_id="ARO:1", seq_id="sequence1", pident=99.0, qstart=1, qend=300, sstart=1001, send=1300,
                evalue=0.0, bitscore=500.0, qlen=300, slen=5000, btop="300")
    base.update(kw)
    return MappingHit(**base)


def test_read_blast_table(tmp_path):
    f = tmp_path / "t.tsv"
    f.write_text("ARO:1\tsequence1\t99.0\t300\t3\t0\t1\t300\t1001\t1300\t0.0\t500\t300\t5000\t300\n")
    hits = read_blast_table(str(f))
    assert len(hits) == 1 and hits[0].ref_id == "ARO:1" and hits[0].btop == "300" and hits[0].slen == 5000


def test_read_blast_table_empty(tmp_path):
    f = tmp_path / "empty.tsv"
    f.write_text("")
    assert read_blast_table(str(f)) == []


def test_region_plus_strand_extrapolates_slack():
    h = _hit(qstart=2, qend=298, sstart=1002, send=1298)
    assert (h.strand, h.left0, h.right) == (1, 1000, 1300)


def test_region_minus_strand_extrapolates_slack():
    h = _hit(qstart=2, qend=298, sstart=1299, send=1003)
    assert (h.strand, h.left0, h.right) == (-1, 1000, 1300)


def test_full_length_slack_and_identity():
    assert is_full_length(_hit(qstart=4, qend=297), 90)
    assert not is_full_length(_hit(qstart=5), 90)
    assert not is_full_length(_hit(qend=296), 90)
    assert not is_full_length(_hit(pident=89.9), 90)


def test_full_length_rejects_region_outside_contig():
    assert not is_full_length(_hit(qstart=3, sstart=1, send=298), 90)        # would start before position 1
    assert not is_full_length(_hit(qend=298, sstart=4703, send=5000), 90)   # would end after the contig


def test_indel_analysis_plus_strand():
    h = _hit(btop="10A-20-C-C30")   # 1 deletion at 1011, then a 2-base insertion at 1031
    assert indel_analysis(h) == (1, [1011, 1031])


def test_indel_analysis_minus_strand():
    h = _hit(sstart=1300, send=1001, btop="10-A290")
    assert indel_analysis(h) == (1, [1290])


def test_indel_analysis_in_frame_gap_has_zero_shift():
    assert indel_analysis(_hit(btop="10A-A-A-287")) == (-3, [1011])


def test_group_by_locus_merges_overlapping_hits_and_sorts_by_score():
    a = _hit(ref_id="A", bitscore=400); b = _hit(ref_id="B", sstart=1010, send=1309, bitscore=500)
    c = _hit(ref_id="C", sstart=3001, send=3300)
    groups = group_by_locus([a, b, c])
    assert [[h.ref_id for h in g] for g in groups] == [["B", "A"], ["C"]]


def test_group_by_locus_separates_strands():
    a = _hit(ref_id="A"); b = _hit(ref_id="B", sstart=1300, send=1001)
    assert len(group_by_locus([a, b])) == 2


def test_group_by_locus_keeps_two_loci_of_one_reference():
    a = _hit(ref_id="A"); b = _hit(ref_id="A", sstart=4001, send=4300)
    assert [[h.sstart for h in g] for g in group_by_locus([a, b])] == [[1001], [4001]]


GENE = "ATG" + "GCT" * 20 + "TAA"


def test_reference_validity():
    assert is_valid_reference(GENE, 11) == (True, True)
    assert is_valid_reference(GENE[:-3], 11) == (True, False)          # no terminal stop
    assert is_valid_reference("CCC" + GENE[3:], 11) == (False, True)   # no start codon
    assert is_valid_reference(GENE[:30] + "TAA" + GENE[30:], 11)[0] is False   # internal stop
    assert is_valid_reference(GENE + "A", 11)[0] is False              # not a multiple of 3
    assert is_valid_reference(GENE.replace("GCT", "GNT", 1), 11)[0] is False   # ambiguous base


def _region(gene, ref_len=None, btop=None, strand=1):
    """Genome = 100 C + gene + 100 C (C flanks: no stop codon in any frame); a hit mapping a reference of
    ref_len bases (default: the gene length) over the gene, on the given strand."""
    seq = Seq("C" * 100 + (gene if strand == 1 else str(Seq(gene).reverse_complement())) + "C" * 100)
    n = ref_len or len(gene)
    s, e = (101, 100 + n) if strand == 1 else (100 + n, 101)
    return seq, _hit(qlen=n, qend=n, sstart=s, send=e, slen=len(seq), btop=btop or str(n))


def test_valid_region():
    seq, hit = _region(GENE)
    r = resolve_region(seq, hit, 11, True)
    assert (r.kind, r.left0, r.right, r.strand) == ("valid", 100, 100 + len(GENE), 1)
    assert r.translation == "M" + "A" * 20 and r.notes == []


def test_valid_region_minus_strand():
    seq, hit = _region(GENE, strand=-1)
    r = resolve_region(seq, hit, 11, True)
    assert (r.kind, r.left0, r.right, r.strand, r.translation) == ("valid", 100, 100 + len(GENE), -1, "M" + "A" * 20)


def test_frameshift_from_btop():
    seq, hit = _region(GENE, btop="10-A" + str(len(GENE) - 10))
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "frameshift" and r.indels == {111} and r.translation == "" and "frameshift at around 111" in r.notes


def test_internal_stop():
    gene = GENE[:30] + "TAA" + GENE[33:]
    seq, hit = _region(gene)
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "internal_stop" and r.stop_codons == {(130, 133, 1)} and r.notes == ["internal stop codon at 131"]


def test_internal_stop_minus_strand():
    gene = GENE[:30] + "TAA" + GENE[33:]
    seq, hit = _region(gene, strand=-1)
    r = resolve_region(seq, hit, 11, True)
    # codon 10 of the gene occupies genome [100 + len - 33, 100 + len - 30)
    start0 = 100 + len(gene) - 33
    assert r.kind == "internal_stop" and r.stop_codons == {(start0, start0 + 3, -1)} and r.notes == [f"internal stop codon at {start0 + 1}"]


def test_stop_in_last_two_codons_shortens():
    gene = GENE[:-6] + "TAA" + GENE[-3:]     # stop one sense codon early
    seq, hit = _region(gene)
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "valid" and r.right == 100 + len(gene) - 3 and r.notes == ["1 codon(s) shorter than the reference"]


def test_stop_in_last_two_codons_shortens_minus_strand():
    gene = GENE[:-6] + "TAA" + GENE[-3:]
    seq, hit = _region(gene, strand=-1)
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "valid" and (r.left0, r.right) == (103, 100 + len(gene)) and r.translation == "M" + "A" * 19


def test_lost_stop_extended_within_limit():
    gene = GENE[:-3] + "GCT" + "GCT" * 3 + "TAG"   # stop lost; next in-frame stop 12 bp later
    seq, hit = _region(gene, ref_len=len(GENE))
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "valid" and r.right == 100 + len(gene) and r.notes == ["stop codon of the reference is not conserved; extended by 12 bp"]


def test_lost_stop_extended_minus_strand():
    gene = GENE[:-3] + "GCT" + "GCT" * 3 + "TAG"
    seq, _ = _region(gene, strand=-1)
    top = 100 + len(gene)   # the gene starts at the highest genome coordinate; the hit covers its first len(GENE) bases
    hit = _hit(qlen=len(GENE), qend=len(GENE), sstart=top, send=top - len(GENE) + 1, slen=len(seq), btop=str(len(GENE)))
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "valid" and (r.left0, r.right) == (100, 100 + len(gene)) and r.notes == ["stop codon of the reference is not conserved; extended by 12 bp"]


def test_lost_stop_extended_silently_without_reference_stop():
    gene = GENE[:-3] + "TAG"
    seq, hit = _region(gene, ref_len=len(GENE) - 3)
    r = resolve_region(seq, hit, 11, False)
    assert r.kind == "valid" and r.right == 100 + len(gene) and r.notes == []


def test_lost_stop_beyond_limit_is_pseudogene():
    seq, hit = _region(GENE[:-3] + "GCT")             # C flank: no stop within 100 bp
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "stop_lost" and r.notes == ["stop codon lost"]


def test_reference_without_stop_and_no_stop_found_is_ignored():
    seq, hit = _region(GENE[:-3])
    assert resolve_region(seq, hit, 11, False) is None


def test_no_start_codon():
    seq, hit = _region("CCC" + GENE[3:])
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "no_start" and r.notes == ["no start codon"] and r.translation.endswith("A" * 20)


def test_alternative_start_codon_translates_to_m():
    seq, hit = _region("GTG" + GENE[3:])
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "valid" and r.translation == "M" + "A" * 20
