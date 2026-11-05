from dfc.utils.mapping_util import (truncated_region, MappingHit, read_blast_table, is_full_length, indel_analysis,
                                    group_by_locus, is_valid_reference, resolve_region,
                                    MappedRegion, classify_cds)
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from dfc.models.bio_feature import ExtendedFeature


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


def test_group_by_locus_keeps_operon_neighbours_apart():
    a = _hit(ref_id="A")                                    # 1000..1300
    for start in (1300, 1290):                              # overlaps of 1 bp and 11 bp (adeA/adeB)
        b = _hit(ref_id="B", sstart=start, send=start + 299)
        assert [[h.ref_id for h in g] for g in group_by_locus([a, b])] == [["A"], ["B"]]


def test_group_by_locus_merges_references_of_different_lengths():
    a = _hit(ref_id="A", qlen=600, qend=600, sstart=1001, send=1600)
    b = _hit(ref_id="B", qlen=612, qend=612, sstart=1001, send=1612, bitscore=600)
    assert [[h.ref_id for h in g] for g in group_by_locus([a, b])] == [["B", "A"]]


def test_group_by_locus_is_transitive():
    x = _hit(ref_id="X")                                    # 1000..1300
    y = _hit(ref_id="Y", sstart=1281, send=1580)            # 1280..1580: overlaps X by 20 bp only
    z = _hit(ref_id="Z", qlen=15, qend=15, sstart=1286, send=1300, bitscore=10)   # 1285..1300: half in X and in Y
    assert [sorted(h.ref_id for h in g) for g in group_by_locus([x, y, z])] == [["X", "Y", "Z"]]


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


def test_tail_stop_before_read_through_reference_stop():
    # case A: stop one codon early, the reference's own stop position reads through, another stop 6 bp downstream
    gene = GENE[:-6] + "TAA" + "GCT" + "GCT" * 1 + "TAA"
    for strand in (1, -1):
        seq, hit = _region(gene, ref_len=len(GENE), strand=strand)
        if strand == -1:
            top = 100 + len(gene)
            hit = _hit(qlen=len(GENE), qend=len(GENE), sstart=top, send=top - len(GENE) + 1, slen=len(seq), btop=str(len(GENE)))
        r = resolve_region(seq, hit, 11, True)
        expected = (100, 163) if strand == 1 else (100 + len(gene) - 63, 100 + len(gene))
        assert r.kind == "valid" and (r.left0, r.right) == expected and r.notes == ["1 codon(s) shorter than the reference"]


def test_tail_stop_in_reference_without_stop():
    # case B: reference without a stop (63 bp); the genome has a stop at its last sense codon
    gene = GENE[:-6] + "TAA"
    seq, hit = _region(gene)
    r = resolve_region(seq, hit, 11, False)
    assert r.kind == "valid" and r.right == 100 + len(gene) and r.notes == ["1 codon(s) shorter than the reference"]
    assert r.translation == "M" + "A" * 19


def test_tail_stop_in_reference_without_stop_ignores_downstream_stop():
    # case C: stop in the second-to-last sense codon of a reference without stop, next stop 30 bp downstream
    gene = GENE[:-9] + "TAA" + "GCT" + "GCT" * 9 + "TAA"
    seq, hit = _region(gene, ref_len=len(GENE) - 3)
    r = resolve_region(seq, hit, 11, False)
    assert r.kind == "valid" and r.right == 100 + len(GENE) - 6 and r.notes == ["2 codon(s) shorter than the reference"]


REGION = MappedRegion("sequence1", 1000, 1300, 1, "valid", "M", [], set(), set())


def _cds(start, end, strand=1):
    return ExtendedFeature(location=FeatureLocation(start, end, strand=strand), type="CDS", id="MGA_1", seq_id="sequence1")


def test_classify_cds():
    assert classify_cds(_cds(1030, 1300), REGION, 30) == "same_gene"          # later start, same stop
    assert classify_cds(_cds(1000, 1300, -1), REGION, 30) == "overlapping"    # other strand, full overlap
    assert classify_cds(_cds(1100, 1200), REGION, 30) == "overlapping"        # fragment inside
    assert classify_cds(_cds(1260, 1600), REGION, 30) == "overlapping"        # 40 bp overlap
    assert classify_cds(_cds(1280, 1600), REGION, 30) == "neighbour"          # 20 bp overlap
    assert classify_cds(_cds(1300, 1600), REGION, 30) is None


def test_classify_cds_minus_strand_same_gene():
    region = MappedRegion("sequence1", 1000, 1300, -1, "valid", "M", [], set(), set())
    assert classify_cds(_cds(1000, 1270, -1), region, 30) == "same_gene"
    assert classify_cds(_cds(1000, 1270, 1), region, 30) == "overlapping"


# Truncated CDS: a frameshift or internal stop codon late in the gene leaves a long 5' ORF (e.g. E. coli yddK)
LONG = "ATG" + "GCT" * 16 + "AAA" * 3 + "TAA"   # 63 bp, 20 sense codons


def test_truncated_region_after_late_frameshift():
    gene = LONG[:51] + "T" + LONG[51:]          # 1-bp insertion: TAA right after codon 16
    seq, hit = _region(gene, ref_len=len(LONG), btop="51-T12")
    region = resolve_region(seq, hit, 11, True)
    assert region.kind == "frameshift"
    t = truncated_region(seq, region, hit, 11, True, 0.7)
    assert (t.kind, t.left0, t.right, t.strand) == ("valid", 100, 154, 1)
    assert t.translation == "M" + "A" * 16
    assert t.notes == ["85% of the reference length; frameshift at around 152"] and t.truncated == 85
    assert not t.indels and not t.stop_codons and region.truncated == 0


def test_truncated_region_minus_strand():
    gene = LONG[:51] + "T" + LONG[51:]
    seq, _ = _region(gene, strand=-1)
    hit = _hit(qlen=len(LONG), qend=len(LONG), sstart=100 + len(gene), send=101, slen=len(seq), btop="51-T12")
    t = truncated_region(seq, resolve_region(seq, hit, 11, True), hit, 11, True, 0.7)
    assert (t.kind, t.left0, t.right, t.strand) == ("valid", len(seq) - 154, len(seq) - 100, -1)
    assert t.translation == "M" + "A" * 16


def test_truncated_region_after_late_internal_stop():
    gene = LONG[:51] + "TAA" + LONG[54:]        # stop at codon 17 of 20
    seq, hit = _region(gene)
    region = resolve_region(seq, hit, 11, True)
    assert region.kind == "internal_stop"
    t = truncated_region(seq, region, hit, 11, True, 0.7)
    assert (t.kind, t.right) == ("valid", 154) and t.notes == ["85% of the reference length; internal stop codon at 152"]


def test_early_defect_stays_a_pseudogene():
    gene = LONG[:12] + "TAA" + LONG[15:]        # stop at codon 4 of 20
    seq, hit = _region(gene)
    region = resolve_region(seq, hit, 11, True)
    assert truncated_region(seq, region, hit, 11, True, 0.7) is region


def test_truncation_keeps_other_kinds():
    seq, hit = _region(LONG)
    region = resolve_region(seq, hit, 11, True)
    assert region.kind == "valid" and truncated_region(seq, region, hit, 11, True, 0.7) is region
