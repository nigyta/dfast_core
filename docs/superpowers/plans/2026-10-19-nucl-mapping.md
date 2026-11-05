# NuclMapping Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Annotate CARD/VFDB genes by mapping the reference nucleotide sequences directly onto the genome with blastn, so that unpredicted genes are rescued, start codons are corrected and frameshifted/interrupted genes become `/pseudogene` CDSs.

**Architecture:** A new functional annotation component `NuclMapping` (placed before the NuclSearch entries) builds a BLAST database of the genome, maps `<database>.nucl.fasta` onto it, keeps full-length high-identity hits, and mutates the genome's CDS features. The pure logic (hit table, full-length test, btop indel analysis, locus grouping, reference validity, region resolution, CDS classification) lives in `dfc/utils/mapping_util.py` and is unit-tested without blastn; the component does I/O and feature mutation. Enabled by `--amr_map` (`enable_amr_map`).

**Tech Stack:** Python 3.10+, Biopython (`Seq`, `FeatureLocation`, `CodonTable`), BLAST+ `makeblastdb`/`blastn` on PATH, pytest.

**Spec:** `docs/superpowers/specs/2026-10-19-nucl-mapping-design.md`

## Global Constraints

- Python >= 3.10, Biopython only; `blastn` and `makeblastdb` (BLAST+ >= 2.13) on PATH.
- Everything in tracked files (code, comments, tests, docs, commit messages) is English.
- No result changes unless `--amr_map` or the config entries enable the component.
- `pident_cutoff` default 90; `overlap_tolerance` default 30 bp; end slack 3 bp at each reference end; stop-codon extension limit 100 bp; stop within the last 2 sense codons shortens the CDS.
- `example/dfc_config/dfast_web_config.py` is kept in sync with `dfc/default_config.py`.
- Component pattern: `BaseAnnotationComponent` subclass registered in `dfc/functionalAnnotation.py:COMPONENTS`; enable helper in `dfc/utils/config_util.py`; CLI flag in `dfast`.
- Run tests with `python -m pytest -q tests` (env `dfast-dev`, with a `mefinder` on PATH or a fake one; see CLAUDE.md).
- Commit after each task; branch `feature/amr-map`.

## Review Focus

1. A reference mapping to the minus strand near position 1 of a contig, or to the plus strand near the contig end, where the extrapolated region would leave the contig: must be treated as not full-length, never produce a negative or out-of-range location. (Task 1: `test_full_length_rejects_region_outside_contig`.)
2. A reference that maps full-length to two loci (two gene copies): both loci must be annotated, each with its own feature. (Task 2: `test_group_by_locus_keeps_two_loci_of_one_reference`.)
3. The genome DB sequence ids contain the renamed ids (`sequence1`); a blastn table whose `saccver` is not a known `seq_id` (e.g. `-parse_seqids` prefixing) must not crash: unknown contigs are skipped with a warning. (Task 5: `test_apply_skips_unknown_contig`.)
4. A predicted CDS overlapping two mapped regions: it must be removed once and never referenced again (no KeyError on the second removal). (Task 6: `test_cds_overlapping_two_regions_removed_once`.)
5. `--amr` and `--amr_map` given together: NuclSearch must be disabled with a warning and every other `--amr` effect kept. (Task 8: `test_amr_map_with_amr_disables_nuclsearch`.)

---

### Task 1: Hit table, full-length test and indel analysis

**Files:**
- Create: `dfc/utils/mapping_util.py`
- Test: `tests/test_mapping_util.py`

**Interfaces:**
- Produces (used by every later task):
  - `BLAST_OUTFMT = "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen btop"`
  - `END_SLACK = 3`, `STOP_EXTENSION = 100`, `TAIL_STOP_CODONS = 2`
  - `@dataclass class MappingHit` with fields `ref_id: str, seq_id: str, pident: float, qstart: int, qend: int, sstart: int, send: int, evalue: float, bitscore: float, qlen: int, slen: int, btop: str` and properties `strand -> int` (1 if `sstart <= send` else -1), `left0 -> int` and `right -> int`: the 0-based start and exclusive end of the region extrapolated to the full reference. Plus strand: `left0 = sstart - 1 - (qstart - 1)`, `right = send + (qlen - qend)`. Minus strand: `left0 = send - 1 - (qlen - qend)`, `right = sstart + (qstart - 1)`.
  - `read_blast_table(file_name: str) -> list[MappingHit]` (one row per line of `-outfmt BLAST_OUTFMT`; empty file -> `[]`).
  - `is_full_length(hit: MappingHit, pident_cutoff: float) -> bool`: `qstart <= 1 + END_SLACK`, `qend >= qlen - END_SLACK`, `pident >= pident_cutoff`, and `0 <= left0 < right <= slen`.
  - `indel_analysis(hit: MappingHit) -> tuple[int, list[int]]`: `(shift, positions)` where `shift` is (bases inserted in the genome − bases deleted from the genome relative to the reference) and `positions` are the 1-based genome coordinates where each run of gaps starts, both derived from `btop`. A btop item is a run of digits (matches), a pair of letters (mismatch) or a pair with `-` (`A-` = base in the query/reference missing in the genome = deletion; `-A` = base in the genome missing in the reference = insertion). The genome coordinate advances with the subject: `sstart + offset` on the plus strand, `sstart - offset` on the minus strand.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_mapping_util.py
from dfc.utils.mapping_util import MappingHit, read_blast_table, is_full_length, indel_analysis

def _hit(**kw):
    base = dict(ref_id="ARO:1", seq_id="sequence1", pident=99.0, qstart=1, qend=300, sstart=1001, send=1300,
                evalue=0.0, bitscore=500.0, qlen=300, slen=5000, btop="300")
    base.update(kw); return MappingHit(**base)

def test_read_blast_table(tmp_path):
    f = tmp_path / "t.tsv"
    f.write_text("ARO:1\tsequence1\t99.0\t300\t3\t0\t1\t300\t1001\t1300\t0.0\t500\t300\t5000\t300\n")
    hits = read_blast_table(str(f))
    assert len(hits) == 1 and hits[0].ref_id == "ARO:1" and hits[0].btop == "300" and hits[0].slen == 5000

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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest -q tests/test_mapping_util.py`
Expected: FAIL (ImportError: no module `dfc.utils.mapping_util`)

- [ ] **Step 3: Implement `dfc/utils/mapping_util.py`: the constants, `MappingHit`, `read_blast_table`, `is_full_length`, `indel_analysis`**

`indel_analysis` walks the btop with a regular expression for `\d+|[A-Z-]{2}`; consecutive gap items of the same kind form one run (one position); `shift` adds the run length for insertions (`-X`) and subtracts it for deletions (`X-`).

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest -q tests/test_mapping_util.py`
Expected: 9 passed

- [ ] **Step 5: Commit**

```bash
git add dfc/utils/mapping_util.py tests/test_mapping_util.py
git commit -m "feat: parse blastn mappings of nucleotide references onto the genome"
```

---

### Task 2: Locus grouping and reference validity

**Files:**
- Modify: `dfc/utils/mapping_util.py`
- Test: `tests/test_mapping_util.py`

**Interfaces:**
- Produces:
  - `group_by_locus(hits: list[MappingHit]) -> list[list[MappingHit]]`: groups hits whose extrapolated regions overlap on the same `seq_id` and `strand`; each group sorted by `(bitscore, pident)` descending; groups ordered by `(seq_id, left0)`. Input is already full-length-filtered.
  - `is_valid_reference(nucl_seq: str, transl_table: int) -> tuple[bool, bool]`: `(valid, has_terminal_stop)`. Valid: length multiple of 3, only `ACGT`, first codon in `CodonTable.unambiguous_dna_by_id[transl_table].start_codons`, no internal stop in the translation. `has_terminal_stop`: the last codon is a stop codon of the table. A sequence that is valid but lacks the terminal stop is `(True, False)`.

- [ ] **Step 1: Write the failing tests**

```python
from dfc.utils.mapping_util import group_by_locus, is_valid_reference

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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest -q tests/test_mapping_util.py -k "locus or validity"`
Expected: FAIL (ImportError)

- [ ] **Step 3: Implement `group_by_locus` and `is_valid_reference`**

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest -q tests/test_mapping_util.py`
Expected: 13 passed

- [ ] **Step 5: Commit**

```bash
git add dfc/utils/mapping_util.py tests/test_mapping_util.py
git commit -m "feat: group reference mappings by locus and check reference validity"
```

---

### Task 3: Region resolution (valid CDS, defects, 3' tolerances)

**Files:**
- Modify: `dfc/utils/mapping_util.py`
- Test: `tests/test_mapping_util.py`

**Interfaces:**
- Produces:
  - `@dataclass class MappedRegion`: `seq_id: str, left0: int, right: int, strand: int, kind: str, translation: str, notes: list[str], stop_codons: set[tuple[int, int, int]], indels: set[int]`. `kind` is one of `"valid"`, `"no_start"`, `"frameshift"`, `"internal_stop"`, `"stop_lost"`. `translation` is the protein without the stop (`""` for frameshift/internal_stop/stop_lost). `stop_codons` uses the `PseudoGene` format `(start0, end, strand)` with genome coordinates; `indels` are 1-based genome positions from `indel_analysis`.
  - `resolve_region(record_seq, hit: MappingHit, transl_table: int, ref_has_stop: bool) -> MappedRegion | None`, applying the spec rules in this order:
    1. `shift, positions = indel_analysis(hit)`; `shift % 3 != 0` → `kind="frameshift"`, `notes=["frameshift at around <p>", ...]`, region = extrapolated hit region.
    2. Take the region sequence (reverse-complemented on the minus strand). If it does not end with a stop codon: look for the next in-frame stop within `STOP_EXTENSION` bp beyond it (within the contig). Found and `ref_has_stop` → extend, `notes=["stop codon of the reference is not conserved; extended by N bp"]`. Found and not `ref_has_stop` → extend silently. Not found and `ref_has_stop` → `kind="stop_lost"`, `notes=["stop codon lost"]`. Not found and not `ref_has_stop` → return `None`.
    3. Translate (without the terminal stop). Internal stops: if all of them lie within the last `TAIL_STOP_CODONS` sense codons → shorten the region to the first of them, `notes=["N codon(s) shorter than the reference"]`; otherwise `kind="internal_stop"`, `notes=["internal stop codon at <1-based genome position of the first base of the codon>", ...]`, `stop_codons` filled.
    4. First codon not a start codon of the table → `kind="no_start"`, `notes=["no start codon"]` (translation kept).
    5. Otherwise `kind="valid"`; the translation starts with `M` (as `StructuralAnnotation._get_translation` does for alternative start codons).

- [ ] **Step 1: Write the failing tests**

```python
from Bio.Seq import Seq
from dfc.utils.mapping_util import resolve_region

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

def test_stop_in_last_two_codons_shortens():
    gene = GENE[:-6] + "TAA" + GENE[-3:]     # stop one sense codon early
    seq, hit = _region(gene)
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "valid" and r.right == 100 + len(gene) - 3 and r.notes == ["1 codon(s) shorter than the reference"]

def test_lost_stop_extended_within_limit():
    gene = GENE[:-3] + "GCT" + "GCT" * 3 + "TAG"   # stop lost; next in-frame stop 12 bp later
    seq, hit = _region(gene, ref_len=len(GENE))
    r = resolve_region(seq, hit, 11, True)
    assert r.kind == "valid" and r.right == 100 + len(gene) and r.notes == ["stop codon of the reference is not conserved; extended by 12 bp"]

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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest -q tests/test_mapping_util.py -k "region or stop or start or frameshift"`
Expected: FAIL (ImportError)

- [ ] **Step 3: Implement `MappedRegion` and `resolve_region`**

Use `Bio.Data.CodonTable.unambiguous_dna_by_id[transl_table]` for start and stop codons and `Seq.translate(table=transl_table)` for the translation; shortening/extension adjust `right` on the plus strand and `left0` on the minus strand.

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest -q tests/test_mapping_util.py`
Expected: 22 passed

- [ ] **Step 5: Commit**

```bash
git add dfc/utils/mapping_util.py tests/test_mapping_util.py
git commit -m "feat: classify a mapped region as a CDS or a pseudogene"
```

---

### Task 4: Classification of predicted CDSs

**Files:**
- Modify: `dfc/utils/mapping_util.py`
- Test: `tests/test_mapping_util.py`

**Interfaces:**
- Produces: `classify_cds(cds, region: MappedRegion, tolerance: int) -> str | None` where `cds` is an `ExtendedFeature` on the same contig. Returns `"same_gene"` (same strand and same stop coordinate: `int(cds.location.end) == region.right` on the plus strand, `int(cds.location.start) == region.left0` on the minus strand), `"overlapping"` (overlap length `min(ends) - max(starts) > tolerance`, either strand), `"neighbour"` (overlap of 1..tolerance bases), or `None` (no overlap).

- [ ] **Step 1: Write the failing tests**

```python
from Bio.SeqFeature import FeatureLocation
from dfc.models.bio_feature import ExtendedFeature
from dfc.utils.mapping_util import MappedRegion, classify_cds

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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest -q tests/test_mapping_util.py -k classify`
Expected: FAIL (ImportError)

- [ ] **Step 3: Implement `classify_cds`**

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest -q tests/test_mapping_util.py`
Expected: 24 passed

- [ ] **Step 5: Commit**

```bash
git add dfc/utils/mapping_util.py tests/test_mapping_util.py
git commit -m "feat: classify predicted CDSs against a mapped reference"
```

---

### Task 5: Component: references, BLAST commands, hit selection

**Files:**
- Create: `dfc/components/NuclMapping.py`
- Modify: `dfc/functionalAnnotation.py` (import and `COMPONENTS["NuclMapping"]`)
- Test: `tests/test_nucl_mapping.py`

**Interfaces:**
- Consumes: Task 1–2 functions; `dfc.components.NuclSearch.models` (`{"CARD": CARD, "VFDB": VFDB, ...}`), `NucRefBase.load_from_tsv_file`, `Blastn` through `BaseAnnotationComponent` (`options["aligner"] = "blastn"` gives `self.aligner`, whose `format_db_command(fasta, db_name)` builds `makeblastdb`).
- Produces:
  - `class NuclMapping(BaseAnnotationComponent)` with `instances = 0`; `__init__(self, genome, options, workDir, CPU)` reads `database`, `db_name`, `pident_cutoff` (90), `evalue_cutoff` (1e-10), `overlap_tolerance` (30), `transl_table` (11), `aligner_options` (""), sets `options["aligner"] = "blastn"` before calling `super().__init__`, logs the DB version like `NuclSearch.get_db_version`, and calls `self.load_references()`.
  - `load_references(self) -> None`: `self.references: dict[str, model]` from `<database>.nucl.ref`; `self.ref_validity: dict[str, tuple[bool, bool]]` from `is_valid_reference(model.nucl_seq, transl_table)`.
  - `genome_fasta(self) -> str`: `os.path.join(<WORK_DIR>, "input", "genome.fna")` (the component's `workDir` is `<WORK_DIR>/NuclMapping[_N]`; keep the pipeline work dir as `self.pipeline_dir`).
  - `create_db(self) -> str`: runs `makeblastdb` on the genome FASTA into `<workDir>/genome` and returns the prefix.
  - `blast_command(self, db_prefix: str) -> list[str]`: `["blastn", "-task", "dc-megablast", "-query", <database>.nucl.fasta, "-db", db_prefix, "-out", self.result_file, "-outfmt", "'" + BLAST_OUTFMT + "'", "-evalue", str(evalue_cutoff), "-perc_identity", str(pident_cutoff), "-max_target_seqs", "1000", "-max_hsps", "50", "-num_threads", str(self.CPU)]` plus `aligner_options.split()`; `self.result_file = <workDir>/mapping.out`. (`dc-megablast` is the starting choice; Task 10 measures it.)
  - `select_hits(self) -> list[list[MappingHit]]`: reads `self.result_file`, keeps `is_full_length` hits whose `seq_id` is in `self.genome.seq_records` (others logged with a warning and dropped), returns `group_by_locus`.
  - `run(self)`: `create_db`, run the blast command through `self.executeCommand(cmd, shell=True)`, `self.apply(self.select_hits())` (Task 6), `self.write_hits()` (Task 7).

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_nucl_mapping.py
import os, shutil, logging
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest -q tests/test_nucl_mapping.py`
Expected: FAIL (ImportError)

- [ ] **Step 3: Implement `dfc/components/NuclMapping.py` (`__init__`, `load_references`, `genome_fasta`, `create_db`, `blast_command`, `select_hits`, a `run` that calls `apply` and `write_hits` to be added in Tasks 6–7) and register it in `dfc/functionalAnnotation.py`**

`apply` and `write_hits` may be empty methods in this task.

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest -q tests/test_nucl_mapping.py`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add dfc/components/NuclMapping.py dfc/functionalAnnotation.py tests/test_nucl_mapping.py
git commit -m "feat: NuclMapping component maps CARD/VFDB references onto the genome"
```

---

### Task 6: Feature actions (rescue, annotate, start correction, pseudogene, removal)

**Files:**
- Modify: `dfc/components/NuclMapping.py`
- Test: `tests/test_nucl_mapping.py`

**Interfaces:**
- Consumes: Task 3 `resolve_region`, Task 4 `classify_cds`, `NuclHit(id_, model, db_name, e_value, score, identity, q_cov, s_cov, flag)`, `PseudoGene(ref_id, stop_codon, indel)` (`dfc/models/hit.py`).
- Produces (module-level, testable with `FakeGenome`):
  - `new_cds(region: MappedRegion, feature_id: str, transl_table: int) -> ExtendedFeature`: type `CDS`, `FeatureLocation(left0, right, strand)`, qualifiers `translation=[region.translation]`, `transl_table=[str(transl_table)]`, `codon_start=[1]`, `note=list(region.notes)`; `annotations={"partial_flag": "00", "nucl_mapping": "rescued"}`.
  - `new_pseudogene(region, feature_id, transl_table) -> ExtendedFeature`: type `CDS`, qualifiers `pseudogene=["unknown"]`, `transl_table`, `codon_start=[1]`, `note=list(region.notes)`, no `translation`; `annotations={"partial_flag": "00", "nucl_mapping": "pseudogene"}`.
  - `relocate_cds(cds: ExtendedFeature, region: MappedRegion, ref_label: str) -> None`: sets the location to the region, `translation`, `codon_start=[1]`, appends the note `"start codon moved N bp upstream|downstream by mapping of <ref_label>"` (N = difference of the start coordinates; upstream = the CDS became longer), `annotations["nucl_mapping"] = "start corrected"`, keeps the id. `ref_label` is `f"{db_name}:{model.to_dict()['gene']}"` (e.g. `CARD:OXA-48`).
  - `NuclMapping.apply(self, groups: list[list[MappingHit]]) -> None`: for each group, in this order:
    1. `primary = group[0]`; `model = self.references[primary.ref_id]`; `valid, has_stop = self.ref_validity[primary.ref_id]`; `record = self.genome.seq_records[primary.seq_id]`.
    2. If a feature with `annotations.get("nucl_mapping")` overlaps the extrapolated region on the same strand: same location → attach the hit as a secondary hit and continue; otherwise log at DEBUG and continue (earlier instance wins).
    3. `region = resolve_region(record.seq, primary, self.transl_table, has_stop)`; `None` → continue. For an invalid reference, `region` is only used for `classify_cds`: a `same_gene` CDS is annotated (hit primary if the CDS has none, else secondary), nothing else happens; continue.
    4. `classes = {cds.id: classify_cds(cds, region, self.overlap_tolerance)}` over the record's CDS features that have no `nucl_mapping` annotation; `same = first same_gene`.
    5. Action by `region.kind` (table in the spec): valid → rescue / annotate / start correction; no_start → annotate if `same` else pseudogene; frameshift / internal_stop / stop_lost → pseudogene (removing `same`). A pseudogene feature also gets `PseudoGene(primary.ref_id, region.stop_codons, region.indels)` appended to `secondary_hits`.
    6. Remove every `overlapping` CDS from `record.features` (a CDS already removed by an earlier group is skipped).
    7. The hit `NuclHit(primary.ref_id, model, self.db_name, primary.evalue, primary.bitscore, primary.pident, 100.0, 100.0, "")` becomes `feature.primary_hit` (replacing any) except for the annotate action, where it follows NuclSearch (primary if none, else secondary). Other hits of the group become secondary `NuclHit`s with flag `"alternative hit"`.
    8. `self.hits_json[feature.id] = d` where `d = hit.to_dict()` and, unless the action is annotate, `d["model"]["note"] = action + "; " + d["model"]["note"]` (actions: `"rescued"`, `"start corrected"`, `"pseudogene: frameshift"`, `"pseudogene: internal stop codon"`, `"pseudogene: stop codon lost"`, `"pseudogene: no start codon"`).
    9. After all groups: `self.genome.sort_features()`, `self.genome.set_feature_dictionary()`, INFO log of the counts (rescued, start corrected, pseudogene, annotated, removed CDSs).
  - Rescued/pseudogene feature ids: `"MAP_{n}"` with a class-level counter shared by the instances.

- [ ] **Step 1: Write the failing tests**

Use `GENE` from Task 2 embedded in a `FakeGenome` record with C flanks (`"C" * 100 + gene + "C" * 100`), a reference model built with `CARD("ARO:1", "OXA-48", "OXA beta-lactamase", "P1", "", "N1", gene, "E. coli", "")` registered in `comp.references`, `comp.ref_validity = {"ARO:1": (True, True)}`, and groups built from `_hit(...)` of Task 1 with the matching coordinates (`sstart=101`, `send=100+len(gene)`, `qlen=len(gene)`, `qend=len(gene)`, `slen=len(record)`).

```python
from Bio.SeqFeature import FeatureLocation
from dfc.models.bio_feature import ExtendedFeature
from dfc.models.card import CARD
from dfc.models.hit import NuclHit, PseudoGene
from dfc.components.NuclMapping import new_cds, new_pseudogene, relocate_cds
from tests.test_mapping_util import _hit, GENE

def _setup(tmp_path, gene, cds=(), validity=(True, True), btop=None, ref_len=None):
    """One contig: 100 C + gene + 100 C, with predicted CDSs given as (id, start0, end, strand)."""
    rec = SeqRecord(Seq("C" * 100 + gene + "C" * 100), id="sequence1")
    rec.features = [ExtendedFeature(location=FeatureLocation(s, e, strand=st), type="CDS", id=i, seq_id="sequence1",
                                    qualifiers={"translation": ["M"], "transl_table": ["11"], "codon_start": [1]},
                                    annotations={"partial_flag": "00"}) for i, s, e, st in cds]
    comp = _component(tmp_path, FakeGenome([rec]))
    comp.references = {"ARO:1": CARD("ARO:1", "OXA-48", "OXA beta-lactamase", "P1", "", "N1", gene, "E. coli", "")}
    comp.ref_validity = {"ARO:1": validity}
    n = ref_len or len(gene)
    hit = _hit(qlen=n, qend=n, sstart=101, send=100 + n, slen=len(rec), btop=btop or str(n))
    return comp, rec, hit

def _cds_features(rec):
    return {f.id: f for f in rec.features if f.type == "CDS"}

def test_rescue_adds_cds_with_primary_hit(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    comp.apply([[hit]])
    f = _cds_features(rec)["MAP_1"]
    assert (int(f.location.start), int(f.location.end), f.location.strand) == (100, 100 + len(GENE), 1)
    assert f.qualifiers["translation"] == ["M" + "A" * 20] and f.primary_hit.id == "ARO:1"
    assert comp.hits_json["MAP_1"]["model"]["note"].startswith("rescued; ")
    assert comp.genome.features["MAP_1"] is f

def test_annotate_keeps_identical_cds(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 100, 100 + len(GENE), 1)])
    comp.apply([[hit]])
    assert list(_cds_features(rec)) == ["MGA_1"]
    f = rec.features[0]
    assert f.primary_hit.id == "ARO:1" and "nucl_mapping" not in f.annotations
    assert not comp.hits_json["MGA_1"]["model"]["note"].startswith(("rescued", "start corrected", "pseudogene"))

def test_start_correction_relocates_same_gene(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 130, 100 + len(GENE), 1)])
    comp.apply([[hit]])
    f = _cds_features(rec)["MGA_1"]
    assert int(f.location.start) == 100 and f.qualifiers["translation"] == ["M" + "A" * 20]
    assert "start codon moved 30 bp upstream by mapping of CARD:OXA-48" in f.qualifiers["note"]
    assert f.annotations["nucl_mapping"] == "start corrected" and f.primary_hit.id == "ARO:1"
    assert comp.hits_json["MGA_1"]["model"]["note"].startswith("start corrected; ")

def test_frameshift_makes_pseudogene_and_removes_fragments(tmp_path):
    gene = GENE[:10] + GENE[11:]   # 1-bp deletion at position 111 of the genome
    comp, rec, hit = _setup(tmp_path, gene, ref_len=len(GENE), btop="10A-" + str(len(GENE) - 11),
                            cds=[("MGA_1", 100, 140, 1), ("MGA_2", 140, 100 + len(gene), 1)])
    comp.apply([[hit]])
    assert list(_cds_features(rec)) == ["MAP_1"]
    f = rec.features[0]
    assert f.qualifiers["pseudogene"] == ["unknown"] and "translation" not in f.qualifiers
    assert "frameshift at around 111" in f.qualifiers["note"]
    assert any(isinstance(h, PseudoGene) for h in f.secondary_hits) and f.primary_hit.id == "ARO:1"
    assert comp.hits_json["MAP_1"]["model"]["note"].startswith("pseudogene: frameshift; ")

def test_overlapping_cds_other_strand_removed_neighbour_kept(tmp_path):
    right = 100 + len(GENE)
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 110, right - 10, -1), ("MGA_2", right - 20, right + 100, 1)])
    comp.apply([[hit]])
    assert sorted(_cds_features(rec)) == ["MAP_1", "MGA_2"]

def test_no_start_with_same_gene_only_annotates(tmp_path):
    gene = "CCC" + GENE[3:]
    comp, rec, hit = _setup(tmp_path, gene, cds=[("MGA_1", 130, 100 + len(gene), 1)])
    comp.apply([[hit]])
    f = _cds_features(rec)
    assert list(f) == ["MGA_1"] and int(f["MGA_1"].location.start) == 130 and f["MGA_1"].primary_hit.id == "ARO:1"

def test_no_start_without_same_gene_is_pseudogene(tmp_path):
    comp, rec, hit = _setup(tmp_path, "CCC" + GENE[3:])
    comp.apply([[hit]])
    f = _cds_features(rec)["MAP_1"]
    assert f.qualifiers["pseudogene"] == ["unknown"] and "no start codon" in f.qualifiers["note"]

def test_invalid_reference_only_annotates(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE, validity=(False, True))
    comp.apply([[hit]])
    assert _cds_features(rec) == {}
    comp, rec, hit = _setup(tmp_path, GENE, validity=(False, True), cds=[("MGA_1", 130, 100 + len(GENE), 1)])
    comp.apply([[hit]])
    f = _cds_features(rec)["MGA_1"]
    assert int(f.location.start) == 130 and f.primary_hit.id == "ARO:1"

def test_cds_overlapping_two_regions_removed_once(tmp_path):
    gene2 = GENE + "C" * 30 + GENE
    comp, rec, hit1 = _setup(tmp_path, gene2, ref_len=len(GENE), cds=[("MGA_1", 100, 100 + len(gene2), -1)])
    hit2 = _hit(qlen=len(GENE), qend=len(GENE), sstart=101 + len(GENE) + 30, send=100 + len(gene2), slen=len(rec), btop=str(len(GENE)))
    comp.apply([[hit1], [hit2]])
    assert sorted(_cds_features(rec)) == ["MAP_1", "MAP_2"]

def test_second_instance_does_not_touch_mapped_feature(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    mapped = ExtendedFeature(location=FeatureLocation(100, 100 + len(GENE), strand=1), type="CDS", id="MAP_9", seq_id="sequence1",
                             qualifiers={"translation": ["M"]}, annotations={"partial_flag": "00", "nucl_mapping": "rescued"})
    mapped.primary_hit = NuclHit("VFG1", comp.references["ARO:1"], "VFDB", 0.0, 1.0, 99.0, 100.0, 100.0, "")
    rec.features.append(mapped); comp.genome.set_feature_dictionary()
    comp.apply([[hit]])
    assert list(_cds_features(rec)) == ["MAP_9"] and [h.id for h in mapped.secondary_hits] == ["ARO:1"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest -q tests/test_nucl_mapping.py`
Expected: FAIL (ImportError for `new_cds` etc.)

- [ ] **Step 3: Implement `new_cds`, `new_pseudogene`, `relocate_cds` and `NuclMapping.apply`**

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest -q tests/test_nucl_mapping.py`
Expected: 13 passed

- [ ] **Step 5: Commit**

```bash
git add dfc/components/NuclMapping.py tests/test_nucl_mapping.py
git commit -m "feat: rescue, relocate and pseudogene CDSs from full-length reference mappings"
```

---

### Task 7: Outputs and guards (amr_summary, pseudogene_summary, prepareQueries, fix_partial_CDS)

**Files:**
- Modify: `dfc/components/NuclMapping.py` (`write_hits`)
- Modify: `dfc/utils/summarize_amr.py:8` (glob)
- Modify: `dfc/utils/summarize_pseudo.py:17-19` (accept `NuclHit`)
- Modify: `dfc/components/baseComponent.py:76-86` (`prepareQueries` skips empty translations)
- Modify: `dfc/utils/feature_util.py:383-384` (`fix_partial_CDS` last branch)
- Test: `tests/test_nucl_mapping.py`, `tests/test_overlap_rules.py` (or a new `tests/test_partial_cds_fix.py`)

**Interfaces:**
- Produces:
  - `NuclMapping.write_hits(self) -> None`: writes `self.hits_json` to `<workDir>/nucl_hits.json` (same format as NuclSearch).
  - `summarize_amr`: `glob` pattern becomes both `NuclSearch*/*.json` and `NuclMapping*/*.json`.
  - `summarize_pseudo`: the reference hit may be a `NuclHit`; its description column uses `hit.model.to_dict()["product"]` when the hit has no `description`.
  - `prepareQueries`: a feature whose `translation` is empty is skipped (not written to the query file).
  - `fix_partial_CDS`: the final `elif isinstance(hit, NuclHit)` applies only when `feature.annotations.get("partial_flag", "00") != "00"`.

- [ ] **Step 1: Write the failing tests**

```python
import json
from dfc.utils.summarize_amr import summarize_amr
from dfc.utils.summarize_pseudo import summarize_pseudo

def test_write_hits_json(tmp_path):
    comp, rec, hit = _setup(tmp_path, GENE)
    comp.apply([[hit]]); comp.write_hits()
    data = json.load(open(tmp_path / "nucl_hits.json"))
    assert data["MAP_1"]["model"]["note"].startswith("rescued; ") and data["MAP_1"]["id"] == "ARO:1"

def test_summarize_amr_reads_nucl_mapping(tmp_path):
    work = tmp_path / "NuclMapping"; work.mkdir()
    comp, rec, hit = _setup(work, GENE)
    comp.apply([[hit]]); comp.write_hits()
    rec.features[0].qualifiers["locus_tag"] = ["LOCUS_010"]
    summarize_amr(comp.genome, str(tmp_path), {}, str(tmp_path / "amr_summary.tsv"))
    rows = [l.split("\t") for l in open(tmp_path / "amr_summary.tsv") if not l.startswith("#")]
    assert rows[0][0] == "LOCUS_010" and rows[0][-1].startswith("rescued; ")

def test_summarize_pseudo_accepts_nucl_hit(tmp_path):
    gene = GENE[:10] + GENE[11:]
    comp, rec, hit = _setup(tmp_path, gene, ref_len=len(GENE), btop="10A-" + str(len(GENE) - 11))
    comp.apply([[hit]])
    rec.features[0].qualifiers["locus_tag"] = ["LOCUS_010"]
    summarize_pseudo(comp.genome, str(tmp_path / "pseudogene_summary.tsv"))
    rows = [l.split("\t") for l in open(tmp_path / "pseudogene_summary.tsv")]
    assert rows[1][1] == "LOCUS_010" and rows[1][3] == "frameshift" and rows[1][6] == "ARO:1"

def test_prepare_queries_skips_cds_without_translation(tmp_path):
    from dfc.components.baseComponent import BaseAnnotationComponent
    comp, rec, hit = _setup(tmp_path, GENE, cds=[("MGA_1", 0, 30, 1)])
    rec.features.append(new_pseudogene(MappedRegion("sequence1", 100, 100 + len(GENE), 1, "frameshift", "", [], set(), set()), "MAP_1", 11))
    comp.genome.set_feature_dictionary()
    base = BaseAnnotationComponent.__new__(BaseAnnotationComponent)
    base.genome, base.workDir, base.CPU, base.skipAnnotatedFeatures, base.logger = comp.genome, str(tmp_path), 1, False, comp.logger
    base.query_files, base.query_sequences = {}, {}
    base.prepareQueries(split_query=False)
    assert open(base.query_files[0]).read().count(">") == 1

def test_fix_partial_cds_keeps_complete_nucl_hit_cds(tmp_path):
    from tests.test_overlap_rules import _Genome, _Config, _feature
    from dfc.utils.feature_util import FeatureUtil
    comp, rec, hit = _setup(tmp_path, GENE)
    nucl_hit = NuclHit("ARO:1", comp.references["ARO:1"], "CARD", 0.0, 500.0, 99.0, 100.0, 100.0, "")
    complete = _feature("CDS", 0, 66, "MAP_1", qualifiers={"product": ["x"], "translation": ["M"]})
    complete.annotations["partial_flag"] = "00"; complete.primary_hit = nucl_hit
    partial = _feature("CDS", 0, 66, "MGA_1", qualifiers={"product": ["x"], "translation": ["M"]})
    partial.annotations["partial_flag"] = "10"; partial.primary_hit = nucl_hit   # shorter than 500 bp: not fixable
    fu = FeatureUtil(_Genome([complete, partial]), _Config())
    fu.fix_partial_CDS(complete); fu.fix_partial_CDS(partial)
    assert complete.type == "CDS" and partial.type == "misc_feature"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest -q tests/test_nucl_mapping.py tests/test_overlap_rules.py`
Expected: the 5 new tests FAIL

- [ ] **Step 3: Implement the five changes**

- [ ] **Step 4: Run the whole suite**

Run: `python -m pytest -q tests`
Expected: all passed (no change in the existing tests)

- [ ] **Step 5: Commit**

```bash
git add dfc/components/NuclMapping.py dfc/utils/summarize_amr.py dfc/utils/summarize_pseudo.py dfc/components/baseComponent.py dfc/utils/feature_util.py tests/
git commit -m "feat: report mapped genes in the AMR and pseudogene summaries"
```

---

### Task 8: Config entries, `enable_amr_map`, `--amr_map`

**Files:**
- Modify: `dfc/default_config.py:238` (insert two `NuclMapping` entries before the first NuclSearch entry)
- Modify: `example/dfc_config/dfast_web_config.py:225` (same two entries)
- Modify: `dfc/utils/config_util.py` (`enable_amr_map`, `set_genetic_code`)
- Modify: `dfast:93-94` (argument) and `dfast:268-272` (call)
- Test: `tests/test_nucl_mapping_config.py`

**Interfaces:**
- Produces:
  - Config entries (both files), `enabled: False`, options exactly: `database` (`@@DB_ROOT@@/card/CARD` / `@@DB_ROOT@@/vfdb/VFDB`), `db_name` (`CARD` / `VFDB`), `pident_cutoff: 90`, `evalue_cutoff: 1e-10`, `overlap_tolerance: 30`, `transl_table: 11`, `aligner_options: ""`, with a comment pointing to `--amr_map`.
  - `enable_amr_map(config) -> None`: enables every `NuclMapping` entry and everything `enable_amr` enables except `NuclSearch` (reuse `enable_amr` then set the NuclSearch CARD/VFDB entries to `False`, logging a warning if any was enabled).
  - `set_genetic_code`: also sets `transl_table` of `NuclMapping` entries.
  - `dfast`: `group_workflow.add_argument("--amr_map", help="[Preliminary implementation] Like --amr, but CARD/VFDB references are mapped directly onto the genome (NuclMapping) instead of being searched with the predicted CDSs", action="store_true")`; after the `--amr` block: `if args.amr_map: enable_amr_map(config)`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_nucl_mapping_config.py
import os
from dfc.utils.config_util import load_config, enable_amr, enable_amr_map, set_genetic_code
APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _enabled(config, name):
    return [s["enabled"] for s in config.FUNCTIONAL_ANNOTATION if s["component_name"] == name]

def _config(path="dfc/default_config.py"):
    return load_config(APP_ROOT, os.path.join(APP_ROOT, path))

def test_entries_disabled_by_default():
    for path in ["dfc/default_config.py", "example/dfc_config/dfast_web_config.py"]:
        assert _enabled(_config(path), "NuclMapping") == [False, False]

def test_amr_alone_leaves_mapping_disabled():
    config = _config(); enable_amr(config)
    assert _enabled(config, "NuclMapping") == [False, False]

def test_amr_map_enables_mapping_and_plasmidfinder():
    config = _config(); enable_amr_map(config)
    assert _enabled(config, "NuclMapping") == [True, True]
    assert all(not e for e in _enabled(config, "NuclSearch"))
    assert [t["enabled"] for t in config.CONTIG_ANNOTATION if t["tool_name"] == "Plasmidfinder"] == [True]

def test_amr_map_with_amr_disables_nuclsearch(caplog):
    config = _config(); enable_amr(config); enable_amr_map(config)
    assert all(not e for e in _enabled(config, "NuclSearch")) and "NuclSearch" in caplog.text

def test_genetic_code_propagates():
    config = _config(); set_genetic_code(config, 4)
    assert [s["options"]["transl_table"] for s in config.FUNCTIONAL_ANNOTATION if s["component_name"] == "NuclMapping"] == [4, 4]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest -q tests/test_nucl_mapping_config.py`
Expected: FAIL

- [ ] **Step 3: Add the config entries, `enable_amr_map`, the `set_genetic_code` change and the CLI flag**

- [ ] **Step 4: Run tests and the CLI help**

Run: `python -m pytest -q tests/test_nucl_mapping_config.py && ./dfast -h | grep -A2 amr_map`
Expected: 5 passed; the help shows `--amr_map`

- [ ] **Step 5: Commit**

```bash
git add dfc/default_config.py example/dfc_config/dfast_web_config.py dfc/utils/config_util.py dfast tests/test_nucl_mapping_config.py
git commit -m "feat: --amr_map enables direct mapping of CARD/VFDB references"
```

---

### Task 9: Integration test with blastn

**Files:**
- Test: `tests/test_nucl_mapping_blast.py`

**Interfaces:**
- Consumes: the whole component (`NuclMapping.__init__` with a real `.nucl.ref`/`.nucl.fasta` pair written by the test, `run`).

- [ ] **Step 1: Write the test**

`@pytest.mark.skipif(shutil.which("blastn") is None or shutil.which("makeblastdb") is None, reason="BLAST+ is not on PATH")`.
The test writes, under `tmp_path`: a reference table `REF.nucl.ref` with one CARD line (`CARD.to_tabular()` of a model whose `nucl_seq` is a 600 bp gene generated from a seeded `random.Random(1)` with no in-frame stop: `"ATG" + 198 codons from {GCT, GGT, CTG, AAA, GAA, CGT, TTC} + "TAA"`), `REF.nucl.fasta` (`CARD.to_nucl_fasta()`), `REF.version`; a genome FASTA `input/genome.fna` with seven contigs (`sequence1`..`sequence7`, 300 bp `C`/`G` flanks) carrying the gene (a) intact, (b) with a 1-bp deletion at codon 50, (c) with codon 100 replaced by `TAA`, (d) intact, (e) truncated after 400 bp at the contig end, (f) with the stop codon replaced by `GCT` and `TAA` 30 bp later, (g) with `GTG` instead of `ATG`; a `FakeGenome` with `SeqRecord`s of those contigs and predicted CDSs: (a) at 30 bp downstream of the start to the stop, (b) two fragments (start..codon 50, codon 52..end), (c) start..codon 100, (d) none, (e) start..contig end (partial), (f) start..the later stop, (g) at the GTG start to the stop.
It instantiates `NuclMapping(genome, options={...}, workDir=str(tmp_path), CPU=1)` with `database=<tmp>/REF`, `db_name="CARD"`, and runs `run()`.

Assertions (by contig): (a) one CDS at the full gene, note "start codon moved 30 bp upstream by mapping of CARD:OXA-48"; (b) one pseudogene CDS, `pseudogene == ["unknown"]`, note contains "frameshift", both fragments gone; (c) one pseudogene, note contains "internal stop codon"; (d) one rescued CDS with the hit as primary; (e) the partial CDS untouched and no hit; (f) one CDS ending at the later stop with the "extended by 30 bp" note; (g) one CDS starting at the GTG, translation starting with `M`, no pseudogene. `nucl_hits.json` exists with 6 entries (a, b, c, d, f, g; e is untouched).

- [ ] **Step 2: Run the test**

Run: `python -m pytest -q tests/test_nucl_mapping_blast.py -v` (in an env with BLAST+, e.g. `dfast15-test`)
Expected: PASS (or SKIPPED without blastn). If `dc-megablast` splits the 1-bp-deletion gene into two HSPs, record it and keep the test failing until Task 10 decides the task/gap settings.

- [ ] **Step 3: Commit**

```bash
git add tests/test_nucl_mapping_blast.py
git commit -m "test: map a synthetic CARD gene onto genomes with blastn"
```

---

### Task 10: BLAST settings, validation runs and documentation

**Files:**
- Modify: `dfc/components/NuclMapping.py` (BLAST task/gap options if the measurements require it)
- Modify: `README.md:344-364` (AMR section: `--amr_map` paragraph), the CLI option list in README (`--amr` block around line 273)
- Create (untracked, under `dev/`): `dev/nucl_mapping/` with the validation runs

- [ ] **Step 1: Measure the BLAST task**

Run `NuclMapping` (through `./dfast --amr_map --cpu 8`) on `dev/new_tools/genomes/abau_aye.fna` and `ecoli_mg1655.fna` with `-task dc-megablast`, `megablast` and `blastn` (via `aligner_options` is not possible for `-task` since it is fixed: edit `blast_command` temporarily or set the task from an option `blast_task` added to the component with default `dc-megablast`). Record wall time and the number of full-length hits and of references split into several HSPs. Keep the fastest task that gives one HSP across the 1-bp deletion of Task 9. Add `blast_task` to the options (default = the chosen task) in `NuclMapping`, both config files and the spec's option list if it differs from `dc-megablast`.

- [ ] **Step 2: Validation runs**

For `example/pOXA-48.fa`, `example/sample.genome.fna`, `dev/new_tools/genomes/abau_aye.fna` and `ecoli_mg1655.fna`: run `./dfast -g <genome> -o dev/nucl_mapping/<name>_map --amr_map --cpu 8 --force` and `--amr` into `<name>_search`; compare `amr_summary.tsv` (genes only in one of them, pseudogenes, start corrections); run `ddbj-validator ddbj <out>/ddbj -o <out>/validator -f -j 1` on the `_map` outputs and confirm no error other than missing metadata. Write the comparison to `dev/nucl_mapping/README_ja.md`.
Expected: AYE shows at least one rescued or corrected gene and no validator error about `/pseudogene`; pOXA-48 reports OXA-48 as before.

- [ ] **Step 3: README**

Add to the AMR section: what `--amr_map` does (full-length mapping, rescue, start correction, `/pseudogene` with the reasons), that partial mappings are ignored, the `overlap_tolerance`/`pident_cutoff` options, and that it is planned to replace NuclSearch under `--amr`. Add `--amr_map` to the option list copied from `dfast -h`.

- [ ] **Step 4: Run the whole suite**

Run: `python -m pytest -q tests`
Expected: all passed

- [ ] **Step 5: Commit**

```bash
git add dfc/components/NuclMapping.py dfc/default_config.py example/dfc_config/dfast_web_config.py README.md
git commit -m "docs: describe --amr_map and set the blastn task from the validation runs"
```
