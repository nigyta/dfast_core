# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

DFAST (DDBJ Fast Annotation and Submission Tool) is a prokaryotic genome annotation pipeline that also generates DDBJ (MSS) submission files. Pure Python (3.10+, Biopython only) that shells out to bundled binaries in `bin/{Linux,Darwin}/` (MGA, Aragorn, Barrnap, CRT, GHOSTX). BLAST+ (`blastp`, `blastn`, `makeblastdb`, `blastdbcmd`, `rpsblast`, >= 2.13), `rpsbproc` (>= 0.5), HMMER 3 and LAST (>= 1180; `lastal -G <transl_table>`) are not bundled and must be on `PATH` (Bioconda); the bundled binaries are being phased out in favor of Bioconda. BLAST+ creates version 5 databases; existing version 4 databases remain readable.

**Language:** Write everything in tracked files in English: code, comments, docstrings, tests, commit messages and docs. Japanese is allowed only in dedicated locations such as `docs/report_ja/`. Put new Japanese documents in a `*_ja` directory like that, not alongside the English files.

## Commands

```bash
# Tests (pytest + biopython)
python -m pytest -q tests
python -m pytest -q tests/test_mag.py::test_boolean_render_true
# 5 tests in tests/test_dfast_mge.py fail unless `mefinder` (MobileElementFinder) is on PATH.

# Smoke run (needs the default protein DB: scripts/dfast_file_downloader.py --protein dfast)
./dfast --config example/test_config.py          # writes RESULT_TEST/
./dfast -g example/sample.genome.fna -o OUT --force
./dfast -g genome.fna --show_config               # print the resolved config and exit
./dfast ... --debug                               # keep temp files, dump genome.pickle
```

There is no build step and no linter config. The `dfast` script resolves `APP_ROOT` from its own real path, so it must stay at the repo root next to `dfc/`, `bin/`, and `db/`.

## Architecture

**Config is a Python file.** `dfc/utils/config_util.load_config` reads the config (default `dfc/default_config.py`, or `--config`), substitutes `@@APP_ROOT@@` / `@@DB_ROOT@@` (`--dbroot` > `$DFAST_DB_ROOT` > `APP_ROOT/db`), then `exec`s it to get the `Config` class. The workflow is the lists `STRUCTURAL_ANNOTATION`, `FUNCTIONAL_ANNOTATION`, and `CONTIG_ANNOTATION` of dicts with `tool_name`/`component_name`, `enabled`, and `options`. CLI flags in `dfast` mutate this config through helpers in `config_util.py` (`enable_amr`, `enable_mge`, `enable_mag`, `set_aligner`, `disable_*`, ...). They mostly flip `enabled` on entries already defined in `default_config.py`. A new optional step therefore needs a disabled entry in the default config plus an enable helper and a CLI flag. `example/dfc_config/dfast_web_config*.py` are the web-service configs and must be kept in sync when entries are added.

**Pipeline order** (`dfc/pipeline.py`):
1. `StructuralAnnotation`: tools run in parallel threads.
2. `FeatureUtil.execute()`: location-only overlap rules (rRNA vs assembly_gap / contig ends, rRNA vs rRNA), optional CDS merge.
3. `FunctionalAnnotation`: components run sequentially in config order.
4. `ContigAnnotation`: PlasmidFinder and MobileElementFinder.
5. `FeatureUtil.execute_after_annotation()`: overlap rules that need CDS products (rRNA/tRNA vs CDS), then `resolve_overlap()` as the fallback for overlaps no rule handles, then partial-feature removal. Locus tags are assigned after this, so removals never leave gaps in the numbering.

**Overlap rules** (`dfc/utils/feature_util.py`) target the DDBJ validator checks (ANN5310: same-strand rRNA vs CDS/rRNA overlap; ANN5320: tRNA inside a same-strand CDS). Add new rules in feature type priority order (assembly_gap > CRISPR > tRNA/tmRNA/rRNA > CDS), before `resolve_overlap()`. Removals by rules are summarized at INFO, fallback removals at WARNING, and each removal at DEBUG.
6. Add source and contig features.
7. Write outputs: GenBank, GFF, FASTA, stats, DDBJ `.ann`/`.fasta`, `pseudogene_summary.tsv`, `amr_summary.tsv`.

**Registries.** Each stage looks up classes by name in a module-level dict: `TOOLS` in `structuralAnnotation.py` / `contigAnnotation.py`, `COMPONENTS` in `functionalAnnotation.py`, and `ALIGNERS` in `components/baseComponent.py` (ghostx/ghostz/blastp/diamond/blastn). An unregistered name in the config aborts with `exit(1)`. Register any new tool or component there.

- `dfc/tools/`: wrappers around external programs, subclassing `base_tools.Tool` (structural tools subclass `StructuralAnnotationTool`, contig tools subclass `ContigAnnotationTool`, aligners subclass `Aligner`). On first instantiation each class runs `VERSION_CHECK_CMD` and matches `VERSION_PATTERN`. A missing binary calls `exit(1)`; it does not raise. That is why tests that instantiate tools need the real executable or a fake one on PATH.
- `dfc/components/`: functional-annotation steps, subclassing `BaseAnnotationComponent`. Each gets its own work dir (`WORK_DIR/<ClassName>[_N]`). Components that search references (`DBsearch`, `OrthoSearch`, ...) expose `references`, and `PseudoGeneDetection` re-aligns against those to detect frameshifts and internal stops. `skipAnnotatedFeatures` lets later searches skip CDSs that an earlier search already annotated, so config order matters.
- `dfc/genome.py`: `Genome` holds Biopython `SeqRecord`s with `ExtendedFeature` (`models/bio_feature.py`). All stages mutate this shared object.
- `dfc/models/`: parsers and models for reference hits and databases (CARD, VFDB, CDD, plasmid DB, nucleotide refs).
- `dfc/utils/ddbj_submission.py` + `metadata_util.py` + `metadata_definition.tsv`: MSS submission output. Metadata fields, their MSS feature/qualifier mapping, and validation patterns are defined in the TSV. `--metadata_file` values are loaded via `set_values_from_metadata`. MAG/ENV handling depends on `GENOME_CONFIG["project_type"]`.
- `dfc/dev/`: experimental entry points (`dfast_gff`, `dfast_re` reannotation), not part of the main pipeline.
- `scripts/dfast_file_downloader.py`: downloads and indexes reference DBs (protein, CDD, HMM, CARD/VFDB, PlasmidFinder, MobileElementFinder). It dispatches on argv at import time, so tests run it as a subprocess. `reference_util*.py` builds custom reference DBs.

## Releases

- Version string: `dfc/__init__.py` (`dfast_version`).
- Changelog: prepend an entry to `docs/history.txt` in the existing `Ver X.Y.Z (YYYYMMDD) ---` format. Update README when CLI options change.
- Pushing a bare semver git tag (e.g. `1.4.2`) triggers `.github/workflows/docker-publish.yml`, which builds and pushes `nigyta/dfast_core:<version>` and `:latest`.
- Bioconda packages the repo layout as-is.
