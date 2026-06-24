# MAG CLI Options Design

**Date:** 2026-06-24

## Goal

Add `--mag` and `--isolate` CLI options to `dfast` so that MAG/MAG-WGS submission mode can be activated without a metadata file.

## Architecture

Two files change. All existing MAG rendering logic (genome.py, metadata_util.py, ddbj_submission.py, metadata_definition.tsv) is untouched — it already branches on `project_type` in `GENOME_CONFIG`.

**Tech Stack:** Python argparse, existing config_util helper pattern.

## Global Constraints

- `--mag` alone → `project_type=mag-wgs`, `complete` unchanged
- `--mag` + `--complete t` → `project_type=mag`, `complete=True`
- `--isolate` and `--strain` are mutually exclusive (argparse-enforced)
- `--isolate` is optional; omitting it still produces a placeholder row (existing TSV `required_for` mechanism)
- CLI overrides metadata file (existing convention)
- No changes to genome.py, metadata_util.py, ddbj_submission.py, metadata_definition.tsv, or any test file other than adding new tests

---

## Design

### File: `dfc/utils/config_util.py`

Add one function, following the `enable_amr` / `enable_mge` pattern:

```python
def enable_mag(config, complete=False):
    project_type = "mag" if complete else "mag-wgs"
    config.GENOME_CONFIG["project_type"] = project_type
    if complete:
        config.GENOME_CONFIG["complete"] = True
    logger.info("MAG mode enabled (project_type={})".format(project_type))
```

Add `enable_mag` to the import line in `dfast`.

### File: `dfast` (CLI entry point)

#### Argument definitions

1. Replace the standalone `--strain` definition in `group_basic` with a mutually exclusive group containing both `--strain` and `--isolate`:

```python
group_identifier = group_basic.add_mutually_exclusive_group()
group_identifier.add_argument("--strain", help="Strain name", metavar="STR", default="")
group_identifier.add_argument("--isolate", help="Isolate name (for MAG/MAG-WGS; mutually exclusive with --strain)", metavar="STR", default="")
```

2. Add `--mag` to `group_source` (the INSDC submission section), since it only affects submission output:

```python
group_source.add_argument("--mag",
    help="Enable MAG submission mode (project_type=mag-wgs). Use with --complete t for project_type=mag.",
    action="store_true")
```

#### Processing logic

Add after the existing `--strain` processing block:

```python
if args.isolate:  # --isolate
    config.GENOME_SOURCE_INFORMATION["isolate"] = args.isolate
```

Add after the `--complete` processing block (so `complete` is already resolved in GENOME_CONFIG):

```python
if args.mag:  # --mag
    enable_mag(config, complete=config.GENOME_CONFIG.get("complete", False))
```

### Behaviour matrix

| CLI | project_type | complete |
|-----|-------------|---------|
| `--mag` | `mag-wgs` | unchanged |
| `--mag --complete t` | `mag` | `True` |
| `--metadata_file` with `project_type=mag` | `mag` | `True` |
| `--mag` + `--metadata_file` with `project_type=wgs` | `mag-wgs` (CLI wins) | unchanged |
| `--isolate foo` (no `--mag`) | project_type unchanged; isolate set | unchanged |

### Impact

- Existing 46 MAG tests: unaffected (all use metadata-file path)
- `--strain` help text position: minor cosmetic change (moves into exclusive group)
- `--metagenome` and `--mag` are orthogonal; both can be specified together

---

## Testing

New tests in `tests/test_mag.py` (or a new `tests/test_mag_cli.py`):

1. `enable_mag(config, complete=False)` → `project_type="mag-wgs"`, `complete` unchanged
2. `enable_mag(config, complete=True)` → `project_type="mag"`, `complete=True`
3. `enable_mag` after metadata sets `project_type="wgs"` → overwrites to `"mag-wgs"`
4. argparse rejects `--strain foo --isolate bar` (mutual exclusion)
5. `--isolate foo` alone → `config.GENOME_SOURCE_INFORMATION["isolate"] == "foo"`
