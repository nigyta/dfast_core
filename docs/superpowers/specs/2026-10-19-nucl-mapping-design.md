# Direct mapping of CARD/VFDB references onto the genome (NuclMapping)

**Date:** 2026-10-19
**Branch:** feature/amr-map (from master 1.5.3)
**Status:** design, reviewed (builder check 2026-10-19), awaiting approval

## Intent

CARD and VFDB annotation (`--amr`) currently runs `NuclSearch`: the nucleotide sequence of each
predicted CDS is searched with blastn against the CARD/VFDB nucleotide references, and a CDS gets
the reference's product when identity, query coverage and subject coverage are all above 90%.
Because the query is the predicted CDS, the result depends on CDS prediction:

- a gene that the gene finder did not predict is never annotated;
- when the gene finder picked the wrong start codon among several candidates, the wrong start is kept;
- a gene broken by a frameshift or an in-frame stop codon is predicted as a truncated CDS or as
  several fragments, each of which fails the coverage cutoffs or is annotated as a complete gene.

This design maps the reference nucleotide sequences directly onto the genome with blastn instead.
A reference that maps over its full length at high identity gives start-to-stop coordinates that do
not depend on gene prediction, so the gene can be rescued, its start codon taken from the reference,
and a frameshift or in-frame stop recognized and reported as a pseudogene. Direct nucleotide mapping
is only meaningful for very close sequences with trusted references, which is already the condition
under which NuclSearch accepts a hit (90% identity), so the mapping is expected to replace NuclSearch
for CARD/VFDB.

The same approach was used for the SARS-CoV-2 pipeline (`dfc/tools/blastfeature_n.py` on the
`covid19` branch: reference features are mapped onto the genome with blastn and turned into features
with strand, partial flags, translation and transferred qualifiers). That code is the model for the
mapping and parsing part; it is not reused as is because it is a structural annotation tool that
knows nothing about predicted CDSs.

## Goals

1. Rescue CARD/VFDB genes that have no predicted CDS.
2. Correct the start codon of a predicted CDS when a reference maps over its full length.
3. Report a gene with a frameshift or an in-frame stop codon as one CDS with `/pseudogene`, in place
   of the truncated or fragmented predictions.
4. Keep `amr_summary.tsv` and the GenBank/DDBJ output for CARD/VFDB hits in the current form.
5. Run with a new option `--amr_map` during development. Once validated, `--amr` switches to the
   mapping and NuclSearch stays available through the config file.

## Non-goals

- Partial mappings. A reference that does not map over its full length is ignored by this component:
  no feature is created, changed or removed. (The predicted CDS in that region keeps whatever other
  components give it.) Detecting IS insertions in a gene from a split mapping is a possible later
  extension and is noted at the end; it is not part of this work.
- Removing or changing NuclSearch. It stays as a component and may be used for other purposes.
- Changing PseudoGeneDetection.
- Databases other than CARD and VFDB. The component takes any reference prepared by
  `reference_util_for_nucl.py` (a `NucRefBase` model), but only CARD and VFDB entries are added to
  the configs.

## Anti-goals

- Never two CDS features for the same gene (a mapped one next to the predicted one).
- No change at all to the results when neither `--amr_map` nor the config entries enable it.
- No lowering of the identity cutoff. The user's premise is that mapping is only used where
  NuclSearch would already accept the hit.

## Constraints

- Pure Python 3.10+ with Biopython; `blastn` and `makeblastdb` (BLAST+ >= 2.13) on PATH, as for
  NuclSearch.
- Follows the component pattern: a `BaseAnnotationComponent` subclass in `dfc/components/`,
  registered in `functionalAnnotation.COMPONENTS`, with disabled entries in `dfc/default_config.py`
  and `example/dfc_config/dfast_web_config.py` (kept in sync), an `enable_*` helper in
  `dfc/utils/config_util.py` and a CLI flag in `dfast`.
- Everything in tracked files is English.
- Tests with pytest; the parsing and feature logic must be testable without blastn.

## Decisions made by the user

- Map CARD/VFDB references directly onto the genome with blastn; the mapping is intended to replace
  NuclSearch for these databases.
- Only full-length mappings at high identity act on the annotation; partial mappings do nothing.
- A full-length mapping with a frameshift or an in-frame stop codon is annotated as a CDS with
  `/pseudogene` (with this evidence the call is considered justified; the existing
  PseudoGeneDetection only writes a "possible pseudo" note).
- A lost stop codon is tolerated up to a 100 bp extension (then a pseudogene), and an in-frame stop
  in the last two sense codons shortens the CDS instead of making it a pseudogene.
- A gene annotated by the mapping has priority over predicted CDSs that overlap it: the predicted
  CDS is removed (or relocated, when it is the same gene).
- Development flag `--amr_map`; `--amr` switches to it after validation.

## Design

### Pipeline position

`NuclMapping` is a functional annotation component. In `FUNCTIONAL_ANNOTATION` its two entries
(CARD, VFDB) come right before the NuclSearch entries, i.e. before DBsearch, so that:

- a CDS created or annotated by the mapping has a primary hit and is skipped by the later
  DBsearch/HMMscan (`skipAnnotatedFeatures`), exactly as a NuclSearch hit is today;
- CDDsearch, PseudoGeneDetection and the overlap rules (`FeatureUtil.execute_after_annotation`)
  run afterwards and see the final CDS set;
- locus tags are assigned after functional annotation, so rescued CDSs are numbered in place.

The component mutates the genome: it adds, relocates and removes CDS features of
`genome.seq_records[...].features` and refreshes `genome.features` (`Genome.set_feature_dictionary`,
`Genome.sort_features`). This is new for a functional annotation component; the structural stage and
`FeatureUtil` already do it, so the methods exist.

### Options

Config entry (same shape as NuclSearch):

```python
{
    "component_name": "NuclMapping",
    "enabled": False,
    "options": {
        "database": "@@DB_ROOT@@/card/CARD",   # prefix of <db>.nucl.fasta / <db>.nucl.ref / <db>.version
        "db_name": "CARD",                       # CARD or VFDB (models in NuclSearch.models)
        "pident_cutoff": 90,                     # minimum identity of a mapping (%)
        "evalue_cutoff": 1e-10,
        "overlap_tolerance": 30,                 # predicted CDSs overlapping a mapped gene by more than this (bp) are removed
        "aligner_options": "",                   # extra blastn arguments
    },
},
```

`pident_cutoff` 90 is the current NuclSearch value; the coverage cutoffs are not options because
full length is required.

CLI: `--amr_map` calls `enable_amr_map(config)`, which enables the `NuclMapping` entries for CARD and
VFDB and everything else `--amr` enables except NuclSearch (PlasmidFinder today; DefenseFinder once
`feature/integron-defense` is merged, by adding it to the helper). If `--amr` and `--amr_map` are both
given, NuclSearch is disabled with a warning, so that a CDS is not reported twice in
`amr_summary.tsv`. Switching `--amr` to the mapping later is a change in `enable_amr` only.

### Mapping

1. Build a BLAST database of `WORK_DIR/input/genome.fna` in the component's work directory
   (`makeblastdb -dbtype nucl`); it takes a second. Whether the two instances (CARD, VFDB) share one
   database is the builder's choice.
2. `blastn` with the reference `<database>.nucl.fasta` as query and the genome as database, `-outfmt 6`
   with `qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen
   slen btop`, `-num_threads CPU`, enough alignments to report all loci of a reference
   (several contigs and several copies: `-max_target_seqs` and `-max_hsps` well above 1), and
   `-perc_identity pident_cutoff` as a prefilter. The BLAST task is the builder's choice
   (see "Left to the builder").
3. Full-length test: a hit is full-length when the alignment reaches to within one codon of both
   ends of the reference (`qstart <= 4` and `qend >= qlen - 3`) and `pident >= pident_cutoff`. The
   slack is needed because a BLAST alignment never ends on a mismatch: a GTG start where the
   reference has ATG, or a TAG stop where it has TAA, leaves the first or last bases unaligned. The
   mapped region is extrapolated to the full extent of the reference (`qstart - 1` bases before the
   alignment and `qlen - qend` after it), and the validity test below looks at the genome's own
   codons there. If a reference is covered by several collinear HSPs on the same contig and strand
   (an indel longer than the gap penalties allow in one HSP), they may be chained into one
   full-length mapping; see "Left to the builder".
4. The mapped region is the extrapolated `[min(sstart, send), max(sstart, send)]` on the strand
   given by `sstart <= send`. If the reference itself ends without a stop codon (about 3% of CARD
   and 1% of VFDB references), the region is extended to the next in-frame stop codon in the genome;
   if that is more than 100 bp away, the hit is not full-length.
5. One reference per locus: full-length hits whose mapped regions overlap on the same contig and
   strand are grouped, the highest bit score (then identity) is the primary hit, the others become
   secondary hits of the same feature ("alternative hit" notes at verbosity >= 2, as NuclSearch does).
   A reference may map to several loci; each locus is handled separately.

### Reference validity

The user's premise is that direct mapping is only justified when the reference is known to be
correct. A reference is **valid** when its length is a multiple of 3, it starts with a start codon
of the translation table, its translation has no internal stop codon and it has no ambiguous bases;
a terminal stop codon is optional (step 4 above). About 2% of CARD and 3% of VFDB references are not
valid (no start codon, length not a multiple of 3, an internal stop, or Ns). A full-length mapping of
an invalid reference is **annotate-only**: it attaches its hit to a predicted CDS it overlaps in the
same frame (as NuclSearch does today) and never rescues, relocates or calls a pseudogene. Validity is
computed once when the reference table (`<database>.nucl.ref`) is loaded.

### What a full-length mapping does

The genome sequence of the mapped region is translated with `transl_table` 11 (the value used by
PseudoGeneDetection; taken from the config if the pipeline carries one). The region is a **valid CDS**
when its length is a multiple of 3, its first codon is a start codon of the table, the translation has
no internal stop and ends with a stop. Otherwise the region carries a **defect**: a frameshift (the
alignment has gaps whose total length is not a multiple of 3, or the chained HSPs shift frame), an
internal stop codon, or a lost stop codon beyond the tolerance below.

Two tolerances at the 3' end (decided by the user):

- **Lost stop codon.** When the genome has no stop codon at the reference's stop position (same
  frame, read-through), the region is extended to the next in-frame stop codon. An extension of at
  most 100 bp gives a valid CDS with a longer C terminus and the note "stop codon of the reference is
  not conserved; extended by N bp"; the same-gene test then uses the extended stop. If there is no
  stop within 100 bp, the region is a pseudogene with the note "stop codon lost".
- **Stop in the last codons.** An in-frame stop within the last two sense codons of the reference
  (a protein shorter by at most two residues) is not a defect: the region is shortened to that stop
  and is a valid CDS with the note "N codon(s) shorter than the reference". An earlier stop is an
  internal stop codon (pseudogene).

Predicted CDSs on the same contig that overlap the mapped region are classified. The mapping has
priority (decided by the user): a predicted CDS gives way to the mapped gene unless the overlap is
small enough to be two adjacent genes.

- **same gene**: on the same strand, same reading frame and the same stop codon (same right end on
  the plus strand, same left end on the minus strand);
- **overlapping**: any other predicted CDS, on either strand, whose overlap with the mapped region is
  longer than `overlap_tolerance` (default 30 bp): a truncated or frame-shifted piece of the gene, a
  misprediction on the gene, or a spurious antisense ORF; removed;
- **neighbour**: an overlap of at most `overlap_tolerance` (adjacent genes of an operon overlapping by
  a few bases, convergent genes); untouched.

Then:

| mapped region | same-gene CDS exists | action |
|---|---|---|
| valid CDS | no | **rescue**: a new CDS at the mapped region with the reference's product/gene/note (the `NuclHit` becomes its primary hit) |
| valid CDS | yes, same location | annotate: the hit becomes the primary hit, nothing else changes (today's NuclSearch result) |
| valid CDS | yes, different start | **start correction**: the CDS location is set to the mapped region, the translation is recomputed, a note records the change ("start codon moved N bp upstream/downstream by mapping of CARD:..."); the feature keeps its id |
| defect | any | **pseudogene**: one CDS at the mapped region with `/pseudogene="unknown"`, the reference's product/gene/note, a note with the reason and position(s) ("frameshift at around 1234", "internal stop codon at 567", "stop codon lost"); no `/translation`; the same-gene CDS, if any, is removed |

In every case the **overlapping** CDSs are removed. A region that is intact except that its first codon
is not a start codon is handled like this: if a same-gene CDS exists, it is only annotated (no start
correction, since the reference start is not usable); if none exists, the region is a pseudogene with
the note "no start codon".

Primary hit: when the mapping changes the feature (rescue, start correction, pseudogene) its
`NuclHit` becomes the primary hit, replacing an earlier one (e.g. from OrthoSearch with
`--references`), whose coverage values would be stale after a relocation. When it only annotates an
existing CDS, it follows NuclSearch: primary if the CDS has none, otherwise secondary.

`FeatureUtil.fix_partial_CDS` currently turns any CDS within 10 bp of a contig end whose primary hit
is a `NuclHit` into a `misc_feature`, even a complete one (its first branches handle genuinely partial
CDSs; the last `elif` has no partial check). A rescued or corrected CDS near a contig end would lose
its product. That branch is limited to CDSs whose `partial_flag` is not "00"; this also stops the
same demotion of complete NuclSearch-annotated CDSs, which is a pre-existing quirk.

RNA features overlapping a mapped gene are handled by the existing overlap rules
(`FeatureUtil.execute_after_annotation`): a CDS with a product wins over an rRNA/tRNA on the same
strand, so a mapped CDS or pseudogene (both carry a product) also has priority there. Removals,
relocations and rescues are logged at INFO as counts and at DEBUG per feature, like the overlap rules.

### Representation

- Rescued and corrected CDSs are ordinary `ExtendedFeature`s of type `CDS` with `translation`,
  `transl_table`, `codon_start` and the reference qualifiers set through `NuclHit.assign`
  (`model.set_info_to_feature`, which gives CARD its gene family as product and the "similar to ..."
  note). They get locus tags as any CDS.
- A pseudogene is a `CDS` feature with `/pseudogene="unknown"` (INSDC controlled vocabulary; the
  reason is in `/note`), no `/translation`, and the reference qualifiers. It gets a locus tag. It is
  written to GenBank, GFF and the DDBJ `.ann` (the writers copy qualifiers as they are) and counted as
  a CDS in `statistics.txt`.
- Components that build protein queries skip a CDS without `translation`
  (`BaseAnnotationComponent.prepareQueries`, a one-line guard). PseudoGeneDetection only examines
  CDSs with a `ProteinHit`, so it ignores mapped CDSs and pseudogenes without any change.
- `amr_summary.tsv`: the component writes `nucl_hits.json` in its work directory in the NuclSearch
  format (`feature_id -> NuclHit.to_dict()`), and `summarize_amr` also globs `NuclMapping*`. The
  `note` column of a pseudogene line starts with "pseudogene: frameshift" / "pseudogene: internal stop
  codon". A rescued CDS is marked "rescued" and a corrected one "start corrected" in that column.
- `pseudogene_summary.tsv`: a mapping pseudogene also gets a `PseudoGene` secondary hit (reference id,
  internal stop positions, indel positions), so `summarize_pseudo` lists it; the writer, which today
  looks up a matching `ProteinHit` for the reference columns, is extended to accept the `NuclHit`.

### Testing

- Unit tests (no blastn): the hit table parser and full-length test; grouping of overlapping hits;
  classification of a mapped region (valid / frameshift / internal stop / no start) from a synthetic
  genome; classification of predicted CDSs (same gene / overlapping / neighbour, both strands); the
  feature actions.
  Fixtures are a small synthetic genome and a blastn `-outfmt 6` table written by hand or cut from a
  real run on pOXA-48 (`example/pOXA-48.fa`, which has OXA-48 in CARD).
- One integration test that runs blastn on a synthetic genome containing a CARD reference (a) intact
  with a downstream ATG as the predicted start, (b) with a 1-bp deletion, (c) with an internal stop
  codon, (d) with no predicted CDS, (e) truncated at the contig end, (f) with the stop codon mutated
  and the next in-frame stop within 100 bp, (g) with a GTG start where the reference has ATG; skipped
  when blastn is not on PATH (as other tool tests do).
- `--amr_map` enables the entries; `--amr` alone leaves them disabled.
- Full runs for validation (not tests): `example/pOXA-48.fa`, `example/sample.genome.fna`, A.
  baumannii AYE (CU459141.1) and E. coli MG1655 (dev genomes under `dev/new_tools/genomes/`),
  compared with `--amr` (NuclSearch) on `amr_summary.tsv`, plus ddbj-validator on the `.ann`.

### Documentation

README (a paragraph under the AMR section describing `--amr_map`, what full-length mapping does and
that it is planned to replace NuclSearch under `--amr`), `dfast -h`, `docs/history.txt` at release.

## Left to the builder

- The blastn task (`megablast`, `dc-megablast` or `blastn`): megablast is fastest and enough at 90%
  identity, but may fragment an alignment around an indel; measure on AYE and MG1655 (41,000
  references against a 4 Mb genome must stay well under a minute) and pick the task and gap settings
  that give one full-length HSP across a 1-bp indel.
- Whether to chain collinear HSPs into one mapping, and the maximum gap allowed between them. Start
  without chaining; add it if the validation runs show references split by larger indels.
- How to locate frameshift positions (from `btop` or from `qseq`/`sseq`) and internal stop positions
  (from the translation); positions in the note are genome coordinates.
- The handling of a predicted CDS that overlaps two mapped regions, and of a predicted CDS that
  already has a primary hit from an earlier component (it is removed like any other; the hit is
  dropped with it).
- Coordination of the two instances (CARD then VFDB): a feature created or relocated by the first
  instance is not changed by the second; a second full-length mapping with the same location is
  attached as a secondary hit, one with a different location on an overlapping region is logged and
  skipped.
- Whether the genome BLAST database is built once and shared by both instances.

## Future extension (not in scope)

A reference that maps in two collinear parts with a large insertion between them (an IS element
inside the gene) could be reported as a pseudogene with the note "interrupted by an insertion", using
the MobileElementFinder/IntegronFinder features to name the insertion. This needs the chaining logic
above and a decision on how far apart the parts may be; it is left for a later design.
