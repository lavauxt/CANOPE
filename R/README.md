# CANOPE — Enhanced Edition

This is your CANOPE pipeline with bugs fixed and the ECHO feature set ported
over. **`hmm_engine.R` is byte-for-byte unchanged** — every fix and new
feature works *around* the HMM math, never inside it.

> **Update:** two more bugs found and fixed after the first round — the
> report couldn't find its own template (`canope_report.R` looked in the
> wrong installed-package subdirectory), and the NB predictive-interval /
> z-score panels showed the test sample far outside the interval even
> without a real CNV (a depth-scale mismatch between raw and internally-
> normalised counts). See "Bugs fixed" and "Round 2" below for details.

## Files

| File | Status |
|---|---|
| `hmm_engine.R` | **Unchanged** — the engine core, as requested |
| `call_cnvs.R` | Bug-fixed (round 1 + round 2) |
| `data_utils.R` | Bug-fixed |
| `generate_plots.R` | Bug-fixed (round 1 + round 2 — z-score/CI panels) |
| `canope_report.R` | Bug-fixed (round 2 — template path) |
| `CANOPE_report.Rmd` | Bug-fixed (round 2 — same CI/z-score fix as generate_plots.R) |
| `run_canope.R` | Bug-fixed (round 1 + round 2) |
| `qc_reference_utils.R` | Minor robustness only, logic untouched |
| `run_canope.R` | Bug-fixed + orchestrates all new features |
| `canope_utils.R` | **New** — shared helpers (PCA plot, exon numbering, parsers...) |
| `canope_qc_metrics.R` | **New** — QC metrics table (ECHO's `metrics.R` ported) |
| `canope_confidence.R` | **New** — HIGH/MEDIUM/LOW confidence scoring |
| `canope_vcf_export.R` | **New** — VCF export |
| `canope_bed_process.R` | **New** — BED preprocessing (ECHO's `process_bed.R` ported) |
| `canope_report.R` | **New** — renders the HTML report |
| `canoes_legacy_engine.R` | **New** — optional faithful port of the original CANOES HMM core |
| `CANOPE_report.Rmd` | **New** — interactive HTML report template |

## Bugs fixed

**`run_canope.R`** (6 bugs):
1. `reads_file` ingestion renamed columns starting at index 4 instead of 5,
   silently destroying the `GENE` column for any pre-computed counts table
   using the standard `chromosome,start,end,GENE,<samples>` layout.
2. Auto-coverage extraction deduplicated `coords` and then sliced
   `counts_df` by row position — any duplicate BED interval misaligned
   every coordinate with the wrong count row from that point on.
3. `modechrom`/`removeY` filtering hardcoded `"chrX"`/`"chrY"`, returning
   zero rows for BED files using bare `"X"`/`"Y"` names.
4. GC content wasn't normalised to a 0–1 fraction before the
   `gc_extreme_filter` step — a 0–100-scale external GC table would have
   nearly every target filtered out.
5. GC column auto-detection assumed the value was always column 4.
6. `save(..., bed_file = bed_file_out, ...)` — **this never actually
   worked**. R's `save()` does not rename objects via a named argument; the
   RData never contained an object literally named `bed_file`, which is why
   `generate_plots()` always had to "reconstruct" it from `counts`. Now an
   actual `bed_file` variable is created and saved properly.

**`call_cnvs.R`** (2 bugs):
1. **The big one** — CNV calls were grouped by row-index contiguity only,
   never checking whether the chromosome changed between two adjacent rows.
   Since the HMM models the whole genome as one sequence, an abnormal-state
   run could straddle the last target of one chromosome and the first of
   the next, producing one fake CNV record spanning two chromosomes
   (mislabelled with only the first row's chromosome). Verified fixed with
   a real HMM run injecting a deletion exactly at a chromosome boundary —
   it now correctly splits into two separate calls.
2. Unparseable chromosome labels (e.g. `chrM`) became `NA` after the
   numeric remap and were silently carried through sorting/distance
   calculations. Now dropped with a warning.

**`data_utils.R`** (2 bugs):
1. `targets_to_rows()` had a dangerous fallback: if a target ID didn't
   match, it treated the ID itself as a row index — which is wrong as soon
   as any QC filtering has removed rows. Now returns `NA` with a warning,
   like every caller already expects.
2. `get_coverage_from_bams()` always called `summarizeOverlaps(singleEnd =
   TRUE)` regardless of how the BAM was actually sequenced. Most
   exome/panel BAMs are paired-end. Now auto-detected (overridable).

**`generate_plots.R`** (1 bug): the z-score panel divided by `sd()` of a
single reference sample, which is undefined (`NA`) — with only one
reference, the panel silently rendered blank. Now falls back to a robust
MAD-based variance estimate and logs why.

**`canope_qc_metrics.R`** (found during testing, now fixed): the QC-metrics
helper attached `._gene`/`._chrom` columns to `bed_df` for its own use, then
passed that same `bed_df` into `assign_exon_numbers_per_gene()`, which uses
those exact names internally — a silent `data.table` column collision. Now
uses local variables instead of mutating the data frame.

## Round 2 bug fixes

Two more bugs turned up after actually trying to use the report and the
per-call plots.

### 1. Report couldn't be found (`canope_report.R`)

`generate_canope_report()`'s fallback template search looked for the
installed package's copy at `system.file("rmarkdown/CANOPE_report.Rmd",
package = "CANOPE")`. The template actually ships at `inst/rmd/`, which
`system.file()` resolves to `.../rmd/CANOPE_report.Rmd` once installed —
not `.../rmarkdown/...`. `system.file()` returns `""` (not an error) for a
subdirectory that doesn't exist, so this candidate silently never matched
anything, and unless a copy of the Rmd happened to sit next to your
`rdata_output` or in the working directory, report generation failed with
"Could not locate CANOPE_report.Rmd". Fixed to look under `rmd/`.

### 2. NB predictive interval / z-score panel showed the test sample way
   outside the interval even with no CNV

**Root cause:** `call_cnvs()` normalises coverage before fitting anything —
first a global median normalisation across *all* samples (test included),
then it rescales each *reference* sample's already-normalised counts to
match the test sample's own median. `counts$mean` (→ `models[[sample]]
$mean`) and the MAD-based `var_estimate` are computed from that doubly-
normalised data. But all of that normalisation only ever happened to
`call_cnvs()`'s local copy of `counts` — the caller's original data frame,
which `run_canope()` saves verbatim to the RData workspace as `counts`,
never saw it.

`generate_plots.R` (and the identical logic duplicated in
`CANOPE_report.Rmd`) built the "observed" side of both the NB
predictive-interval ratio and the z-score panel by pulling straight from
that raw, un-normalised `counts` object — comparing it against
`model_mean`/`model_var`, which are on a different depth scale. Since
`call_cnvs()`'s normalisation exists specifically to correct for
sample-to-sample sequencing-depth differences, skipping it here reintroduced
exactly that depth bias into both panels: a sample sequenced deeper or
shallower than its references would show a large, systematic, entirely
spurious deviation — a fake "CNV signal" from a normal region, which is
what you were seeing.

Confirmed with a simulation of 5 samples sharing the same true coverage but
different sequencing depths (0.5×–2×) and no injected CNV: the old
raw-vs-model comparison averaged a **log2 ratio of +1.03** (≈2× — a
convincing but entirely fake duplication signal) across a flat region;
using values on a consistent scale gives ~0, as it should.

**Fix:** `call_cnvs(..., full_output = TRUE)` now also returns
`test_counts` and `ref_matrix` — the test sample's counts and the reference
matrix exactly as they stood when `mean`/`var_estimate` were computed.
`run_canope()` stores these in `models[[sample]]` alongside the existing
`mean`/`var_estimate`. Both `generate_plots.R`'s CI/z-score panels and the
equivalent code in `CANOPE_report.Rmd` now read `test_counts`/`ref_matrix`
from the model (via a new `model_matrix_lookup()` helper, matching the
existing target-ID-based `model_lookup()`) instead of re-deriving from raw
counts — so the compared values are always on the scale the model was
actually fit on. The raw-coverage panel at the top of each plot (the plain
log2 coverage traces) is untouched; it already does its own local
median-ratio scaling for display and wasn't reported as wrong.

**Caveat:** an RData workspace written by the *old* `run_canope()` won't
have `test_counts`/`ref_matrix` in its `models`. `generate_plots()` now
checks for this up front and stops with a clear message telling you to
re-run `run_canope()`, rather than failing deep inside a per-call
`tryCatch` with a vague NULL-subscript error.

## New features (ECHO parity)

- **BED preprocessing** (`canope_bed_process.R`) — `STANDARD`/`REGEN`/`NO`
  modes, gene/exon annotation from panel files or RefSeq.
- **QC metrics table** (`canope_qc_metrics.R`) — per-sample and per-exon
  flags (low correlation, low depth, high CV, missing sex chromosomes),
  written to TSV and shown in the report.
- **Confidence scoring** (`canope_confidence.R`) — HIGH/MEDIUM/LOW based on
  `Q_SOME`, `NUM_REFS`, `NUM_TARG`, and MLCN consistency, with a built-in
  list of pseudogene/homology-prone genes always scored LOW.
- **PCA of coverage profiles** (`canope_utils.R::plot_coverage_pca`).
- **VCF export** (`canope_vcf_export.R`) — multi-sample or per-sample.
- **Interactive HTML report** (`CANOPE_report.Rmd` /
  `canope_report.R::generate_canope_report()`) — QC table, PCA, per-sample
  CNV tables colour-coded by confidence, and three plotly panels per call:
  coverage, log2(observed/expected) with NB confidence interval, **and a
  new z-score-vs-references panel** (the same z-score concept that was
  already in your static PDFs, now also interactive).

All new features are wired into `run_canope()` as toggleable steps
(`run_qc_metrics`, `score_confidence`, `pca_plot`, `export_vcf`, `report`,
...) so existing calling code keeps working — just call `run_canope()` as
before and the extras run automatically with sensible output paths next to
your `output_file`.

## Fast coverage extraction (megadepth, Windows-native)

`data_utils.R` also has `get_coverage_from_bams_megadepth()` — a much faster
alternative to the default `summarizeOverlaps()`-based extraction, for when
`run_canope()` is pulling coverage straight from BAMs (i.e. no `reads_file`).

Of the usual fast-depth tools (mosdepth / megadepth / PanDepth), **megadepth
is the only one with an official native Windows binary** — confirmed in its
Bioinformatics paper and Bioconductor's own Windows build reports. mosdepth
and PanDepth are Linux/macOS-only upstream and need WSL2 on Windows. The
Bioconductor `megadepth` package's `install_megadepth()` downloads the right
binary for whatever OS R is on, Windows included, so there's nothing to
compile or hand-place on `PATH`.

Switch backends with one argument:

```r
run_canope(
  ...,
  coverage_backend = "megadepth",   # default is "bioconductor" (unchanged behaviour)
  megadepth_op = "sum",             # or "mean"
  megadepth_threads = 4
)
```

Note megadepth reports base-level *coverage* (summed/averaged depth), not a
fragment count the way `summarizeOverlaps()` does. That's fine for CNV
calling — everything downstream (GC correction, reference normalisation, the
HMM) works on relative deviations, not literal Poisson counts — but don't mix
the two backends for the same target across one analysis.

**Before trusting it on a full cohort**, validate on one BAM:

```r
BiocManager::install("megadepth")  # one-time
library(megadepth); install_megadepth()
source("data_utils.R")
test_counts <- get_coverage_from_bams_megadepth("one_sample.bam", "panel.bed")
head(test_counts)
```

The output-parsing logic was unit-tested against megadepth's documented
`chr,start,end,score` schema; the binary itself wasn't exercised end-to-end
in the environment this was written in (no network path to GitHub's
release-asset host from there), and megadepth's exact CLI flags have shifted
across releases before — hence the one-BAM sanity check.

## Optional legacy HMM engine (original CANOES)

`call_cnvs()` and `run_canope()` now take `engine = c("new", "legacy_canoes")`
(default `"new"`, unchanged behaviour). `"legacy_canoes"` switches to
`canoes_legacy_engine.R` — a faithful port of the original published CANOES
HMM core (Backenroth et al. 2014).

**What's actually different between the two**, since the underlying
statistical model turned out to already be the same one (the "new" engine in
`hmm_engine.R` is itself a refactor of CANOES, not a different model):

| | `"new"` (default) | `"legacy_canoes"` |
|---|---|---|
| Transition/emission formulas | identical | identical |
| Inter-target distance | start-to-end, floored at 0 | start-to-start, no floor (original) |
| Stationary transition mode | available | not available (always distance-decayed, as published) |
| Phred score | floored at `1e-10` to avoid NaN | unguarded, exactly `-10*log10(1-prob)` — can emit `NaN` on edge-case probabilities, faithfully |
| Emission all-impossible-state fallback | soft (`log(1e-12)`/`0`) | hard (`-Inf`/`-0.01`, the original's choice) |
| Initial-state abnormal prior | tunable via `hmm_params(prior_abnormal=...)` | hardcoded `0.0075`, as published |
| Reference selection, NNLS weighting, variance estimation | shared, unchanged | shared, unchanged |

`canoes_legacy_engine.R` was checked function-by-function against your
actual uploaded `CANOES.R` (`GetDistances`, `GetTransitionMatrix`,
`EmissionProbs`, `Viterbi`, `GetForwardMatrix`/`GetBackwardMatrix`,
`GetModifiedLikelihood`, and the inline Phred formula) on identical synthetic
inputs — every one matched exactly, including the `NaN`-producing Phred edge
case. A full `call_cnvs()` run with an injected deletion and duplication was
also run under both engines end-to-end and detected the same calls.

```r
run_canope(..., engine = "legacy_canoes")   # or per-sample: call_cnvs(..., engine = "legacy_canoes")
```

If you request `decode_method = "stationary"` together with
`engine = "legacy_canoes"`, you'll get a warning and it falls back to
distance-based transitions — the original never had a stationary mode.

Note this only swaps the HMM core (transition/emission/decoding/Phred).
Reference-sample selection, NNLS bootstrap weighting, and per-target variance
estimation are shared between both engines and weren't touched — the original
CANOES used a GAM-smoothed `var ~ s(mean) + s(gc)` variance estimator instead
of the MAD-based one CANOPE uses now; that's a separate, identifiable piece I
didn't fold into this toggle. Happy to add it as a further option if you want
a fully byte-for-byte original pipeline rather than just the original HMM.

## Quick example

```r
source("hmm_engine.R")
source("canope_utils.R")
source("data_utils.R")
source("qc_reference_utils.R")
source("call_cnvs.R")
source("canope_confidence.R")
source("canope_qc_metrics.R")
source("canope_vcf_export.R")
source("canope_bed_process.R")
source("generate_plots.R")
source("canope_report.R")
source("run_canope.R")

results <- run_canope(
  fasta_file  = "hg38.fa",
  bed_file    = "panel.bed",
  samples     = "bams/",
  output_file = "results/CNVCall.csv",
  rdata_output = "results/canope_workspace.RData"
)
```

This single call now also writes, next to `CNVCall.csv`:
`CANOPE_QC_metrics.tsv`, `CANOPE_PCA.pdf`, `CANOPE_calls.vcf`,
`canope_pipeline.log`, the per-call PDFs in `plots/`, and an interactive
`report/CANOPE_CANOPE_report.html`.

## Testing performed

Everything was actually exercised, not just read — `r-cran-{data.table,
dplyr,matrixstats,nnls,ggplot2}` were installed to run real synthetic data
through the **unmodified** HMM engine:

- Injected a deletion and duplication into synthetic coverage data → both
  correctly detected and confidence-scored.
- Injected a deletion straddling a chromosome boundary → confirmed it now
  produces two separate, correctly-labelled calls instead of one fake
  cross-chromosome call.
- Exercised `targets_to_rows()`'s unmatched-ID path, `get_distances()`'s
  chromosome-break sentinel, `assign_exon_numbers_per_gene()`'s ordering and
  de-duplication, the single-reference z-score fallback, QC metrics
  generation, and VCF formatting.
- Every file parses cleanly and sources without error.

**Round 2:** all R files (including the code chunks inside
`CANOPE_report.Rmd`) were re-parsed cleanly after the fixes above.
`model_lookup()`/`model_matrix_lookup()` were unit-tested directly —
reordering, subsetting, and median-imputing a missing target ID both
behave correctly for a vector field (`test_counts`) and a matrix field
(`ref_matrix`). The depth-scale mismatch itself was reproduced and
confirmed fixed with the 5-sample/no-CNV simulation described above (old
approach: +1.03 mean log2 ratio; fixed approach: ~0). The actual
`rmarkdown::render()` call and a full `run_canope()` → `generate_plots()`
round-trip weren't exercised end-to-end in this environment (no BAM/FASTA
test data or `rmarkdown`/`ggplot2`/etc. available here beyond base R) —
worth a quick real run before trusting it on a full cohort, same as the
megadepth note above.
