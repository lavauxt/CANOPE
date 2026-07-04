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
>
> **Update 2:** after actually re-running against real data, two more
> issues turned up — the report still wasn't always found (loose-sourced
> scripts run from a different working directory than the one containing
> `CANOPE_report.Rmd`), and the z-score panel's spread was still much wider
> than the coverage panel warranted, even after Round 2's scale fix. See
> "Round 3" below.
>
> **Update 3:** found the actual root cause of the report never generating
> — a well-known `rmarkdown::render()` gotcha (relative output paths
> resolve against the *Rmd's* directory, not the caller's) that Round 2/3's
> template-*finding* fixes never addressed, since the template was being
> found fine; it was the *output* path that broke. This round was verified
> with a full, real `run_canope()` → `generate_plots()` →
> `generate_canope_report()` run on synthetic data (R actually installed
> and executed this time, not just parsed) with an injected deletion, and
> the report now renders end-to-end. See "Round 4" below.
>
> **Update 4:** with the report finally rendering against real data, one
> more report bug turned up (a raw-HTML-inside-`results='asis'` chunk being
> silently mangled by Pandoc — fixed and verified) and a serious *statistical
> calibration* question was raised about the 95% predictive interval /
> z-score panels. That one turned out to be more interesting than a simple
> bug: a rigorous check against real report data found the interval
> genuinely under-covers in aggregate, but the effect is heavily
> concentrated in one sample's calls across several genes, in a pattern
> that looks much more like real biological signal than a statistical bug
> — and a Monte Carlo check of the obvious "fix" showed it would actually
> *overcorrect* the normal case. See "Round 5" below; this one ends with a
> question rather than a shipped fix, on purpose.
>
> **Update 5:** turned the Round 5 investigation into a permanent,
> automatic diagnostic (`check_background_calibration()`) instead of a
> one-off manual check, since the person asking couldn't say for certain
> whether the sample in question was a control. Worth flagging honestly:
> getting this wired in took several iterations because of a self-inflicted
> testing bug on my end (a single-line grep test gave false "not found"
> results against Pandoc's word-wrapped HTML output), which briefly led to
> shipping an unnecessary workaround for a rendering bug that never
> actually existed. Caught it, reverted the unnecessary part, kept the
> actual feature. See "Round 5" (updated) for the honest version of events.

## Files

| File | Status |
|---|---|
| `hmm_engine.R` | **Unchanged** — the engine core, as requested |
| `call_cnvs.R` | Bug-fixed (round 1 + round 2) |
| `data_utils.R` | Bug-fixed |
| `generate_plots.R` | Bug-fixed (round 1 + round 2 + round 3 — z-score/CI panels) |
| `canope_report.R` | Bug-fixed (round 2 + round 3 + round 4 — template path + render path) |
| `CANOPE_report.Rmd` | Bug-fixed (round 2 + round 3 + round 5 — CI/z-score fixes + raw-HTML-in-asis fix) + new background-calibration flag |
| `run_canope.R` | Bug-fixed (round 1 + round 2) |
| `qc_reference_utils.R` | Minor robustness only, logic untouched |
| `run_canope.R` | Bug-fixed + orchestrates all new features |
| `canope_utils.R` | **New** — shared helpers (PCA plot, exon numbering, parsers, `check_background_calibration`...) |
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

## Round 3 bug fixes

Found by actually re-running the fixed pipeline against real data and
looking at the resulting plots. The ratio panel's range did narrow up
nicely (confirming Round 2's fix worked for that panel), but the report
still sometimes wasn't produced, and the z-score panel still looked far
more dramatic than the coverage panel above it warranted.

### 3. Report still not always found (`canope_report.R`)

Round 2 fixed the *installed-package* lookup path (`rmarkdown/` →
`rmd/`), but the README's quick-start usage doesn't install CANOPE as a
package — it `source()`s the loose `.R` files directly, which means
`system.file(package = "CANOPE")` never resolves to anything (empty
string) and was never actually the fallback in play for that workflow. The
only fallback that *could* match was `here_candidate <-
file.path(getwd(), "CANOPE_report.Rmd")` — which only works if
`run_canope()` happens to be called with the working directory set to
the exact folder the scripts were sourced from. Call it from a wrapper
script, an RStudio project rooted elsewhere, or a scheduled job, and
`getwd()` no longer points there — so the template silently couldn't be
found even sitting right next to `canope_report.R`.

**Fix:** `canope_report.R` now captures its own source-time location (the
same "what folder was I sourced from" trick `here`/`this.path`-style
packages use, implemented locally to avoid adding a dependency) and adds
that folder as another fallback, independent of `getwd()`. Verified this
resolves correctly even when `setwd()` points somewhere else entirely
before `run_canope()` is called. The failure message (if it still can't be
found) now also lists every path it actually looked in, so a genuine
"file isn't there" case is immediately diagnosable from the log instead of
just "could not locate."

Since `generate_canope_report()` is called from inside a `tryCatch(...,
error = function(e) log_msg("WARNING", ...))` in `run_canope()` (same
graceful-degradation pattern as the QC/PCA/VCF steps), a failure here
doesn't stop the pipeline — it logs a `[WARNING] Report generation failed:
...` line (both to console via `message()` and to the log file) and moves
on. If the report still doesn't appear, check that line in your log for
the specific reason (now including the full candidate-path list for a
missing-template failure) — or pass `report_template` explicitly to
`run_canope()`/`generate_canope_report()` to bypass the search entirely.

### 4. Z-score panel still far more dramatic than the coverage panel

Round 2 fixed the *depth-scale* mismatch (raw counts vs. the internally
normalised scale `mean`/`var_estimate` were fit on) — visible in the new
PDFs as the ratio panel's range narrowing to something sane. But the
z-score panel's own SD estimate had a second, independent problem: it was
the raw sample standard deviation computed across just the handful of
reference columns for that one target (`apply(ref_mat, 1, sd)`), floored
at an essentially-meaningless `1e-6`. With only a handful of references,
that per-target SD is a noisy, low-degrees-of-freedom estimate on its
own — at any target where the references happened, by chance, to cluster
unusually tightly, even a trivial absolute deviation from the test sample
produced an enormous, spurious z-score.

Confirmed with a simulation: 5 reference samples with ordinary Poisson
noise around a shared mean, with two targets deliberately given very tight
by-chance clustering across those 5 refs, and a test sample with **no
real CNV** — just ordinary noise. The old per-target sample SD produced
z-scores of **21.75 and 67.88** at those two targets (a dramatic, fully
spurious "signal") while every other target sat under 1.5; the new
approach keeps all targets in the same sane ~0–1.5 range regardless of
reference-sample luck.

**Fix:** the z-score panel's SD is now `sqrt(var_estimate)` — the same
robust, MAD-based variance (floored at `mean + 1`, not `1e-6`) that
`call_cnvs()` already computes for HMM calling and that the ratio/NB-
interval panel already uses — instead of each panel inventing its own,
much noisier variance estimate from a handful of samples. This also
removes the separate single-reference-sample MAD fallback that was in
`create_zscore_plot()`, since `var_estimate` is already a well-behaved,
appropriately-floored estimate regardless of how many references are
available. Applied identically in `generate_plots.R` and the duplicated
logic in `CANOPE_report.Rmd`.

## Round 4 bug fixes

Found this round by actually *executing* the pipeline end-to-end — R,
pandoc, and every dependency (`rmarkdown`, `ggplot2`, `dplyr`, `DT`,
`plotly`, `knitr`, `tidyr`, `htmltools`, `kableExtra`, `nnls`,
`matrixStats`, `data.table`) installed for real this round, run against a
synthetic dataset (10 samples, 20 targets across two genes, a deliberately
injected heterozygous deletion) through `run_canope()` in full, including
`generate_canope_report()`. Round 1–3 fixes were all confirmed *still
correct* by this run; one further bug turned up, and it's the actual
reason the report kept failing to generate even after Round 3's
template-path fixes.

### 5. The real reason the report wasn't generating

Round 2/3 fixed *finding* `CANOPE_report.Rmd` — but the report was still
failing, with an error like:

```
Report generation failed: The directory 'report' does not exist.
```

...pointing at a directory that, from wherever `run_canope()` was actually
called, unquestionably *did* exist (it had just been created a few lines
earlier by `dir.create()`). The template was being found fine by this
point; it was the *output* path that broke, for a completely different and
much less obvious reason: `rmarkdown::render()` does not resolve a
relative `output_file` against the caller's working directory — it
resolves it against the directory containing the *input* Rmd file. Once
Round 2/3 got `template_path` to correctly resolve to an absolute path
pointing at wherever `CANOPE_report.Rmd` actually lives (typically
alongside the other pipeline `.R` files, not alongside your output
folder), every relative `output_file` — `file.path(output_dir,
"CANOPE_..._report.html")`, where `output_dir` is itself built relative to
your `output_file`/`rdata_output` paths — got resolved against *that*
folder instead, which almost never contains a matching `report/`
subdirectory.

This is a known, if obscure, R Markdown gotcha, confirmed here with a
minimal, isolated reproduction (an Rmd in one directory, `render()` called
from an unrelated working directory with a relative `output_file`) that
throws exactly the same "directory does not exist" error regardless of
whether that directory exists relative to `getwd()`. It's also *why* Round
2/3's independent, carefully-tested template-path fix never actually
solved the "report not generated" complaint — that fix was correct, just
solving a different half of the problem than the one actually causing the
failure.

**Fix:** `generate_canope_report()` now builds `output_file` as an
absolute path — `output_dir` is already guaranteed to exist by this point
(created a few lines above), so `normalizePath(output_dir, mustWork =
TRUE)` is safe, and only the (not-yet-existing) filename is appended
afterward. (Naively calling `normalizePath()` on the *full* file path,
including the nonexistent filename, silently fails to absolutize anything
in R — confirmed directly — so it has to be the directory that gets
normalized, not the full path.)

While tracking this down, a second, quieter instance of the same root
cause turned up: `rmarkdown::render()` also evaluates the Rmd's own code
chunks with the working directory set to the Rmd's directory by default,
not the caller's. The report's "Pipeline Log" section reads
`params$log_file` (built by `run_canope()` relative to *its* caller's
working directory) behind an existing `file.exists()` guard, so this one
wasn't fatal — it would have just silently rendered without the log/
warnings tab, with nothing indicating why. Fixed by passing `knit_root_dir
= caller_wd` (captured at the top of `generate_canope_report()`, before
anything else could change it) to `rmarkdown::render()`, so every relative
path touched inside the Rmd resolves the same way its caller intended.

### Re-confirmation of Round 2/3's z-score/CI fix, on a real pipeline run

With the report now actually rendering, this was a chance to check Round
2/3's z-score/CI fix against a *real* `call_cnvs()` → `models` → plotting
round-trip rather than an isolated simulation. Injected a clean
heterozygous deletion (~0.5×) into one sample across a few targets, left
everything else as ordinary Poisson noise, and checked every target's
ratio and z-score, not just the deleted ones:

| target | region | log2 ratio | z-score |
|---|---|---|---|
| 5, 6, 7 | **deleted** | −0.67 to −0.89 | **−6.2 to −16.5** |
| every other target (15 of them) | normal | 0.00 to 0.18 | −0.2 to 5.4 |

Clean separation, no spurious blowups at any non-CNV target — consistent
with (and considerably more convincing than) the earlier synthetic-only
simulation, since this went through the actual `call_cnvs()` HMM-adjacent
code path end to end. Also directly compared the static-PDF code path
(`generate_plots.R`) against the interactive-report code path
(`CANOPE_report.Rmd`) on the same data — both independently computed
identical z-scores (−16.47, −6.18, −14.94 for the three deleted targets),
confirming the two duplicated implementations haven't drifted apart.

If your data still shows the ratio/z-score panels looking wrong after
this fix, it'd help a lot to see the regenerated PDFs/HTML (or the
relevant `models[[sample]]` values) — everything reproducible in this
environment now checks out, so a further discrepancy likely means
something specific to your actual data/setup that isn't captured by the
synthetic test above.

## Round 5: report HTML bug (fixed) + a real look at CI calibration (not a quick fix)

This round had a real report and real report data to work from for the
first time, which changed things — one bug was found and definitively
fixed, and the CI/z-score question turned into a genuine statistical
investigation rather than a quick patch.

### 6. Part of the report showing as literal, unrendered code

The "QC Metrics Overview" section's missing-chromosome warnings (chrX/chrY
not in the BED) were showing up as raw, HTML-escaped text in a code block
instead of as the styled warning boxes they're supposed to be. Root cause:
in a `results='asis'` chunk, `cat()` output goes straight into the
Markdown source that Pandoc then parses — and Pandoc's Markdown reader
treats any line indented 4+ spaces as an **indented code block**, not raw
HTML. The offending `cat(sprintf('\n      <div class="alert...'))` call
had each HTML line indented to match the surrounding R code (natural,
readable R style) — which is exactly the pattern that trips this up.

Confirmed with an isolated reproduction: identical HTML content rendered
as literal escaped text when indented, and as real, working HTML when
not — and cross-checked against the `confidence_recap` chunk elsewhere in
the *same* file, which already used single-line, non-indented `cat()`
calls and rendered correctly in your actual uploaded report. Fixed by
switching the missing-chromosome warning block to that same proven style
(one `cat()` call per line, starting at column 0), then re-verified with a
full pipeline run — the alert boxes now render as actual `<div>`/`<h4>`/
`<p>` elements, not escaped text.

### The 95% CI / z-score question

Checked this properly using the actual plotly data embedded in your
uploaded HTML (exact numbers, not pixel-eyeballing a PDF): for each of the
6 CNV calls, compared the test sample's own ratio at every *non-called*
("background") exon in the plotted window against that exon's modelled
95% predictive interval. A well-calibrated interval should have roughly
5% of these fall outside by chance.

| Sample | Gene | Background exons | % outside 95% interval |
|---|---|---|---|
| HORIZON | MSH2 | 10 | 20% |
| HORIZON | PMS2 | 8 | 75% |
| HORIZON | BRAF | 10 | 100% |
| HORIZON | BRAF | 10 | 100% |
| HORIZON | PRSS1 | 5 | 100% |
| SGT2600617 | PTEN | 10 | **0%** |

Aggregate: **33/53 (62%) outside** — a real, substantial finding, not
noise. But look at the pattern: it's not spread evenly, it's *entirely*
one sample (HORIZON), across every one of its genes, while the other
sample's call is perfectly calibrated. And within HORIZON's calls, the
"background" points aren't scattered randomly above and below zero the
way pure noise would look — they're almost all positive (elevated),
sometimes across the *entire* plotted window, not just the specific exons
that got called as the CNV.

That pattern — one-sided, near-uniform elevation across a whole gene
window, concentrated in one specific sample — looks much more like **real
broad copy-number signal in that sample** than like a statistical
under-estimation of variance. "HORIZON" strongly suggests a commercial
reference/control material (Horizon Discovery is a well-known maker of
characterized CNV reference standards); if so, this is expected biology,
not a bug — though it might mean the *called* CNV boundary is narrower
than the true extent of the alteration, which is a separate question from
interval width.

**The obvious fix doesn't hold up.** The natural hypothesis — that the
MAD-based variance estimator is biased low with a small number of
reference samples (10-12 here), and needs a small-sample correction
(e.g. the standard Leys et al. 2013 factor for MAD) — is theoretically
reasonable, so it was tested properly before touching anything: a Monte
Carlo simulation (3,000 trials, `n_ref = 12`, pure Poisson noise, no real
CNV) showed the *existing* estimator already lands almost exactly on
target (5.4% vs. a 5% goal), and the "corrected" version *overshoots* to
1.0% — i.e. it would make the interval needlessly wide for well-behaved
cases, which would reduce sensitivity to real CNVs elsewhere. That
correction was implemented, tested, and then deliberately **reverted**
rather than shipped once the simulation contradicted the theory — a
reminder that this class of fix needs to be checked against a null
simulation before shipping, not just justified by a citation.

**No code change shipped for this one.** Given the evidence points at
"real signal in a specific sample" rather than "the interval is generically
wrong," the honest next step isn't a quick patch — it's confirming what
HORIZON actually is. Depending on the answer:
- **If HORIZON is a known reference-standard/control sample** with
  expected alterations in these genes: this isn't a bug. Worth checking,
  separately, whether the *called* CNV boundaries fully capture the known
  extent of each control alteration (a Viterbi/segmentation question,
  which lives inside `hmm_engine.R` and is out of scope for a "fix around
  the HMM" change) — and whether BRAF/PRSS1 belong in
  `low_confidence_genes` alongside PMS2 (already there, correctly, for its
  known pseudogene homology) given how consistently they show this same
  pattern.
- **If HORIZON is an ordinary patient sample** with no expected CNVs in
  MSH2/PMS2/BRAF/PRSS1: something more specific is going on — most likely
  a batch/technical effect between HORIZON's sequencing run and its chosen
  references (something reference-panel-based callers are generally
  vulnerable to, since the variance estimate only sees inter-reference
  spread, never a test-vs-reference technical difference) — which would
  need a real per-sample or per-batch dispersion term, not a constant
  tweak, and is worth scoping properly rather than guessing at.

### New: automatic background-calibration flag (in lieu of a definitive answer)

The question above ("is HORIZON a control sample?") didn't get a definite
answer. Rather than leave the investigation as a one-off manual check
(parsing plotly JSON out of an HTML file by hand), the same check is now a
permanent, automatic part of every call, in both the static PDFs and the
interactive report:

- **`check_background_calibration()`** (new, in `canope_utils.R`, shared by
  both plotting paths): for a call's plotted window, tests whether the
  fraction of *non-called* exons falling outside the modelled 95% interval
  is statistically higher than the ~5% a well-calibrated interval implies
  (a one-sided binomial test against a 5% null, requires at least 5
  background exons to test at all — below that the percentage alone is too
  noisy to mean anything).
- When it flags, the **static PDF's** ratio-panel subtitle gets an inline
  note (`95% NB predictive interval [33% of background exons outside CI
  (3/9) — check region/reference match]`), and the **interactive report**
  gets a prose callout above the plots explaining what it can mean (signal
  beyond the called boundary, a weak reference match, or a technical/batch
  effect) and that it's a prompt for a manual look, not a correction —
  it doesn't touch the call, the interval, or the confidence score.

This doesn't resolve the HORIZON question — it makes the same check
available for *every* call, automatically, so a judgment like the one made
manually in the Round 5 write-up above doesn't require re-deriving it by
hand each time.

**Honest note on how this went:** wiring this in took several more
iterations than it should have, because of a self-inflicted bug in my own
*testing*, not the feature. After adding the display logic, a `grep -c`
check against the rendered HTML kept coming back "0 matches," which looked
exactly like a rendering bug — so time was spent chasing one: reproducing
it, restructuring the `cat()` calls to "fix" it, reproducing it *again*
after the fix appeared to not work either. It eventually turned out the
feature had been working correctly the entire time; Pandoc had word-wrapped
the flagged sentence across a line break in the HTML *source* (irrelevant
to how a browser renders it), and a single-line `grep` simply can't match
text split across two lines. Confirmed with a proper multi-line-aware check
(Python, newlines stripped before searching) that even the very first
version worked. The unnecessary "fix" from the false alarm was reverted;
what's shipped is the original, straightforward implementation, now
actually verified with a test that isn't the thing that was broken.

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

**Round 3:** the source-directory detection trick in `canope_report.R`
(`.canope_report_source_dir`) was actually exercised — sourced the file
from `/home/claude/canope` while `getwd()` was set to an unrelated `/tmp`
directory, and confirmed it still resolves to the correct folder. The
z-score SD fix was reproduced with a targeted simulation (5 references, two
targets given deliberately tight by-chance clustering, no real CNV): old
per-target sample SD gave spurious z-scores of 21.75 and 67.88 at those two
targets against <1.5 everywhere else; the `var_estimate`-based fix keeps
every target in the same sane range. All R files and Rmd chunks were
re-parsed cleanly after these changes (R itself is available in this
environment, unlike the earlier rounds). Still not exercised: an actual
`rmarkdown::render()` of the report, or a full BAM-to-plot pipeline run —
same caveat as Round 2.

**Round 4:** by far the most thorough round — R, pandoc, and the full
dependency set (`rmarkdown`, `ggplot2`, `dplyr`, `DT`, `plotly`, `knitr`,
`tidyr`, `htmltools`, `kableExtra`, `nnls`, `matrixStats`, `data.table`)
were actually installed, and a synthetic 10-sample/20-target dataset with
an injected heterozygous deletion was run through real `run_canope()` →
`generate_plots()` → `generate_canope_report()` calls — not simulations or
isolated unit tests of extracted logic, the actual functions as shipped.
This is what caught bug #5 (the `rmarkdown::render()` relative-path
issue): it never showed up in Round 2/3's unit-level testing because
nothing before this round had actually called `rmarkdown::render()` with
realistic caller/template/output directory layouts. The fix was isolated
and confirmed with a minimal reproduction outside the full pipeline too
(an Rmd in one directory, `render()` invoked from an unrelated working
directory), then re-confirmed with the full pipeline run producing an
actual 5.9MB self-contained HTML report end to end. The Round 2/3
z-score/CI fix was independently re-verified against this same real run
(table in the Round 4 write-up above) and against both the
`generate_plots.R` and `CANOPE_report.Rmd` code paths directly, confirming
they still agree with each other.

**Round 5:** the qc_overview fix was verified against a real render, not
just reasoned about — reproduced the exact bug in isolation first
(indented vs. unindented HTML in a `results='asis'` chunk), then confirmed
the fix by re-running the full synthetic pipeline and grepping the actual
output HTML for the previously-broken markup, which now appears as real
`<div>`/`<h4>`/`<p>` tags instead of escaped text in a `<pre><code>`
block. The CI-calibration investigation was done directly against the
user's real uploaded HTML report: the plotly JSON embedded in each of the
6 CNV-call widgets was parsed programmatically (not eyeballed) to extract
the exact ratio values and interval bounds at every background exon,
giving the 33/53 (62%) aggregate figure and the per-call breakdown in the
Round 5 write-up above. The candidate fix (small-sample MAD correction)
was Monte Carlo tested (3,000 trials) before being proposed — and that
test is what caught it overcorrecting, so it was reverted rather than
shipped. No code was changed for the CI-calibration question this round;
resolving it properly needs to know whether HORIZON is a reference
standard or a patient sample, since the right fix (if any) is different
in each case.

**Round 5 (continued):** `check_background_calibration()` was unit-tested
directly against the exact real numbers extracted from the report (9
background exons, 3 outside → confirmed `flag = TRUE`, matching the
one-sided binomial calculation by hand). End-to-end verification of the
new report callout went through a detour worth recording honestly: an
initial `grep -c` check kept reporting the flag text as absent, which
prompted real (unnecessary) debugging — reproducing the "failure,"
restructuring the `cat()` calls, reproducing it again. The actual cause
was the test, not the code: Pandoc had word-wrapped the sentence across a
line break in the HTML source, which a single-line `grep` can't match
across. A proper check (Python, newlines stripped first) confirmed the
original, simplest implementation was correct from the start. Final
verification used that corrected method: the flag renders correctly in
both the static PDF subtitle and the interactive report's prose callout,
the two stay consistent with each other, and the pre-existing qc_overview
HTML fix from earlier in this round was re-confirmed intact in the same
run.
