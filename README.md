<div align="center">
  <img src="assets/logo.png" alt="CANOPE Logo" width="250"/>

  # CANOPE
  **Copy-number Analysis using Normalized Observation Profiling for Exomes**

  [![R Build Status](https://github.com/TL/CANOPE/workflows/R-CMD-check/badge.svg)](https://github.com/TL/CANOPE/actions)
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
  [![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-brightgreen.svg)](https://cran.r-project.org/)
</div>

---

## Overview

**CANOPE** is an R package designed for the highly accurate detection of Copy Number Variations (CNVs) in targeted panel and whole-exome sequencing (WES) data. 

By employing low-variance reference sample selection, Negative Binomial emission modeling, and a Viterbi-decoded Hidden Markov Model (HMM), CANOPE robustly distinguishes true deletions and duplications from baseline sequencing noise.

## Key Features

* **YAML-Driven Pipeline:** A single `canope("config.yaml")` command handles BED preprocessing, BAM coverage extraction, GC correction, QC, HMM calling, confidence scoring, VCF export, and HTML reporting.
* **Smart Normalization:** Automatically selects the most highly correlated reference samples (up to a user-defined threshold) to build an expected baseline via Non-Negative Least Squares (NNLS).
* **Robust HMM Engine:** Implements a Hidden Markov Model with dynamically computed transition probabilities (based on genomic distance) and Negative Binomial emission probabilities.
* **Auto-Extraction:** Built-in auto-extraction of per-target read counts directly from BAM files (via `GenomicAlignments` or `megadepth`) and GC content directly from `BSgenome` packages or indexed FASTA files.
* **Interactive, Multi-Sample HTML Report:** one standalone `CANOPE_<prefix>_report.html` (no server, no Shiny app, works offline) with sample/call selectors, sortable and filterable tables, linked plotly panels for every call, PCA / correlation heatmaps and a light/dark theme. It handles any number of samples and renders even when nothing was called.
* **Warnings Everywhere It Matters:** QC thresholds are checked per sample, per chromosome and per call, and every breach is flagged in the metrics tables *and* collected on a **Warnings** tab. In particular, a run that did **not** analyse both autosomes and gonosomes says so in red on the Overview and Chromosomes tabs (see [Warnings](#warnings-in-the-report)).
* **Native Visualization:** static PDF plots (tracking plots, coverage ratios, gene annotations) for every call, independent of heavy third-party plotting frameworks.
* **Terminal-Exon Padding:** `pad_terminal_exons` (bp, default `0` = off)
  extends the outward-facing edge of each gene's first and last exon, where
  there's no neighbouring exon to carry the signal if coverage thins out
  right at the boundary. Padding is applied "if possible" — it never
  overlaps a neighbouring target or crosses a contig boundary, clamping
  short instead. Ported from ECHO.
* **Gene-Boundary Spacing in Plots:** `plot_gene_gap` (default `1`) inserts
  blank x-axis space between a gene's last exon and the next gene's first
  exon in every PDF and HTML-report panel, so a window that spans more
  than one gene doesn't read as one continuous, unbroken feature. Ported
  from ECHO.
* **Off-Target/Filler Region Handling:** BED panels sometimes include
  normalization "backbone" probes that are never a real gene exon (often
  named things like `HorsROI`). Rows matching `off_target_pattern`
  (default `"^HorsROI"`, configured under `bed_process_args`) are handled
  per `off_target_handling`: `"na"` (default), `"remove"`, or `"merge"`.
  Ported from ECHO.

## Installation

You can install the development version of CANOPE directly from GitHub using `devtools`.

```R
# Install devtools if you haven't already
if (!require("devtools")) install.packages("devtools")

# Install CANOPE
devtools::install_github("lavauxt/CANOPE")
```

### Dependencies
Core: `dplyr`, `nnls`, `stats`, `utils`, `tools`, `data.table`, `matrixStats`,
`yaml`, `HMM`, `ggplot2`, `jsonlite`, and the Bioconductor packages `Biostrings`, `rtracklayer`,
`GenomicRanges`, `GenomeInfoDb`, `BiocGenerics`, `Rsamtools`,
`GenomicAlignments`, `SummarizedExperiment`, `IRanges`, `S4Vectors`.

Optional (only needed for specific features): `megadepth` (fast BAM
coverage), `GenomicFeatures` + a `TxDb.Hsapiens.UCSC.*` package + `org.Hs.eg.db`
(BED `REGEN` mode), `ggrepel` (nicer PCA labels), and `plotly` + `htmltools`
(only used to find the bundled plotly.js so it can be inlined in the HTML
report; without them the report loads plotly.js from a CDN instead).
`rmarkdown`, `knitr`, `DT`, `kableExtra` and `tidyr` are **no longer needed**.

## Quick Start

### 1. Data Preparation
Ensure you have your BAM files and a BED file of target regions. If you don't have pre-computed count matrices, CANOPE will automatically extract coverage from your BAMs and compute GC content using a `BSgenome` package or an indexed FASTA file.

### 2. Create a Configuration File
Create a `config.yaml` file to define your inputs, outputs, and HMM parameters.

```yaml
input:
  bed: "data/targets.bed"
  bamdir: "data/bams/"
  # sample_table: "data/samples.tsv"           # sample_name + gender (M/F): PCA colours and report sex check
  # reads_file / gc_file: pre-computed counts and GC (see "Input File Formats")

settings:
  bsgenome_pkg: "BSgenome.Hsapiens.UCSC.hg19"  # Or use fasta_file instead
  modechrom: "A"                               # "A" | "XX" | "XY": one run analyses ONE of them (see Warnings)
  p_value: 1e-08
  Tnum: 6
  numrefs: 30
  coverage_backend: "bioconductor"             # or "megadepth"
  bed_process: "STANDARD"                      # or "NO" / "REGEN"
  bed_process_args:
    off_target_pattern: "^HorsROI"             # e.g. filter out backbone/filler probes
    off_target_handling: "na"                  # "na" | "remove" | "merge"
  pad_terminal_exons: 0                        # e.g. 10 to pad each gene's first/last exon by 10bp
  plot_gene_gap: 1                             # blank x-axis space between adjacent genes in plots
  run_qc_metrics: TRUE
  export_vcf: TRUE
  report: TRUE

output:
  dir: "results"
  prefix: "CANOPE"
```

### 3. Run the Pipeline
Execute the entire pipeline by passing the config file to the `canope()` wrapper:

```R
library(CANOPE)

# Run the full end-to-end pipeline
canope(config_path = "config.yaml")
```

### 4. Visualizing Results
If `report: TRUE` is set (the default), `CANOPE_<prefix>_report.html` is written to the output directory: open it in any browser. Static PDF plots for each CNV call are also saved in `results/plots/`.

The configuration can be flat (as in `inst/config.yaml`) or nested under `input:` / `settings:` / `output:` (as above); every option is found in either layout.

## The Interactive Report

The report is a single self-contained HTML file. It is built from a JSON blob that R embeds in the page; the page itself (tabs, tables, plots) is plain JavaScript with plotly.js inlined, so it can be e-mailed, archived or opened offline. No Shiny app or R session is needed to read it.

| Tab | What it shows |
| :--- | :--- |
| **Overview** | Run status tiles, all red (fail) warnings, which chromosome classes were analysed, calls per sample, sample health |
| **Warnings** | Every warning, filterable by level (fail / warn / note), category, sample and free text; CSV download |
| **Samples** | One row per sample (multi-sample): depth, total reads, best correlation, noise, number of references, calls, inferred vs annotated sex, flags. Cells that breach a threshold are coloured; click a row for its calls |
| **QC metrics** | Bar charts of depth / reads / correlation / noise with the threshold line, plus the sample- and exon-level QC tables |
| **Chromosomes** | Analysed / not analysed per chromosome class, per-chromosome depth relative to the autosomes, sex inference |
| **PCA & similarity** | PCA of the coverage profiles (up to 5 PCs, coloured by sex, warning level or number of calls) and a clustered sample-to-sample Spearman heatmap |
| **CNV calls** | Sortable / filterable table of all calls (sample, type, gene, Q_SOME, confidence, ...) with a genome-wide map; CSV export of the filtered rows |
| **Call viewer** | For the selected sample and call: coverage vs. its reference samples, log2(observed/expected) with the negative-binomial interval, and z-score; the three panels zoom together. Hidden when there are no calls |
| **Method & thresholds** | The parameters and thresholds actually used, and a glossary |
| **Run & log** | Run information, session info, pipeline log |

Interactive plots are embedded for up to `max_call_plots` calls (an argument of `generate_canope_report()`, default 300) to keep the file small; the table lists every call and the PDF plots are unaffected.

You can rebuild the report from an existing workspace without re-running anything:

```R
generate_canope_report(
  rdata_output    = "results/CANOPE_TEST_workspace.RData",
  qc_metrics_file = "results/CANOPE_TEST_QC_metrics.tsv",
  output_dir      = "results",
  prefix          = "TEST",
  sample_table    = "data/samples.tsv",     # optional
  plotly_source   = "inline"                # "cdn" gives a ~120 kB file that needs internet access
)
```

Workspaces written by older CANOPE versions (without `run_info`) still render; the sections that need the missing information say "Not available for this workspace".

### Warnings in the report

Every warning has a level (**fail** = red, results may be incomplete or wrong; **warn** = amber, interpret with care; **info**) and appears next to the metric that triggered it (coloured cells in the Samples and QC tabs) and on the Warnings tab.

| Category | Raised when |
| :--- | :--- |
| **Chromosomes** | **fail**: autosomes and gonosomes were **not both analysed** (`modechrom = "A"` leaves chrX/chrY out, `"XX"`/`"XY"` leave the autosomes out, and chrY is dropped unless `removeY = FALSE`) -- "no call" on those chromosomes does not mean "no CNV"; re-run with the other `modechrom` to cover the rest. **warn**: a chromosome class is absent from the BED, or was in the BED but not analysed (with the reason) |
| **Samples** | samples that **failed** calling (fail) or were **excluded** before calling, e.g. noisy outliers (warn); fewer than 4 samples analysed (too few references) |
| **Sample metrics** | per-sample flags: median depth < `qc_min_cov`; total reads < `qc_min_total_reads`; best Spearman correlation < `qc_min_corr`; noise robust z > `qc_zscore`; fewer than 5 references; genotyping skipped because of too many calls (every call LOW); unusually many calls compared with the cohort; runs of >= 3 targets with ~no reads removed before calling (`min_sample_count`), where a homozygous deletion cannot be called; annotated sex (`sample_table`) different from the sex inferred from chrX/chrY depth, or a sample whose sex does not fit the `modechrom` used |
| **Calls** | calls whose non-called neighbouring exons are not calibrated around the call (background flag); calls without an embedded plot (plot limit). Calls in homology-prone genes are always scored LOW |
| **QC metrics** | exons flagged by exon-level QC (high MAD, low mean, extreme GC) |
| **Method / Pipeline log** | informational notes (legacy engine used, warning lines found in the pipeline log) |

## Input File Formats

**BED File**
Standard 0-based BED file. If `bed_process` is set to `"STANDARD"` or `"REGEN"`, the input BED will be automatically annotated with gene/exon information before coverage extraction.

**Reads File (`canope.reads.tsv`)** *(Optional)*
If you already have pre-computed read counts, you can provide a tab-separated file instead of letting CANOPE extract from BAMs. It must contain exactly four metadata columns (`chromosome`, `start`, `end`, `GENE`) followed by one column of raw read counts per sample -- `run_canope()` matches columns by position after the fourth column, so an extra or missing metadata column will shift every sample's counts by one.

| chromosome | start | end | GENE | Sample1 | Sample2 | Sample3 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| chr1 | 10000 | 10150 | BRCA1 | 45 | 50 | 48 |
| chr1 | 10500 | 10620 | BRCA1 | 12 | 15 | 11 |

**GC file** *(optional, `gc_file`)*
`chromosome`, `start`, `end`, `GENE`, `GC_CONTENT` (fraction 0-1 or percent). Rows are matched to the targets **by coordinates**, so the order does not matter; if no coordinate matches but they all do after a +1 shift on `start` (0-based BED vs 1-based GC table), the shift is applied and reported.

**Sample table** *(optional, `sample_table`)*
Tab-separated, with `sample_name` and `gender` (`M`/`F`, also `male`/`female`/`XY`/`XX`). Used to colour the PCA and, in the report, to compare the annotated sex with the sex inferred from chrX/chrY coverage.

## Output Files

For `output_dir: "results"` and `output_prefix: "TEST"`:

| File | Content |
| :--- | :--- |
| `CANOPE_TEST_CNVCall.csv` | all calls (header only when nothing was called) |
| `CANOPE_TEST_calls.vcf` | multi-sample VCF (`vcf_per_sample: true` also writes one per sample) |
| `CANOPE_TEST_QC_metrics.tsv` | sample / exon / chromosome QC flags |
| `CANOPE_TEST_report.html` | the interactive report |
| `CANOPE_TEST_workspace.RData` | counts, models, references, calls and `run_info` (run summary used by the report) |
| `CANOPE_TEST_PCA.pdf`, `plots/` | static PCA and per-call PDF plots |
| `CANOPE_TEST_pipeline.log` | log |

## Output Format (`CNVCall.csv`)

The output is a tabular file containing the detected variations:

* `CHR`, `INTERVAL`, `MID_BP`: Genomic coordinates of the CNV.
* `CNV`: Type of variant (`DEL` for deletion, `DUP` for duplication).
* `SAMPLE`: The sample in which the event was detected.
* `KB`: Size of the variation in kilobases.
* `TARGETS`: The target ID range spanning the event (e.g. `42..47`).
* `NUM_TARG`: Total number of contiguous targets involved.
* `GENE`: Gene name(s) overlapping the call.
* `MLCN`: Most Likely Copy Number (e.g., 1 for heterozygous deletion, 3 for duplication); `NA` when it could not be estimated.
* `Q_SOME`: Phred-scaled quality score of the CNV call.
* `NUM_REFS`, `REF_SAMPLES`: Number and names of reference samples used.
* `Confidence`, `CN_label`: Added when `score_confidence = TRUE` (default).

## Additional Features

- **BED preprocessing** (`canope_bed_process.R`) — `STANDARD`/`REGEN`/`NO`
  modes, gene/exon annotation from panel files or RefSeq. Exon-number
  extraction from the BED name (e.g. `..._ex12_...`) requires a genuine
  `ex`/`exon` token boundary (not a substring of a larger word) and takes
  the *last* such token in the name rather than the first, so compound
  names like `5utr1-ex1` or `utr2-ex10-extra` resolve to the right exon.
  Off-target/filler intervals (e.g. `HorsROI` backbone probes) are caught
  by name via `off_target_pattern` before any exon numbering happens —
  see `handle_off_target_regions()` below. Ported from ECHO.
- **Off-target/filler region handling** (`canope_utils.R::handle_off_target_regions()`) —
  rows whose parsed gene matches `off_target_pattern` (default
  `"^HorsROI"`, set under `bed_process_args`) are, per
  `off_target_handling`: kept with `GENE = NA` and excluded from numbering
  (`"na"`, default), dropped (`"remove"`), or attached to the nearest
  neighbouring real gene on the same chromosome and numbered as one of its
  exons (`"merge"`). Without this, such rows get numbered as if the filler
  label were a gene in its own right, and show up interleaved with the
  real gene's exons in plots. Ported from ECHO.
- **Terminal-exon padding** (`canope_utils.R::pad_gene_terminal_exons()` /
  `pad_bed_file()`) — extends the outward edge of each gene's first/last
  exon by `pad_terminal_exons` bp (default `0`, off) before coverage/GC
  extraction, so a thin sliver of low coverage right at a gene boundary is
  less likely to drag the whole exon's count down. Applied "if possible":
  clamped short of the requested amount rather than overlapping a
  neighbouring target or crossing a contig boundary. Internal exon
  boundaries are left alone. Only takes effect when extracting fresh from
  BAMs (a supplied `reads_file` is left untouched, since it was already
  extracted over some fixed window). Ported from ECHO.
- **Gene-boundary spacing in plots** (`canope_utils.R::compute_gene_gap_positions()`) —
  shared by both the PDF panels (`canope_generate_plots.R`) and the HTML
  report (`canope_report.R`), so a CNV window spanning more than one
  gene shows real blank space, in every panel, exactly between a gene's
  last exon and the next gene's first exon — lines and ribbons break there
  rather than visually bridging the gap. Controlled by `plot_gene_gap`
  (default `1`; `0` disables it). Ported from ECHO.
- **Sample/exon QC exclusion** (`canope_qc_reference_utils.R`) —
  `sample_qc` drops outlier samples (robust z-score on cross-target noise)
  and `exon_qc` drops problematic exons (high cross-sample MAD, low mean
  coverage, or GC content outside `gc_extreme_filter`) *before* HMM calling.
  Both default `TRUE`. This runs independently of, and before, the
  low-variance reference selection described above.
- **QC metrics table** (`canope_qc_metrics.R`) — per-sample and per-exon
  flags (low correlation, low depth, high CV, missing sex chromosomes),
  written to TSV and shown in the report.
- **Confidence scoring** (`canope_confidence.R`) — HIGH/MEDIUM/LOW based on
  `Q_SOME`, `NUM_REFS`, `NUM_TARG`, and MLCN consistency, with a built-in
  list of pseudogene/homology-prone genes always scored LOW.
- **PCA of coverage profiles** (`canope_utils.R::plot_coverage_pca`).
- **VCF export** (`canope_vcf_export.R`) — multi-sample or per-sample.
- **Interactive HTML report** (`canope_report.R::generate_canope_report()`,
  assets in `inst/report/`) — see [The Interactive Report](#the-interactive-report).

## Fast coverage extraction (megadepth, Windows-native)

`data_utils.R` also has `get_coverage_from_bams_megadepth()` — a much faster
alternative to the default `summarizeOverlaps()`-based extraction, for when
`run_canope()` is pulling coverage straight from BAMs (i.e. no `reads_file`).

Switch backends with one argument in your `config.yaml` or directly in `run_canope()`:

```r
run_canope(
  ...,
  coverage_backend = "megadepth",   
  megadepth_op = "sum",             # or "mean"
  megadepth_threads = 4
)
```

## Optional legacy HMM engine (original CANOES)

```r
run_canope(..., engine = "legacy_canoes")   # or per-sample: call_cnvs(..., engine = "legacy_canoes")
```

If you request `decode_method = "stationary"` together with
`engine = "legacy_canoes"`, you'll get a warning and it falls back to
distance-based transitions — the original never had a stationary mode.

**Tests:** `devtools::test()` runs unit tests (chromosome status, sex inference, GC alignment, exon numbering, calibration) and a small end-to-end run on synthetic data (spiked DEL / DUP found, calls identical for different row orders, report built, "autosomes and gonosomes not both analysed" warning present).

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).