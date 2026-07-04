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

* **Smart Normalization:** Automatically selects the most highly correlated reference samples (up to a user-defined threshold) to build an expected baseline via Non-Negative Least Squares (NNLS).
* **Robust HMM Engine:** Implements a Hidden Markov Model with dynamically computed transition probabilities (based on genomic distance) and Negative Binomial emission probabilities.
* **Native Visualization:** Generates comprehensive PDF reports with tracking plots, coverage ratios, and specific gene annotations completely independent of heavy third-party plotting frameworks.
* **GC Content Correction:** Built-in tools (`compute_gc_content`) to directly parse GC fraction from FASTA files and correct GC-bias natively.

## Installation

You can install the development version of CANOPE directly from GitHub using `devtools`.

```R
# Install devtools if you haven't already
if (!require("devtools")) install.packages("devtools")

# Install CANOPE
devtools::install_github("lavauxt/CANOPE")
```

### Dependencies
CANOPE is lightweight by design and relies on the following R packages:
`dplyr`, `nnls`, `Hmisc`, `mgcv`, `stats`, `utils`, `ggplot2` (for plotting).

## Quick Start

### 1. Data Preparation
Ensure you have your raw coverage counts and GC content ready. CANOPE expects a `.tsv` format where the first columns represent coordinates (`chromosome`, `start`, `end`), followed by sample counts.

### 2. Running the Pipeline
You can run the entire detection pipeline using the `run_canope` wrapper function:

```R
library(CANOPE)

# Execute the pipeline directly on Autosomes
run_canope(
  gc_file     = "data/gc.tsv",
  reads_file  = "data/canope.reads.tsv",
  samples     = c("Sample1", "Sample2", "Sample3"), # Vector of test samples
  modechrom   = "A",                                 # "A" for Autosomes, "XX" or "XY" for sex chromosomes
  output_file = "results/CNVCall_Autosomes.csv",
  rdata_output= "results/canope_workspace.RData"
)
```

### 3. Visualizing Results
If you provided an `rdata_output` path during the run, you can generate detailed PDF visualizations of the called CNVs:

```R
generate_plots(
  rdata_file = "results/canope_workspace.RData",
  output_dir = "results/plots",
  modechrom  = "A"
)
```

## Input File Formats

**GC File (`gc.tsv`)**
Must contain a column exactly named `GC_CONTENT` as the 4th column. *Tip: Use CANOPE's `compute_gc_content(fasta_file)` function to generate this easily.*

**Reads File (`canope.reads.tsv`)**
Tab-separated file containing the raw read counts for all samples across targeted regions.

| chromosome | start | end | Sample1 | Sample2 | Sample3 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| chr1 | 10000 | 10150 | 45 | 50 | 48 |
| chr1 | 10500 | 10620 | 12 | 15 | 11 |

## Output Format (`CNVCall.csv`)

The output is a tabular file containing the detected variations:

* `Chrom`, `Start`, `End`: Genomic coordinates of the CNV.
* `CNV`: Type of variant (`DEL` for deletion, `DUP` for duplication).
* `SAMPLE`: The sample in which the event was detected.
* `KB`: Size of the variation in kilobases.
* `TARGETS`: The indices of the targeted exons spanning the event.
* `NUM_TARG`: Total number of contiguous targets involved.
* `MLCN`: Most Likely Copy Number (e.g., 1 for heterozygous deletion, 3 for duplication).
* `Q_SOME`: Phred-scaled quality score of the CNV call.


## New features

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


## Fast coverage extraction (megadepth, Windows-native)

`data_utils.R` also has `get_coverage_from_bams_megadepth()` — a much faster
alternative to the default `summarizeOverlaps()`-based extraction, for when
`run_canope()` is pulling coverage straight from BAMs (i.e. no `reads_file`).

```r
BiocManager::install("megadepth")  
library(megadepth); install_megadepth()
source("data_utils.R")
test_counts <- get_coverage_from_bams_megadepth("one_sample.bam", "panel.bed")
head(test_counts)
```

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
  coverage_backend = "megadepth",   
  megadepth_op = "sum",             # or "mean"
  megadepth_threads = 4
)
```
## Optional legacy HMM engine (original CANOES)

**What's actually different between the two**, since the underlying
statistical model turned out to already be the same one (the "new" engine in
`hmm_engine.R` is itself a refactor of CANOES, not a different model):

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

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).