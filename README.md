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
* **Native Visualization:** Generates comprehensive PDF reports with tracking plots, coverage ratios, and specific gene annotations completely independent of heavy third-party plotting frameworks.

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
`yaml`, `HMM`, `ggplot2`, and the Bioconductor packages `Biostrings`, `rtracklayer`,
`GenomicRanges`, `GenomeInfoDb`, `BiocGenerics`, `Rsamtools`,
`GenomicAlignments`, `SummarizedExperiment`, `IRanges`, `S4Vectors`.

Optional (only needed for specific features): `megadepth` (fast BAM
coverage), `GenomicFeatures` + a `TxDb.Hsapiens.UCSC.*` package + `org.Hs.eg.db`
(BED `REGEN` mode), `ggrepel` (nicer PCA labels), and `rmarkdown`, `knitr`,
`DT`, `plotly`, `tidyr`, `htmltools`, `kableExtra` (interactive HTML report).

## Quick Start

### 1. Data Preparation
Ensure you have your BAM files and a BED file of target regions. If you don't have pre-computed count matrices, CANOPE will automatically extract coverage from your BAMs and compute GC content using a `BSgenome` package or an indexed FASTA file.

### 2. Create a Configuration File
Create a `config.yaml` file to define your inputs, outputs, and HMM parameters.

```yaml
input:
  bed: "data/targets.bed"
  bamdir: "data/bams/"

settings:
  bsgenome_pkg: "BSgenome.Hsapiens.UCSC.hg19"  # Or use fasta_file instead
  modechrom: "A"
  p_value: 1e-08
  Tnum: 6
  numrefs: 30
  coverage_backend: "bioconductor"             # or "megadepth"
  bed_process: "STANDARD"                      # or "NO" / "REGEN"
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
If `report: TRUE` is set in your config, an interactive HTML report will be generated in the output directory. Additionally, static PDF plots for each CNV call are saved in the `results/plots/` folder.

## Input File Formats

**BED File**
Standard 0-based BED file. If `bed_process` is set to `"STANDARD"` or `"REGEN"`, the input BED will be automatically annotated with gene/exon information before coverage extraction.

**Reads File (`canope.reads.tsv`)** *(Optional)*
If you already have pre-computed read counts, you can provide a tab-separated file instead of letting CANOPE extract from BAMs. It must contain exactly four metadata columns (`chromosome`, `start`, `end`, `GENE`) followed by one column of raw read counts per sample -- `run_canope()` matches columns by position after the fourth column, so an extra or missing metadata column will shift every sample's counts by one.

| chromosome | start | end | GENE | Sample1 | Sample2 | Sample3 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| chr1 | 10000 | 10150 | BRCA1 | 45 | 50 | 48 |
| chr1 | 10500 | 10620 | BRCA1 | 12 | 15 | 11 |

## Output Format (`CNVCall.csv`)

The output is a tabular file containing the detected variations:

* `CHR`, `INTERVAL`, `MID_BP`: Genomic coordinates of the CNV.
* `CNV`: Type of variant (`DEL` for deletion, `DUP` for duplication).
* `SAMPLE`: The sample in which the event was detected.
* `KB`: Size of the variation in kilobases.
* `TARGETS`: The target ID range spanning the event (e.g. `42..47`).
* `NUM_TARG`: Total number of contiguous targets involved.
* `GENE`: Gene name(s) overlapping the call.
* `MLCN`: Most Likely Copy Number (e.g., 1 for heterozygous deletion, 3 for duplication).
* `Q_SOME`: Phred-scaled quality score of the CNV call.
* `NUM_REFS`, `REF_SAMPLES`: Number and names of reference samples used.
* `Confidence`, `CN_label`: Added when `score_confidence = TRUE` (default).

## Additional Features

- **BED preprocessing** (`canope_bed_process.R`) — `STANDARD`/`REGEN`/`NO`
  modes, gene/exon annotation from panel files or RefSeq.
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
- **Interactive HTML report** (`CANOPE_report.Rmd` /
  `canope_report.R::generate_canope_report()`) — QC table, PCA, per-sample
  CNV tables colour-coded by confidence, and three plotly panels per call:
  coverage, log2(observed/expected) with NB confidence interval, and a
  z-score-vs-references panel (interactive version of the static PDF panel).

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

## Related project: ECHO

CANOPE has a sibling pipeline, **ECHO**, which targets the same use case
(exome/panel CNV calling) but calls variants with `ExomeDepth`'s
beta-binomial model instead of CANOPE's own HMM. The two share the same
overall pipeline shape — BED preprocessing, BAM coverage/GC extraction,
sample/exon QC, confidence scoring, VCF export, static PDF plots, and an
interactive HTML report — and most of that shared scaffolding (BED
preprocessing modes, exon-numbering, PCA, background-calibration
diagnostics, VCF/report layout) is kept in sync between the two so that
switching between them, or running both on the same panel, feels
consistent. Only the actual calling engine and its specific tuning
parameters (`p_value`/`Tnum`/`D`/`numrefs` here vs. ExomeDepth's
`transition.probability`/`expected.CNV.length`/`phi.bins` there) are
pipeline-specific by design.

See `CANOPE_vs_ECHO_COMPARISON.md` for the full feature-parity matrix
between the two pipelines, along with a changelog of bugs found and fixed
and parity features added during the most recent cross-pipeline audit.

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).