#' Run the CANOPE Copy Number Calling Pipeline
#'
#' End-to-end wrapper: reads count/GC data, QC-filters samples and exons,
#' optionally infers reference sets, runs per-sample HMM calling, and writes
#' output CSV + RData + plots. Also orchestrates BED preprocessing, CNV
#' confidence scoring, QC metrics, PCA, VCF export, and an interactive HTML
#' report.
#'
#' @param gc_file            Character or NULL. Path to pre-computed GC table.
#' @param fasta_file         Character or NULL. Path to genome FASTA (used with bed_file).
#' @param bsgenome_pkg      Character or NULL. Name of a BSgenome package (e.g. "BSgenome.Hsapiens.UCSC.hg38") – used if \code{fasta_file} is not provided.
#' @param reads_file         Character or NULL. Path to pre-computed counts table
#'   (expected layout: \code{chromosome, start, end, GENE, <sample columns...>}).
#' @param modechrom          Character. "A" (autosomes), "XX", or "XY".
#' @param removeY            Logical. If TRUE, exclude chrY targets.
#' @param samples            Character. Directory of BAMs or path to BAM list file.
#' @param p_value            Numeric. Prior CNV probability per target.
#' @param Tnum               Numeric. Expected number of targets per CNV.
#' @param D                  Numeric. Expected genomic span of a CNV (bp). Default 100 kb.
#' @param numrefs            Integer. Max reference samples.
#' @param min_cor            Numeric or NULL. Minimum Spearman correlation for references.
#' @param homdel_mean        Numeric. Expected read count under homozygous deletion.
#' @param auto_reference     Logical. Auto-infer reference sets per sample.
#' @param sample_qc          Logical. Exclude outlier samples before calling.
#' @param exon_qc            Logical. Exclude problematic exons before calling.
#' @param qc_zscore          Numeric. Robust z-score threshold for sample QC.
#' @param exon_mad_quantile  Numeric. MAD quantile threshold for exon QC.
#' @param gc_extreme_filter  Numeric[2]. Min/max GC (fraction, 0-1) to retain a target.
#' @param min_exon_mean      Numeric. Minimum mean count to retain a target.
#' @param output_file        Character. Path for CSV output.
#' @param rdata_output       Character. Path for RData workspace output.
#' @param output_prefix      Character. Base prefix for all output files (e.g. "CANOPE_TEST").
#'                           If NULL, derived from output_file (strips "_CNVCall.csv").
#' @param refbams_file       Character or NULL. Optional reference BAM list.
#' @param ref_reads          Character or NULL. Pre-computed reference counts table.
#' @param bed_file           Character or NULL. BED file for BAM coverage extraction
#'   (also used for GC computation when \code{fasta_file} is supplied).
#' @param coverage_backend   \code{"bioconductor"} (default) or \code{"megadepth"}.
#'   Only matters when extracting coverage directly from \code{samples}' BAMs
#'   (i.e. \code{reads_file} is NULL). \code{"bioconductor"} uses
#'   \code{\link{get_coverage_from_bams}} (\code{GenomicAlignments::summarizeOverlaps}).
#'   \code{"megadepth"} uses \code{\link{get_coverage_from_bams_megadepth}}, which
#'   shells out to the compiled megadepth tool — much faster, and (unlike
#'   mosdepth/PanDepth) has an official native Windows binary, auto-downloaded by
#'   the Bioconductor \code{megadepth} package. Note megadepth reports base-level
#'   coverage rather than a fragment count; see that function's docs for what this
#'   does and doesn't change downstream.
#' @param megadepth_op       \code{"sum"} (default) or \code{"mean"}; forwarded to
#'   \code{\link{get_coverage_from_bams_megadepth}} when
#'   \code{coverage_backend = "megadepth"}.
#' @param megadepth_threads  Integer. Threads per megadepth invocation; forwarded
#'   when \code{coverage_backend = "megadepth"}.
#' @param decode_method      \code{"distance"} (default) or \code{"stationary"} HMM transitions.
#'   Ignored under \code{engine = "legacy_canoes"} (forced to \code{"distance"} with a warning).
#' @param engine              \code{"new"} (default) or \code{"legacy_canoes"}. Selects the HMM
#'   core used by \code{\link{call_cnvs}} for every sample in this run.
#'   \code{"legacy_canoes"} is a faithful, numerically-verified port of the
#'   original published CANOES HMM (see \code{canoes_legacy_engine.R} for what
#'   that does and doesn't change vs. the modern \code{"new"} engine in
#'   \code{hmm_engine.R}). Reference selection, NNLS weighting, and variance
#'   estimation are identical either way — only the HMM transition/emission/
#'   decoding/Phred math differs.
#' @param bed_process        \code{"NO"} (default), \code{"STANDARD"}, or \code{"REGEN"}.
#'   If not \code{"NO"}, \code{bed_file} is first run through
#'   \code{\link{process_bed_file}} (see \code{bed_process_args}) before coverage
#'   extraction / GC computation.
#' @param bed_process_args   Named list of extra arguments forwarded to
#'   \code{\link{process_bed_file}} (e.g. \code{list(genome_version = "hg38")}).
#'   Also where off-target/filler-region handling is configured -- add
#'   \code{off_target_pattern}/\code{off_target_handling} keys here; see
#'   \code{\link{process_bed_file}} and \code{\link{handle_off_target_regions}}.
#'   Ported from ECHO.
#' @param pad_terminal_exons Integer >= 0. Bases to extend the outward-facing
#'   edge of each gene's first and last exon (both edges, for a single-exon
#'   gene), applied to \code{bed_file} before coverage/GC extraction --
#'   reduces the chance of a 0/low count right at a gene boundary where
#'   there's no neighbouring exon to carry the signal. Never creates an
#'   overlap with a neighbouring interval or crosses a contig boundary;
#'   padding is applied "if possible" and clamped short otherwise. Only
#'   takes effect when extracting fresh from BAMs (i.e. \code{reads_file}
#'   is \code{NULL}) -- a pre-computed \code{reads_file} was already
#'   extracted over some fixed window, so padding \code{bed_file} alone
#'   afterward would desync GC content from the counts it's paired with.
#'   Ported from ECHO; see \code{\link{pad_gene_terminal_exons}}. Default
#'   \code{0} (disabled).
#' @param plot_gene_gap      Numeric >= 0. Extra x-axis units inserted
#'   between a gene's last exon and the next gene's first exon in every PDF
#'   and HTML-report panel, so a window that spans more than one gene
#'   doesn't read as one continuous, unbroken feature. Forwarded to
#'   \code{\link{generate_plots}} and to the interactive report. Ported
#'   from ECHO; see \code{\link{compute_gene_gap_positions}}. Default
#'   \code{1}.
#' @param run_qc_metrics     Logical. Compute and write CANOPE QC metrics
#'   (see \code{\link{run_canope_qc_metrics}}). Default TRUE.
#' @param qc_output_file     Path for the QC metrics TSV. Default derived from
#'   \code{output_file}'s directory.
#' @param qc_min_corr,qc_min_cov,qc_min_total_reads,qc_max_exon_cv
#'   Thresholds forwarded to \code{\link{run_canope_qc_metrics}}.
#' @param score_confidence   Logical. Annotate \code{Confidence}/\code{CN_label}
#'   via \code{\link{score_canope_confidence}}. Default TRUE.
#' @param confidence_args    Named list of extra arguments forwarded to
#'   \code{\link{score_canope_confidence}}.
#' @param generate_plots_flag Logical. Generate per-call static PDF plots. Default TRUE.
#' @param pca_plot           Logical. Generate a PCA-of-coverage PDF. Default TRUE.
#' @param pca_output_file    Path for the PCA PDF. Default derived from \code{output_file}.
#' @param sample_table       Optional path to a TSV with \code{sample_name} (and
#'   optionally \code{gender}) columns, used to colour the PCA plot.
#' @param export_vcf         Logical. Export CNV calls to VCF. Default TRUE.
#' @param vcf_output         Path for the combined multi-sample VCF. Default derived
#'   from \code{output_file}.
#' @param vcf_per_sample     Logical. Also write one VCF per sample. Default FALSE.
#' @param report             Logical. Render the interactive HTML report
#'   (\code{\link{generate_canope_report}}). Default TRUE.
#' @param report_output_dir  Directory for the rendered HTML. Default derived from
#'   \code{output_file}.
#' @param report_template    Optional explicit path to \code{CANOPE_report.Rmd}.
#' @param log_file           Optional path to write a pipeline log (also shown in
#'   the HTML report). Default derived from \code{output_file}'s directory.
#'
#' @return Data frame of all called CNVs across samples (confidence-scored if
#'   \code{score_confidence = TRUE}).
#' @export
run_canope <- function(
    gc_file = NULL, fasta_file = NULL, bsgenome_pkg = NULL, reads_file = NULL,
    modechrom = "A", removeY = TRUE, samples = NULL,
    p_value = 1e-08, Tnum = 6, D = 100000,
    numrefs = 30, min_cor = NULL, homdel_mean = 0.2,
    auto_reference = TRUE, sample_qc = TRUE, exon_qc = TRUE,
    qc_zscore = 3, exon_mad_quantile = 0.90,
    gc_extreme_filter = c(0.15, 0.85), min_exon_mean = 20,
    output_file = "CNVCall.csv", rdata_output = "canope_workspace.RData",
    output_prefix = NULL,
    refbams_file = NULL, ref_reads = NULL, bed_file = NULL,
    coverage_backend = c("bioconductor", "megadepth"),
    megadepth_op = c("sum", "mean"), megadepth_threads = 1L,
    decode_method = c("distance", "stationary"),
    engine = c("new", "legacy_canoes"),
    bed_process = "NO", bed_process_args = list(),
    pad_terminal_exons = 0, plot_gene_gap = 1,
    run_qc_metrics = TRUE, qc_output_file = NULL,
    qc_min_corr = 0.98, qc_min_cov = 100, qc_min_total_reads = 300000, qc_max_exon_cv = 0.5,
    score_confidence = TRUE, confidence_args = list(),
    generate_plots_flag = TRUE,
    pca_plot = TRUE, pca_output_file = NULL, sample_table = NULL,
    export_vcf = TRUE, vcf_output = NULL, vcf_per_sample = FALSE,
    report = TRUE, report_output_dir = NULL, report_template = NULL,
    log_file = NULL
) {
  decode_method <- match.arg(decode_method)
  engine <- match.arg(engine)
  coverage_backend <- match.arg(coverage_backend)
  megadepth_op <- match.arg(megadepth_op)

  # ---- Derive output prefix if not given ----
  if (is.null(output_prefix)) {
    base <- basename(output_file)
    if (grepl("_CNVCall\\.csv$", base)) {
      output_prefix <- sub("_CNVCall\\.csv$", "", base)
    } else {
      output_prefix <- tools::file_path_sans_ext(base)
    }
  }
  # Ensure prefix starts with "CANOPE_"
  if (!grepl("^CANOPE_", output_prefix)) {
    output_prefix <- paste0("CANOPE_", output_prefix)
  }

  out_dir <- dirname(output_file)
  if (out_dir == "") out_dir <- "."

  # ---- Use prefix for all derived output paths ----
  if (is.null(log_file)) {
    log_file <- file.path(out_dir, paste0(output_prefix, "_pipeline.log"))
  }
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
  log_con <- file(log_file, open = "a")
  log_msg <- function(level, ...) {
    txt <- paste0("[", level, "] ", paste0(..., collapse = ""))
    message(txt)
    writeLines(txt, log_con)
  }
  on.exit(close(log_con), add = TRUE)

  log_msg("INFO", "Preparing count matrix")
  log_msg("INFO", sprintf("Terminal-exon padding: %s bp | Plot gene-gap: %s",
                          pad_terminal_exons %||% 0, plot_gene_gap %||% 1))

  # ── BED preprocessing (optional) ─────────────────────────────────────────
  if (!is.null(bed_file) && !identical(bed_process, "NO")) {
    processed_bed <- file.path(out_dir, paste0(output_prefix, "_processed.bed"))
    log_msg("INFO", sprintf("Running process_bed_file(mode = '%s') on %s", bed_process, bed_file))
    do.call(process_bed_file, c(
      list(input_bed = bed_file, output_bed = processed_bed, bed_process = bed_process),
      bed_process_args
    ))
    bed_file <- processed_bed
  }

  if (length(samples) == 1 && dir.exists(samples)) {
    samplesbams <- list.files(samples, pattern = "\\.bam$", full.names = TRUE)
    if (length(samplesbams) == 0) stop("ERROR: No .bam files found in directory.")
  } else if (length(samples) == 1 && file.exists(samples)) {
    samplesbams <- readLines(samples)
    samplesbams <- samplesbams[nzchar(trimws(samplesbams))]
  } else {
    stop("ERROR: 'samples' must be a directory or a valid file path.")
  }

  samples_to_analyse <- sapply(samplesbams, clean_name, USE.NAMES = FALSE)

  if (!is.null(reads_file) && file.exists(reads_file)) {
    canope.reads_un <- utils::read.table(reads_file, header = TRUE,
                                         check.names = FALSE,
                                         stringsAsFactors = FALSE)
    n_expected_meta <- 4L
    if (ncol(canope.reads_un) - n_expected_meta != length(samples_to_analyse)) {
      stop(sprintf(
        "[ERROR] reads_file has %d columns after the expected %d metadata columns ",
        ncol(canope.reads_un) - n_expected_meta, n_expected_meta),
        sprintf("(chromosome, start, end, GENE), but %d samples were supplied. ",
                length(samples_to_analyse)),
        "Check that reads_file has exactly one column per sample after those 4 metadata columns."
      )
    }
    colnames(canope.reads_un)[seq(n_expected_meta + 1L, ncol(canope.reads_un))] <- samples_to_analyse

  } else {
    # ── Terminal-exon padding (optional; ported from ECHO) ─────────────────
    if (!is.null(pad_terminal_exons) && !is.na(pad_terminal_exons) && pad_terminal_exons > 0) {
      chr_lengths_for_padding <- NULL
      if (!is.null(fasta_file) && file.exists(fasta_file)) {
        chr_lengths_for_padding <- tryCatch(
          GenomeInfoDb::seqlengths(Rsamtools::FaFile(fasta_file)),
          error = function(e) {
            log_msg("WARNING", "Could not read contig lengths from fasta_file for padding clamp: ", conditionMessage(e))
            NULL
          })
      } else if (!is.null(bsgenome_pkg) && requireNamespace(bsgenome_pkg, quietly = TRUE)) {
        chr_lengths_for_padding <- tryCatch(
          GenomeInfoDb::seqlengths(getExportedValue(bsgenome_pkg, bsgenome_pkg)),
          error = function(e) {
            log_msg("WARNING", "Could not read contig lengths from bsgenome_pkg for padding clamp: ", conditionMessage(e))
            NULL
          })
      }
      padded_bed <- file.path(out_dir, paste0(output_prefix, "_padded.bed"))
      log_msg("INFO", sprintf("Applying terminal-exon padding (%d bp) to %s", pad_terminal_exons, bed_file))
      pad_bed_file(bed_file, padded_bed, padding = pad_terminal_exons, chr_lengths = chr_lengths_for_padding)
      bed_file <- padded_bed
    }

    log_msg("INFO", sprintf("Auto-coverage extraction from BAMs (backend = '%s')...", coverage_backend))
    counts_df <- switch(
      coverage_backend,
      "bioconductor" = get_coverage_from_bams(bam_files = samplesbams, bed_input = bed_file),
      "megadepth"    = get_coverage_from_bams_megadepth(
        bam_files = samplesbams, bed_input = bed_file,
        op = megadepth_op, threads = megadepth_threads
      )
    )
    colnames(counts_df) <- samples_to_analyse

    bed_gr      <- rtracklayer::import(bed_file)
    coords_full <- as.data.frame(bed_gr)
    gene_names  <- if ("name" %in% colnames(coords_full))
      coords_full$name else paste0("Target_", seq_len(nrow(coords_full)))

    coords <- data.frame(
      chromosome = coords_full$seqnames,
      start      = coords_full$start,
      end        = coords_full$end,
      GENE       = gene_names,
      stringsAsFactors = FALSE
    )

    dup_mask  <- duplicated(coords[, c("chromosome", "start", "end")])
    coords    <- coords[!dup_mask, , drop = FALSE]
    counts_df <- counts_df[!dup_mask, , drop = FALSE]

    canope.reads_un <- cbind(coords, counts_df)
  }

  meta_cols <- c("chromosome", "start", "end", "GENE", "target", "gc")
  sample_cols <- setdiff(colnames(canope.reads_un), meta_cols)
  for (col in sample_cols) {
    raw_vals <- canope.reads_un[[col]]
    canope.reads_un[[col]] <- suppressWarnings(as.numeric(as.character(raw_vals)))
    lost <- is.na(canope.reads_un[[col]]) & !is.na(raw_vals) & nzchar(as.character(raw_vals))
    if (any(lost)) {
      warning(sprintf("Column '%s' had %d non-numeric value(s) coerced to NA", col, sum(lost)), call. = FALSE)
    }
  }
  # ====================================================

  # ---- GC content --------------------------------------------
  if (!is.null(fasta_file) && !is.null(bed_file)) {
    log_msg("INFO", "Computing GC content from FASTA file")
    datagc <- compute_gc_from_fasta(fasta_file = fasta_file, bed_input = bed_file)
  } else if (!is.null(bsgenome_pkg) && !is.null(bed_file)) {
    if (!requireNamespace(bsgenome_pkg, quietly = TRUE))
      stop("BSgenome package '", bsgenome_pkg, "' is not installed.")
    log_msg("INFO", sprintf("Computing GC content from BSgenome package '%s'", bsgenome_pkg))
    datagc <- compute_gc_from_bed(bsgenome_pkg = bsgenome_pkg, bed_input = bed_file)
  } else if (!is.null(gc_file) && file.exists(gc_file)) {
    log_msg("INFO", "Using pre‑computed GC file")
    datagc <- utils::read.table(gc_file, header = TRUE)
  } else {
    stop("ERROR: Provide either (fasta_file + bed_file) OR (bsgenome_pkg + bed_file) OR a pre‑computed gc_file.")
  }

  if (!"GC_CONTENT" %in% colnames(datagc)) {
    gc_col_guess <- grep("^gc", colnames(datagc), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(gc_col_guess)) {
      colnames(datagc)[colnames(datagc) == gc_col_guess] <- "GC_CONTENT"
    } else if (ncol(datagc) >= 4) {
      colnames(datagc)[4] <- "GC_CONTENT"
    } else {
      stop("[ERROR] Could not identify a GC content column in 'datagc'.")
    }
  }
  gc <- as.numeric(datagc$GC_CONTENT)

  if (max(gc, na.rm = TRUE) > 1) {
    log_msg("INFO", "GC content appears to be on a 0-100 scale; converting to fraction.")
    gc <- gc / 100
  }


  all_sample_names <- samples_to_analyse
  canope.reads_un <- cbind(
    target = seq_len(nrow(canope.reads_un)),
    gc     = gc,
    GENE   = canope.reads_un$GENE,
    canope.reads_un[, !(colnames(canope.reads_un) %in% c("GENE", "target", "gc"))]
  )

  # ---- optional reference BAM panel --------------------------
  refsample_names <- character(0)
  if (!is.null(refbams_file) && file.exists(refbams_file)) {
    if (is.null(ref_reads) || !file.exists(ref_reads))
      stop("[ERROR] 'refbams_file' was supplied but 'ref_reads' is missing or does not exist. ",
           "'ref_reads' must be a pre-computed counts table for the reference panel.")

    rawrefbams <- utils::read.csv(refbams_file, header = TRUE, sep = "\t")
    if ("gender" %in% colnames(rawrefbams)) {
      if (modechrom == "XX")      rawrefbams <- subset(rawrefbams, gender == "F")
      else if (modechrom == "XY") rawrefbams <- subset(rawrefbams, gender == "M")
    }
    refbams         <- apply(rawrefbams, 1, toString)
    refsample_names <- tools::file_path_sans_ext(basename(refbams))

    # Same 4-metadata-column convention as 'reads_file' (chromosome, start,
    # end, GENE) -- previously this dropped only the first 3 columns, which
    # silently kept GENE as if it were the first reference sample's read
    # counts whenever ref_reads used the same layout as reads_file.
    data_ref <- utils::read.table(ref_reads, header = TRUE, check.names = FALSE,
                                  stringsAsFactors = FALSE)
    n_ref_meta <- 4L
    if (ncol(data_ref) <= n_ref_meta)
      stop("[ERROR] ref_reads has no sample columns after the expected ",
           n_ref_meta, " metadata columns (chromosome, start, end, GENE).")
    canope.reads_ref <- data_ref[, seq(n_ref_meta + 1L, ncol(data_ref)), drop = FALSE]

    missing_refs <- setdiff(refsample_names, colnames(canope.reads_ref))
    if (length(missing_refs) > 0)
      stop(sprintf(
        "[ERROR] %d reference sample(s) from refbams_file not found as columns in ref_reads: %s",
        length(missing_refs), paste(missing_refs, collapse = ", ")))

    canope.reads_un <- cbind(canope.reads_un, canope.reads_ref[, refsample_names, drop = FALSE])
  }

  target_samples <- if (length(refsample_names) > 0)
    c(samples_to_analyse, refsample_names) else all_sample_names

  canope.reads <- canope.reads_un[, c("target", "gc", "GENE", "chromosome", "start", "end", target_samples)]

  # ===== SAFETY: ensure all sample columns are numeric (again) =====
  for (nm in target_samples) {
    raw_vals <- canope.reads[[nm]]
    canope.reads[[nm]] <- as.numeric(as.character(raw_vals))
    lost <- is.na(canope.reads[[nm]]) & !is.na(raw_vals) & nzchar(as.character(raw_vals))
    if (any(lost)) {
      warning(sprintf("Column '%s' had %d non-numeric value(s) coerced to NA", nm, sum(lost)), call. = FALSE)
    }
  }
  # ================================================================

  has_chr_prefix <- any(grepl("^chr", as.character(canope.reads$chromosome)))
  chrX_label <- if (has_chr_prefix) "chrX" else "X"
  chrY_label <- if (has_chr_prefix) "chrY" else "Y"

  if (modechrom == "XX")
    canope.reads <- subset(canope.reads, chromosome == chrX_label)
  else if (modechrom == "XY")
    canope.reads <- subset(canope.reads, chromosome %in% c(chrX_label, chrY_label))
  else if (modechrom == "A")
    canope.reads <- subset(canope.reads, !chromosome %in% c(chrX_label, chrY_label))

  if (removeY)
    canope.reads <- subset(canope.reads, chromosome != chrY_label)

  if (nrow(canope.reads) == 0)
    stop(sprintf(
      "[ERROR] No targets remain after modechrom = '%s' filtering. Detected chromosome naming style: %s.",
      modechrom, if (has_chr_prefix) "'chr'-prefixed" else "bare (no 'chr' prefix)"
    ))

  canope.reads$start <- as.numeric(canope.reads$start)
  canope.reads$end   <- as.numeric(canope.reads$end)

  chrom_base   <- c(as.character(1:22), "X", "Y", "M")
  chrom.levels <- if (has_chr_prefix) paste0("chr", chrom_base) else chrom_base
  canope.reads$chromosome <- factor(canope.reads$chromosome, levels = chrom.levels)
  canope.reads <- canope.reads[order(canope.reads$chromosome, canope.reads$start, canope.reads$end), ]
  canope.reads <- canope.reads[!duplicated(canope.reads[, c("chromosome", "start", "end")]), ]
  rownames(canope.reads) <- NULL

  count_matrix <- as.matrix(canope.reads[, target_samples, drop = FALSE])
  storage.mode(count_matrix) <- "numeric"

  keep_gc      <- canope.reads$gc >= gc_extreme_filter[1] & canope.reads$gc <= gc_extreme_filter[2]
  canope.reads <- canope.reads[keep_gc, , drop = FALSE]
  count_matrix <- count_matrix[keep_gc, , drop = FALSE]

  if (nrow(canope.reads) == 0)
    stop(sprintf(
      "[ERROR] No targets remain after gc_extreme_filter = [%.2f, %.2f]. Check GC values are fractions (0-1).",
      gc_extreme_filter[1], gc_extreme_filter[2]
    ))

  count_matrix <- gc_correct_counts(count_matrix, gc = canope.reads$gc, method = "loess")
  count_matrix <- round(pmax(count_matrix, 0))
  canope.reads[, target_samples] <- count_matrix

  if (sample_qc) {
    sqc  <- detect_outlier_samples(count_matrix, z_threshold = qc_zscore)
    excl <- sqc$sample[sqc$is_outlier]
    if (length(excl) > 0) {
      log_msg("WARNING", sprintf("Excluding %d outlier sample(s): %s",
                                 length(excl), paste(excl, collapse = ", ")))
      target_samples     <- setdiff(target_samples, excl)
      samples_to_analyse <- intersect(samples_to_analyse, target_samples)
      count_matrix       <- count_matrix[, target_samples, drop = FALSE]
      canope.reads       <- canope.reads[, c("target", "gc", "GENE", "chromosome", "start", "end", target_samples)]
    }
  }

  if (exon_qc) {
    eqc <- detect_problematic_exons(
      count_matrix,
      chromosomes  = canope.reads$chromosome,
      mad_quantile = exon_mad_quantile,
      min_mean     = min_exon_mean,
      gc           = canope.reads$gc
    )
    good_exons   <- eqc$exon[!eqc$problematic]
    canope.reads <- canope.reads[good_exons, , drop = FALSE]
    count_matrix <- count_matrix[good_exons, , drop = FALSE]
    rownames(canope.reads) <- NULL
  }

  required_first <- c("target", "chromosome", "start", "end", "gc")
  missing_meta <- setdiff(required_first, colnames(canope.reads))
  if (length(missing_meta) > 0)
    stop("Missing required metadata columns: ", paste(missing_meta, collapse = ", "))
  other_cols <- setdiff(colnames(canope.reads), required_first)
  if ("GENE" %in% other_cols) other_cols <- c("GENE", setdiff(other_cols, "GENE"))
  canope.reads <- canope.reads[, c(required_first, other_cols)]

  auto_refs <- NULL
  if (auto_reference) {
    ref_obj   <- infer_reference_samples(count_matrix, top_n = numrefs)
    auto_refs <- ref_obj$refs
  }

  samples_to_analyse <- intersect(samples_to_analyse, colnames(canope.reads))

  all_cnvs <- list()
  models_list <- list()
  refs_list   <- list()

  for (samp in samples_to_analyse) {
    log_msg("INFO", sprintf("Calling CNVs for sample %s ...", samp))
    tryCatch({
      refs_to_use <- if (auto_reference && !is.null(auto_refs)) auto_refs[[samp]] else refsample_names
      if (auto_reference && is.null(refs_to_use)) {
        log_msg("WARNING", sprintf(
          "Sample '%s' was excluded from the auto-inferred reference panel (likely flagged as noisy); falling back to all other samples as candidate references.",
          samp
        ))
      }

      res <- call_cnvs(
        sample_name = samp, counts = canope.reads,
        p = p_value, Tnum = Tnum, D = D, numrefs = numrefs,
        get_dfs = FALSE, homdel_mean = homdel_mean,
        refsample_names = refs_to_use, full_output = TRUE,
        min_cor = min_cor, decode_method = decode_method, engine = engine
      )

      n_called <- if (!is.null(res) && !is.null(res$cnvs)) nrow(res$cnvs) else 0L
      log_msg("INFO", sprintf("  -> %d CNV(s) found for %s", n_called, samp))

      if (!is.null(res) && nrow(res$cnvs) > 0) all_cnvs[[samp]] <- res$cnvs
      refs_list[[samp]] <- res$reference_samples
      models_list[[samp]] <- list(
        target = res$targets, mean = res$mean, var_estimate = res$var_estimate,
        test_counts = res$test_counts, ref_matrix = res$ref_matrix
      )
    }, error = function(e) {
      log_msg("WARNING", sprintf("Sample %s failed: %s", samp, conditionMessage(e)))
    })
  }

  final_cnvs <- if (length(all_cnvs) > 0) do.call(rbind, all_cnvs) else data.frame()
  rownames(final_cnvs) <- NULL

  # ── CNV confidence scoring ───────────────────────────────────────────────
  if (score_confidence && nrow(final_cnvs) > 0) {
    log_msg("INFO", "Scoring CNV confidence")
    final_cnvs <- do.call(score_canope_confidence, c(list(cnv_calls = final_cnvs), confidence_args))
  }

  for (path in c(output_file, rdata_output)) {
    d <- dirname(path)
    if (d != "." && !dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  plot_dir <- file.path(dirname(output_file), "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  if (!is.null(output_file) && nrow(final_cnvs) > 0)
    utils::write.csv(final_cnvs, output_file, row.names = FALSE)

  if (!is.null(rdata_output)) {
    cnv_calls <- final_cnvs
    counts    <- canope.reads
    bed_file_obj <- canope.reads[, c("chromosome", "start", "end", "target", "gc", "GENE")]
    bed_file <- bed_file_obj
    models <- models_list
    refs   <- refs_list
    hmm_engine_used <- engine
    save(cnv_calls, counts, bed_file, models, refs, hmm_engine_used, file = rdata_output)

    if (generate_plots_flag && nrow(final_cnvs) > 0) {
      log_msg("INFO", "Generating static CNV plots")
      generate_plots(
        rdata_file = rdata_output,
        output_dir = plot_dir,
        modechrom = modechrom,
        prefix = output_prefix,
        gene_gap = plot_gene_gap
      )
    }
  }

  # ── QC metrics ───────────────────────────────────────────────────────────
  qc_metrics_path <- NULL
  if (run_qc_metrics) {
    if (is.null(qc_output_file)) {
      qc_output_file <- file.path(out_dir, paste0(output_prefix, "_QC_metrics.tsv"))
    }
    log_msg("INFO", "Computing QC metrics")
    tryCatch({
      bed_for_qc <- canope.reads[, c("chromosome", "start", "end", "GENE")]
      run_canope_qc_metrics(
        counts = canope.reads, sample_names = target_samples, bed_df = bed_for_qc,
        output_file = qc_output_file,
        min_corr = qc_min_corr, min_cov = qc_min_cov,
        min_total_reads = qc_min_total_reads, max_exon_cv = qc_max_exon_cv
      )
      qc_metrics_path <- qc_output_file
    }, error = function(e) log_msg("WARNING", "QC metrics step failed: ", conditionMessage(e)))
  }

  # ── PCA of coverage profiles ─────────────────────────────────────────────
  if (pca_plot && length(target_samples) >= 3) {
    if (is.null(pca_output_file)) {
      pca_output_file <- file.path(out_dir, paste0(output_prefix, "_PCA.pdf"))
    }
    log_msg("INFO", "Generating coverage PCA plot")
    color_by <- NULL
    if (!is.null(sample_table) && file.exists(sample_table)) {
      st <- utils::read.table(sample_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      if (all(c("sample_name", "gender") %in% colnames(st)))
        color_by <- st$gender[match(target_samples, st$sample_name)]
    }
    tryCatch(
      plot_coverage_pca(canope.reads, target_samples, output_pdf = pca_output_file, color_by = color_by),
      error = function(e) log_msg("WARNING", "PCA plot step failed: ", conditionMessage(e))
    )
  }

  # ── VCF export ───────────────────────────────────────────────────────────
  if (export_vcf && nrow(final_cnvs) > 0) {
    if (is.null(vcf_output)) {
      vcf_output <- file.path(out_dir, paste0(output_prefix, "_calls.vcf"))
    }
    log_msg("INFO", "Exporting CNV calls to VCF")
    tryCatch(export_canope_to_vcf(final_cnvs, vcf_output),
            error = function(e) log_msg("WARNING", "VCF export failed: ", conditionMessage(e)))
    if (vcf_per_sample) {
      for (samp in unique(final_cnvs$SAMPLE)) {
        tryCatch(
          export_canope_to_vcf(final_cnvs, file.path(out_dir, paste0(sanitize_filename(samp), ".vcf")),
                               sample_name = samp),
          error = function(e) log_msg("WARNING", sprintf("VCF export failed for %s: %s", samp, conditionMessage(e)))
        )
      }
    }
  }

  # ── Interactive HTML report ──────────────────────────────────────────────
  if (report && nrow(final_cnvs) > 0) {
    if (is.null(report_output_dir)) report_output_dir <- out_dir
    # Strip leading "CANOPE_" to avoid duplication (report will add it)
    report_prefix <- sub("^CANOPE_", "", output_prefix)
    log_msg("INFO", "Rendering interactive HTML report")
    tryCatch(
      generate_canope_report(
        rdata_output = rdata_output,
        qc_metrics_file = qc_metrics_path,
        output_dir = report_output_dir,
        prefix = report_prefix,
        settings = c(list(low_confidence_genes = c("PMS2", "SMN1", "CYP2D6", "HBA1", "HBA2",
                                                    "STRC", "CYP21A2", "GBA1", "CFTR"),
                          qc_min_cov = qc_min_cov, qc_min_total_reads = qc_min_total_reads,
                          plot_gene_gap = plot_gene_gap),
                    confidence_args),
        sample_table = sample_table, log_file = log_file,
        template_path = report_template
      ),
      error = function(e) log_msg("WARNING", "Report generation failed: ", conditionMessage(e))
    )
  }

  log_msg("INFO", sprintf("Pipeline complete: %d CNV(s) called across %d sample(s)",
                          nrow(final_cnvs), length(samples_to_analyse)))

  return(final_cnvs)
}