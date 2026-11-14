# ============================================================================
# CANOPE report: the adapter between a CANOPE workspace and the shared report.
#
# Everything the two pipelines have in common (QC tables, alerts, PCA, sample
# table, chromosome scope, JSON payload, HTML) lives in cnv_shared_report.R and
# cnv_report_assets.R, byte-identical in ECHO. This file only knows what is
# specific to CANOPE: how its workspace is laid out, how a call row maps to the
# unified call table, how a call window is drawn (canope_plot_data()), and the
# per-sample flags that only a CANOPE run can raise (few references, genotyping
# skipped, low-coverage runs, sex mismatch ...).
# ============================================================================

#' Generate the interactive HTML report of a CANOPE run
#'
#' Writes one self-contained HTML file (no server, no Shiny) from the workspace
#' written by \code{\link{run_canope}}. The page is the same as ECHO's: same tabs,
#' colours, filters and plots, only the data differ.
#'
#' @param rdata_output Path to the CANOPE workspace \code{.RData}.
#' @param qc_metrics_file Optional path to the QC-metrics TSV.
#' @param output_dir Directory for the HTML (created if missing).
#' @param settings Optional named list of thresholds / scoring rules (defaults come
#'   from the workspace).
#' @param sample_table Optional path to a sample sheet (\code{sample_name}, \code{gender}).
#' @param log_file Optional pipeline log (shown on the Run tab).
#' @param prefix Output name prefix (\code{CANOPE_} is added to the file name).
#' @param template_path Deprecated and ignored.
#' @param max_call_plots Maximum number of calls whose interactive plot is embedded.
#' @param max_qc_rows Row cap for the per-target QC rows (low depth, high variability,
#'   removed targets); sample-level rows are always kept.
#' @param flank Number of flanking targets drawn on each side of a call.
#' @param title Report title.
#' @param plotly_source \code{"inline"} (default, works offline) or \code{"cdn"}.
#' @param global \code{TRUE}: one report for all samples; \code{FALSE}: one per sample.
#' @param sample_name Samples to report on when \code{global = FALSE} (default: every
#'   sample with a call).
#' @return Invisibly, the path(s) of the HTML file(s) written.
#' @export
generate_canope_report <- function(rdata_output,
                                   qc_metrics_file = NULL,
                                   output_dir,
                                   settings = NULL,
                                   sample_table = NULL,
                                   log_file = NULL,
                                   prefix = NULL,
                                   template_path = NULL,
                                   max_call_plots = 300L,
                                   max_qc_rows = 5000L,
                                   flank = 5L,
                                   title = NULL,
                                   plotly_source = c("inline", "cdn"),
                                   global = TRUE,
                                   sample_name = NULL) {
  plotly_source <- match.arg(plotly_source)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required.")
  if (!is.null(template_path)) message("[INFO] 'template_path' is ignored.")
  stop_if_not_file(rdata_output, paste0("[ERROR] Workspace not found: ", rdata_output))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("[ERROR] Failed to create output directory: ", output_dir)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  settings <- settings %||% list()
  settings$report_plotly <- plotly_source
  settings$report_max_profiles <- as.integer(max_call_plots)
  settings$report_max_qc_rows <- as.integer(max_qc_rows)

  rd <- canope_report_inputs(rdata_output, qc_metrics_file, sample_table)
  base <- sub("^CANOPE_?", "", as.character(prefix %||% rd$run$prefix %||% ""))
  ad <- canope_report_adapter(rd, settings, flank = flank, title = title)

  report_emit(
    build = function(focus) report_build_payload(rd, settings, NULL, log_file, base, focus = focus, ad = ad),
    output_dir = output_dir,
    global_name = if (nzchar(base)) paste0("CANOPE_", base, "_report.html") else "CANOPE_global_report.html",
    sample_name = function(safe) paste0("CANOPE_report_", safe, ".html"),
    samples = if (global) NULL else (sample_name %||% unique(as.character(rd$cnv_calls$SAMPLE))),
    settings = settings, global = global)
}

# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------

# columns of the counts table that are not samples
CANOPE_META_COLS <- c("target", "chromosome", "start", "end", "gc", "GENE", "mean", "exon_number")

#' Read a CANOPE workspace into the shape the shared report builder expects
#' @noRd
canope_report_inputs <- function(rdata_output, qc_metrics_file = NULL, sample_table = NULL) {
  env <- load_rdata(rdata_output, required = c("cnv_calls", "counts", "bed_file", "models", "refs"))
  run <- env$run_info %||% list()
  counts <- env$counts
  bed_file <- canope_target_annotation(env$bed_file, counts)
  sample_cols <- setdiff(colnames(counts), CANOPE_META_COLS)

  calls <- canope_report_calls(env$cnv_calls)
  analysed <- as.character(run$samples_analysed %||% env$sample_names %||% sample_cols)
  requested <- as.character(run$samples_requested %||% analysed)

  # the cohort the sample metrics are computed on: raw counts of every sample, incl.
  # the ones excluded before calling (older workspaces only have the called ones)
  cohort <- env$raw_counts
  if (is.null(cohort) || is.null(colnames(cohort))) cohort <- as.matrix(counts[, sample_cols, drop = FALSE])
  all_samples <- unique(c(requested, analysed, as.character(calls$SAMPLE)))
  all_samples <- all_samples[all_samples %in% colnames(cohort)]
  analysed <- analysed[analysed %in% all_samples]

  sdf <- report_read_table(sample_table, "sample_name")
  if (is.null(sdf) && is.data.frame(run$sample_table) && "sample_name" %in% colnames(run$sample_table)) sdf <- run$sample_table

  list(cnv_calls = calls, counts = counts, bed_file = bed_file, raw_bed = env$raw_bed,
       models = env$models, refs = env$refs, analysed = analysed, all_samples = all_samples,
       all_chromosomes = as.character(env$input_chromosomes %||% env$raw_bed$chromosome %||% bed_file$chromosome),
       qc = report_read_qc(qc_metrics_file), sample_df = sdf, ref_bams = NULL,
       qc_info = env$qc_info, bed_dropped = NULL, cohort = cohort,
       modechrom = as.character(run$modechrom %||% env$modechrom %||% "A")[1],
       removeY = isTRUE(run$removeY %||% env$removeY %||% TRUE),
       engine = as.character(env$hmm_engine_used %||% run$engine %||% "new")[1],
       run = run)
}

#' A call table every later step can rely on (older workspaces lack some columns)
#' @noRd
canope_report_calls <- function(calls) {
  if (is.null(calls) || !ncol(calls) || !nrow(calls)) return(data.frame())
  if (!"CNV.ID" %in% colnames(calls)) calls$CNV.ID <- seq_len(nrow(calls))
  if (!all(c("Start", "End") %in% colnames(calls))) {
    iv <- lapply(as.character(calls$INTERVAL), parse_canope_interval)
    calls$Start <- vapply(iv, function(p) as.integer(p$start), integer(1))
    calls$End <- vapply(iv, function(p) as.integer(p$end), integer(1))
  }
  for (cn in c("Confidence", "CN_label")) if (!cn %in% colnames(calls)) calls[[cn]] <- NA_character_
  calls
}

# ----------------------------------------------------------------------------
# Per-sample flags only a CANOPE run can raise
# ----------------------------------------------------------------------------

#' QC rows (type "Whole sample") for the CANOPE-specific sample flags
#'
#' The text before the first colon is the flag's name: the alert list groups the
#' samples sharing it ("Few references in 3 sample(s)").
#' @noRd
canope_sample_flags <- function(rd, min_refs = 5L) {
  run <- rd$run; calls <- rd$cnv_calls
  st <- run$sample_status
  rows <- list()
  add <- function(sample, details, severity = "warning")
    rows[[length(rows) + 1L]] <<- data.frame(Sample = sample, Exon = "All", Type = "Whole sample",
                                             Details = details, Gene = "All", Severity = severity, stringsAsFactors = FALSE)
  analysed <- rd$analysed

  if (is.data.frame(st) && nrow(st)) {
    called <- st[st$status %in% "called", , drop = FALSE]
    few <- called[!is.na(called$n_refs) & called$n_refs < min_refs, , drop = FALSE]
    for (i in seq_len(nrow(few))) add(few$sample[i], sprintf("Few references: %d (< %d)", few$n_refs[i], min_refs))
    ng <- called[called$genotyped %in% FALSE, , drop = FALSE]
    for (i in seq_len(nrow(ng))) add(ng$sample[i], "Genotyping skipped: too many calls, Q_SOME is NA and every call is LOW", "fail")
  }

  # a sample with far more calls than the cohort is almost always a noisy one
  if (nrow(calls)) {
    n_calls <- vapply(analysed, function(s) sum(calls$SAMPLE == s), numeric(1))
    if (length(n_calls) >= 5L) {
      med <- stats::median(n_calls); cut <- max(10, med + 3 * max(stats::mad(n_calls), 1))
      for (s in names(n_calls)[n_calls > cut]) add(s, sprintf("Many calls: %d (cohort median %g)", as.integer(n_calls[[s]]), med))
    }
  }

  lcr <- run$low_coverage_runs
  if (is.data.frame(lcr) && nrow(lcr)) {
    for (s in unique(as.character(lcr$SAMPLE)))
      add(s, sprintf("Low-coverage runs: %d run(s) of >= 3 targets with ~no reads removed before calling (a homozygous deletion cannot be called there)",
                     sum(lcr$SAMPLE == s)))
  }

  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# ----------------------------------------------------------------------------
# Adapter
# ----------------------------------------------------------------------------

#' Everything CANOPE-specific the shared payload builder needs
#' @noRd
canope_report_adapter <- function(rd, settings, flank = 5L, title = NULL) {
  run <- rd$run; qi <- rd$qc_info; calls <- rd$cnv_calls
  counts <- rd$counts; bed_file <- rd$bed_file
  modechrom <- rd$modechrom; if (!modechrom %in% c("A", "XX", "XY")) modechrom <- "A"
  engine <- rd$engine
  scorer <- formals(score_canope_confidence)
  pick <- function(nm) settings[[nm]] %||% eval(scorer[[nm]])

  mode_label <- switch(modechrom, XX = "chrX, female references (XX)", XY = "chrX/chrY, male references (XY)", "Autosomes (A)")
  genome <- as.character(run$genome_build %||% "unknown")
  set_num <- function(nm, d) as.numeric(settings[[nm]] %||% run$thresholds[[nm]] %||% d)
  thresholds <- list(
    qc_min_cov = set_num("qc_min_cov", 100), min_cov = set_num("qc_min_cov", 100),
    qc_min_total_reads = set_num("qc_min_total_reads", 3e5), min_total_reads = set_num("qc_min_total_reads", 3e5),
    qc_min_corr = set_num("qc_min_corr", 0.98), min_corr = set_num("qc_min_corr", 0.98),
    qc_max_exon_cv = set_num("qc_max_exon_cv", 0.5), max_exon_cv = set_num("qc_max_exon_cv", 0.5),
    qc_zscore = as.numeric(qi$params$qc_zscore %||% run$params$qc_zscore %||% 3),
    outlier_sides = as.character(qi$params$outlier_sides %||% run$params$outlier_sides %||% "upper"))

  low_genes <- as.character(unlist(pick("low_confidence_genes")))
  conf_rules <- list(high_q = pick("high_q_score"), med_q = pick("med_q_score"),
                     high_refs = pick("high_num_refs"), med_refs = pick("med_num_refs"),
                     high_targ = pick("high_num_targ"), med_targ = pick("med_num_targ"),
                     genes = rep_arr(low_genes))
  rm_term <- isTRUE(settings$remove_terminal_only %||% run$params$remove_terminal_only)
  pen_term <- isTRUE(settings$penalize_terminal_only %||% run$params$penalize_terminal_only)
  scoring <- if (rm_term) list(terminal = "Terminal-only calls removed.")
             else if (pen_term) list(terminal = "Terminal-only calls downgraded to LOW.") else NULL

  # per-sample vectors the shared sample table understands
  st <- run$sample_status
  by_sample <- function(col) if (is.data.frame(st) && nrow(st) && col %in% colnames(st)) stats::setNames(st[[col]], st$sample) else NULL
  failed <- if (is.data.frame(st) && nrow(st)) st[st$status %in% "failed", , drop = FALSE] else NULL
  failure_msg <- if (!is.null(failed) && nrow(failed))
    stats::setNames(as.list(paste0("CNV calling failed: ", failed$message)), failed$sample) else NULL

  # ---- extra alerts ----
  extra_alerts <- list()
  if (identical(engine, "legacy_canoes"))
    extra_alerts[[length(extra_alerts) + 1L]] <- list(kind = "method", severity = "info", title = "Legacy CANOES HMM engine used",
      text = "Reference selection, NNLS weighting and variance estimation are the same as the modern engine; only the HMM transition / emission / decoding / Phred maths differ.", tab = "run", samples = NULL)
  if (nrow(calls)) {
    sparse <- which(as.numeric(calls$KB) >= 5000 & as.integer(calls$NUM_TARG) <= 20)
    if (length(sparse))
      extra_alerts[[length(extra_alerts) + 1L]] <- list(kind = "calls", severity = "warning",
        title = sprintf("%d call(s) span > 5 Mb over <= 20 targets", length(sparse)),
        text = "A very wide call over a few, far-apart targets usually reflects a sparse panel region rather than a real event: check the call viewer.",
        tab = "cnvs", samples = rep_arr(unique(as.character(calls$SAMPLE[sparse]))))
  }

  # exon labels come from the numbering recorded BEFORE QC removed any target
  ann <- if (!is.null(rd$raw_bed) && all(c("target", "exon_number") %in% colnames(rd$raw_bed))) rd$raw_bed else bed_file
  exon_of <- function(id) ann$exon_number[match(as.integer(id), as.integer(ann$target))]
  gene_gap <- as.numeric(settings$plot_gene_gap %||% run$params$plot_gene_gap %||% 1)
  labels <- NULL

  prm <- run$params
  all_kv <- utils::modifyList(if (length(prm)) as.list(prm) else list(),
                              settings[!grepl("^report_", names(settings))])
  meta <- list(engine = engine, mode_label = mode_label, modechrom_label = mode_label, modechrom = modechrom,
               engine_label = if (identical(engine, "legacy_canoes")) "Legacy CANOES (original HMM)" else "CANOPE (modern HMM engine)",
               genome_build = genome, genome_version = if (identical(genome, "unknown")) "" else genome)
  if (!is.null(run$canope_version)) meta$version <- run$canope_version

  list(
    source = "CANOPE", pkg = "CANOPE", title = title %||% "CANOPE report", meta = meta,
    modechrom = modechrom,
    analysed_class = canope_analysed_classes(modechrom, rd$removeY),
    scope = function(input_chr, analysed_chr) canope_scope(input_chr, modechrom, rd$removeY, analysed_chromosomes = analysed_chr),
    scope_hints = CANOPE_SCOPE_HINTS,
    thresholds = thresholds, conf_rules = conf_rules, scoring = scoring,
    cols = list(sample = "SAMPLE", type = "CNV", confidence = "Confidence", id = "CNV.ID"),
    ref_cor = by_sample("ref_cor_median"), n_targets = by_sample("n_targets"),
    failure_msg = failure_msg,
    extra_qc = canope_sample_flags(rd), extra_alerts = extra_alerts,
    chrom_depth = run$chrom_depth, sex_check = run$sex_check,
    settings_kv = report_settings_kv(all_kv, modechrom),
    params = if (length(prm)) lapply(prm, function(x) if (length(x) > 1) rep_arr(x) else x) else NULL,
    order = function(calls) {
      rank <- c(HIGH = 3L, MEDIUM = 2L, LOW = 1L)
      r <- unname(rank[toupper(as.character(calls$Confidence))]); r[is.na(r)] <- 0L
      q <- as.numeric(calls$Q_SOME); q[is.na(q)] <- -1
      order(-r, -q, calls$CNV.ID)
    },
    window = function(calls, i) {
      if (is.null(labels)) labels <<- target_axis_labels(bed_file, counts)
      pd <- canope_plot_data(calls[i, , drop = FALSE], counts, bed_file, rd$models, rd$refs,
                             gene_gap = gene_gap, flank = as.integer(flank), labels = labels)
      if (is.null(pd)) NULL else window_from_plot_data(pd)
    },
    call_table = function(calls, call_type, calib, has_plot) {
      col <- function(nm, default = NA) if (nm %in% colnames(calls)) calls[[nm]] else rep(default, nrow(calls))
      tg <- lapply(as.character(calls$TARGETS), parse_canope_targets)
      first <- vapply(tg, function(t) if (length(t) >= 2L) t[1L] else NA_integer_, integer(1))
      last  <- vapply(tg, function(t) if (length(t) >= 2L) t[2L] else NA_integer_, integer(1))
      e1 <- exon_of(first); e2 <- exon_of(last)
      exons <- ifelse(is.na(e1) | is.na(e2), "", ifelse(e1 == e2, as.character(e1), paste0(e1, "-", e2)))
      kb <- round(as.numeric(col("KB")), 2); nt <- as.integer(col("NUM_TARG"))
      conf <- toupper(as.character(col("Confidence"))); conf[is.na(conf)] <- "NA"
      data.frame(
        id = calls$CNV.ID, sample = as.character(calls$SAMPLE), chr = format_chr_label(calls$CHR),
        start = as.numeric(calls$Start), end = as.numeric(calls$End),
        gene = as.character(col("GENE", "Unknown")), type = call_type, interval = as.character(col("INTERVAL", "")),
        exons = exons, nexons = nt, num_targ = nt,
        ratio = as.numeric(col("READS_RATIO")), bf = NA_real_,
        q_some = as.integer(col("Q_SOME")), mlcn = as.integer(col("MLCN")), cn_label = as.character(col("CN_label")),
        kb = kb, num_refs = as.integer(col("NUM_REFS")),
        expected = as.numeric(col("READS_EXPECTED")), observed = as.numeric(col("READS_OBSERVED")),
        confidence = conf, correlation = as.numeric(col("REF_CORR")),
        n_comp = as.integer(col("NUM_REFS")), comparators = gsub(";", ",", as.character(col("REF_SAMPLES", "")), fixed = TRUE),
        calibration = calib, genomic = as.character(col("INTERVAL", "")),
        gidx = as.numeric(col("global_start")),
        carriers = as.integer(col("Carriers", NA_integer_)), recurrent = col("Recurrent", FALSE) %in% TRUE,
        single_target = nt == 1L, sparse_span = !is.na(kb) & kb >= 5000 & !is.na(nt) & nt <= 20L,
        targets = as.character(col("TARGETS", "")),
        has_plot = has_plot, stringsAsFactors = FALSE, check.names = FALSE)
    })
}
