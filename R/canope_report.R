#' Generate the interactive, multi-sample HTML report for CANOPE (unified).
#'
#' Builds one self-contained HTML file from a CANOPE workspace, using the
#' shared renderer so the output is structurally identical to the ECHO report.
#'
#' @param rdata_output Path to the CANOPE workspace .RData.
#' @param qc_metrics_file Optional path to the QC-metrics TSV.
#' @param output_dir Directory for the HTML (created if missing).
#' @param settings Optional named list of thresholds.
#' @param sample_table Optional path to a sample sheet.
#' @param log_file Optional pipeline log.
#' @param prefix Output name prefix.
#' @param template_path Deprecated and ignored.
#' @param max_call_plots Maximum number of calls whose interactive plots are embedded.
#' @param max_qc_rows Row cap for the embedded exon-level QC flag table.
#' @param flank Number of flanking exons shown each side of a call.
#' @param title Report title.
#' @param plotly_source "inline" (default) or "cdn".
#' @return Invisibly returns the path to the HTML file.
#' @export
generate_canope_report <- function(rdata_output,
                                   qc_metrics_file = NULL,
                                   output_dir,
                                   settings = NULL,
                                   sample_table = NULL,
                                   log_file = NULL,
                                   prefix = "CANOPE",
                                   template_path = NULL,
                                   max_call_plots = 300L,
                                   max_qc_rows = 3000L,
                                   flank = 5L,
                                   title = NULL,
                                   plotly_source = c("inline", "cdn")) {

  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required.")
  if (!is.null(template_path)) message("[INFO] 'template_path' is ignored.")
  if (!file.exists(rdata_output)) stop("[ERROR] Workspace not found: ", rdata_output)
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE, showWarnings = FALSE); if (!dir.exists(output_dir)) stop("Failed to create output directory: ", output_dir) }

  env <- new.env()
  load(rdata_output, envir = env)
  required <- c("cnv_calls", "counts", "bed_file", "models", "refs")
  missing_vars <- setdiff(required, ls(env))
  if (length(missing_vars) > 0) stop("[ERROR] Missing variables in RData: ", paste(missing_vars, collapse = ", "))

  settings <- settings %||% list()
  settings$report_plotly <- match.arg(plotly_source)

  data <- canope_build_data(env, qc_metrics_file, settings, sample_table, log_file,
                            prefix, title, max_call_plots, max_qc_rows, flank)

  base <- sub("^CANOPE_?", "", prefix)
  out_file <- file.path(normalizePath(output_dir, mustWork = TRUE),
                        if (nzchar(base)) paste0("CANOPE_", base, "_report.html") else "CANOPE_report.html")
  report_html(data, settings, out_file)
  invisible(out_file)
}

# ----------------------------------------------------------------------------
# CANOPE data builder -- produces the UNIFIED schema expected by report_assets.R
# ----------------------------------------------------------------------------

canope_meta_cols <- c("target", "chromosome", "start", "end", "gc", "GENE", "mean", "exon_number")

canope_arr <- function(x) I(as.vector(x))
canope_num <- function(x, digits = 3) canope_arr(round(as.numeric(x), digits))

canope_chr_order <- function(chr) {
  s <- toupper(sub("^chr", "", as.character(chr), ignore.case = TRUE))
  n <- suppressWarnings(as.numeric(s))
  n[s == "X"] <- 23; n[s == "Y"] <- 24; n[s %in% c("M", "MT")] <- 25; n[is.na(n)] <- 99; n
}
canope_conf_rank <- function(x) { r <- c(HIGH = 3L, MEDIUM = 2L, LOW = 1L); out <- unname(r[toupper(as.character(x))]); out[is.na(out)] <- 0L; out }

canope_build_data <- function(env, qc_metrics_file, settings, sample_table, log_file,
                              prefix, title, max_call_plots, max_qc_rows, flank) {
  cnv <- env$cnv_calls; counts <- env$counts; bed <- env$bed_file; models <- env$models; refs <- env$refs
  run <- if ("run_info" %in% ls(env)) env$run_info else list()
  engine <- if ("hmm_engine_used" %in% ls(env)) env$hmm_engine_used else (run$engine %||% "new")
  settings <- settings %||% list()
  set <- function(name, default) settings[[name]] %||% run$thresholds[[name]] %||% default

  if (is.null(cnv) || ncol(cnv) == 0L) {
    cnv <- data.frame(
      SAMPLE = character(), CHR = character(),
      INTERVAL = character(), GENE = character(),
      CNV = character(), KB = numeric(),
      NUM_TARG = integer(), MLCN = integer(),
      CN_label = character(), Q_SOME = integer(),
      NUM_REFS = integer(), Confidence = character(),
      TARGETS = character(), stringsAsFactors = FALSE)
  }
  if (!"Confidence" %in% colnames(cnv)) cnv$Confidence <- rep(NA_character_, nrow(cnv))
  if (!"CN_label" %in% colnames(cnv)) cnv$CN_label <- rep(NA_character_, nrow(cnv))

  thr <- list(
    qc_min_cov = as.numeric(set("qc_min_cov", 100)),
    min_cov = as.numeric(set("qc_min_cov", 100)),
    qc_min_total_reads = as.numeric(set("qc_min_total_reads", 3e5)),
    min_total_reads = as.numeric(set("qc_min_total_reads", 3e5)),
    qc_min_corr = as.numeric(set("qc_min_corr", 0.98)),
    min_corr = as.numeric(set("qc_min_corr", 0.98)),
    qc_max_exon_cv = as.numeric(set("qc_max_exon_cv", 0.5)),
    max_exon_cv = as.numeric(set("qc_max_exon_cv", 0.5)),
    qc_zscore = as.numeric(run$params$qc_zscore %||% 3))
  gene_gap <- as.numeric(settings$plot_gene_gap %||% 1)

  # ---- samples ----
  sample_cols <- setdiff(colnames(counts), canope_meta_cols)
  requested <- unique(c(run$samples_requested, sample_cols, as.character(cnv$SAMPLE)))
  requested <- requested[!is.na(requested) & nzchar(requested)]
  depth <- reads <- setNames(rep(NA_real_, length(requested)), requested)
  if (length(sample_cols) > 0) {
    depth[sample_cols] <- vapply(sample_cols, function(s) stats::median(counts[[s]], na.rm = TRUE), numeric(1))
    reads[sample_cols] <- vapply(sample_cols, function(s) sum(counts[[s]], na.rm = TRUE), numeric(1))
  }

  corr <- NULL; max_corr <- setNames(rep(NA_real_, length(requested)), requested)
  if (length(sample_cols) >= 2L && length(sample_cols) <= 400L) {
    mat <- as.matrix(counts[, sample_cols, drop = FALSE]); storage.mode(mat) <- "numeric"
    cm <- tryCatch(suppressWarnings(stats::cor(mat, method = "spearman", use = if (anyNA(mat)) "pairwise.complete.obs" else "everything")), error = function(e) NULL)
    if (!is.null(cm)) {
      cm2 <- cm; diag(cm2) <- NA_real_
      max_corr[sample_cols] <- apply(cm2, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
      ord <- if (length(sample_cols) >= 3L && !anyNA(cm)) tryCatch(stats::hclust(stats::as.dist(1 - cm), method = "average")$order, error = function(e) seq_along(sample_cols)) else seq_along(sample_cols)
      if (length(sample_cols) <= 150L) corr <- list(names = canope_arr(sample_cols[ord]), z = lapply(ord, function(i) canope_num(cm[i, ord], 3)))
    }
  }

  noise <- run$sample_noise
  if (is.null(noise) && length(sample_cols) >= 2L) {
    noise <- tryCatch(
      detect_outlier_samples(counts[, sample_cols, drop = FALSE],
                            z_threshold = thr$qc_zscore),
      error = function(e) NULL)
  }
  status_df <- run$sample_status
  get_status <- function(s, col, default = NA) {
    if (is.null(status_df) || nrow(status_df) == 0 || !s %in% status_df$sample)
      return(default)
    v <- status_df[[col]][match(s, status_df$sample)]
    if (is.null(v)) default else v
  }

  n_by <- function(s, f) sum(cnv$SAMPLE == s & f)
  conf <- toupper(as.character(cnv$Confidence))
  sex_check <- run$sex_check

  # ---- per-sample flags ----
  flags <- setNames(vector("list", length(requested)), requested)
  add_flag <- function(s, level, text) flags[[s]][[length(flags[[s]]) + 1L]] <<- list(level = level, text = text)
  call_counts <- vapply(requested, function(s) sum(cnv$SAMPLE == s), numeric(1))
  analysed_calls <- call_counts[requested %in% (run$samples_analysed %||% sample_cols)]
  mad_calls <- if (length(analysed_calls) >= 5) stats::mad(analysed_calls) else NA_real_
  med_calls <- if (length(analysed_calls) >= 5) stats::median(analysed_calls) else NA_real_
  many_calls_cut <- if (!is.na(mad_calls)) max(10, med_calls + 3 * max(mad_calls, 1)) else NA_real_
  low_genes <- as.character(unlist(settings$low_confidence_genes %||% eval(formals(score_canope_confidence)$low_confidence_genes)))
  lcr <- run$low_coverage_runs

  for (s in requested) {
    st <- get_status(s, "status", if (s %in% sample_cols) "called" else "unknown")
    msg <- get_status(s, "message", "")
    if (identical(st, "failed")) add_flag(s, "fail", paste0("CNV calling FAILED (", msg, ")."))
    if (identical(st, "excluded")) add_flag(s, "fail", paste0("Excluded before calling: ", msg, "."))
    if (!is.na(depth[[s]]) && depth[[s]] < thr$qc_min_cov) add_flag(s, "warn", sprintf("Low depth: %.1f (< %g).", depth[[s]], thr$qc_min_cov))
    if (!is.na(reads[[s]]) && reads[[s]] < thr$qc_min_total_reads)
      add_flag(s, "warn", sprintf("Low reads: %s (< %s).",
              format(round(reads[[s]]), big.mark = ","),
              format(thr$qc_min_total_reads, big.mark = ",")))
    if (!is.na(max_corr[[s]]) && max_corr[[s]] < thr$qc_min_corr)
      add_flag(s, "warn", sprintf("Low correlation: %.3f (< %g).",
              max_corr[[s]], thr$qc_min_corr))
    if (!is.null(noise) && s %in% noise$sample) {
      z <- noise$robust_z[match(s, noise$sample)]
      if (!is.na(z) && z > thr$qc_zscore && !identical(st, "excluded"))
        add_flag(s, "warn", sprintf("Noisy (z = %.2f > %g).", z, thr$qc_zscore))
    }
    nr <- get_status(s, "n_refs", NA_integer_)
    if (!is.na(nr) && nr < 5)
      add_flag(s, "warn", sprintf("Only %d refs.", as.integer(nr)))
    if (isFALSE(get_status(s, "genotyped", NA)))
      add_flag(s, "fail", "Genotyping skipped.")
    if (!is.na(many_calls_cut) && call_counts[[s]] > many_calls_cut)
      add_flag(s, "warn", sprintf("Many calls (%d vs median %g).",
              as.integer(call_counts[[s]]), med_calls))
    if (!is.null(lcr) && s %in% lcr$SAMPLE) {
      r <- lcr[lcr$SAMPLE == s, , drop = FALSE]
      add_flag(s, "warn", sprintf("%d low-coverage run(s) removed.", nrow(r)))
    }
    if (!is.null(sex_check) && s %in% sex_check$sample) {
      k <- match(s, sex_check$sample)
      if (isTRUE(sex_check$mismatch[k]))
        add_flag(s, "warn", sprintf("Sex mismatch (annotated %s, looks %s).",
                sex_check$annotated_sex[k], sex_check$inferred_sex[k]))
    }
  }
  level_of <- function(f) {
    if (any(vapply(f, function(x) x$level == "fail", logical(1)))) "fail"
    else if (any(vapply(f, function(x) x$level == "warn", logical(1)))) "warn"
    else "ok"
  }
  sample_level <- vapply(flags, function(f) if (length(f)) level_of(f) else "ok", character(1))

  st_chr <- vapply(requested, function(s) as.character(get_status(s, "status", if (s %in% sample_cols) "called" else "unknown")), character(1))
  samples_df <- data.frame(
    sample = requested, status = st_chr, level = unname(sample_level),
    n_cnvs = as.integer(call_counts),
    n_del = vapply(requested, function(s) n_by(s, cnv$CNV == "DEL"), numeric(1)),
    n_dup = vapply(requested, function(s) n_by(s, cnv$CNV == "DUP"), numeric(1)),
    n_high = vapply(requested, function(s) n_by(s, conf == "HIGH"), numeric(1)),
    n_medium = vapply(requested, function(s) n_by(s, conf == "MEDIUM"), numeric(1)),
    n_low = vapply(requested, function(s) n_by(s, conf == "LOW"), numeric(1)),
    median_depth = round(as.numeric(depth), 1), total_reads = round(as.numeric(reads)),
    max_corr = round(as.numeric(max_corr), 4),
    noise_z = if (!is.null(noise)) round(noise$robust_z[match(requested, noise$sample)], 2) else NA_real_,
    n_refs = vapply(requested, function(s) as.numeric(get_status(s, "n_refs", length(refs[[s]]) %||% NA)), numeric(1)),
    ref_cor = round(vapply(requested, function(s) as.numeric(get_status(s, "ref_cor_median", NA_real_)), numeric(1)), 3),
    n_targets = vapply(requested, function(s) as.numeric(get_status(s, "n_targets", NA_integer_)), numeric(1)),
    n_dropped = vapply(requested, function(s) as.numeric(get_status(s, "n_low_count_dropped", NA_integer_)), numeric(1)),
    sex_inferred = if (!is.null(sex_check)) sex_check$inferred_sex[match(requested, sex_check$sample)] else NA_character_,
    stringsAsFactors = FALSE, row.names = NULL)

  # ---- chromosomes ----
  if (!is.null(run$chromosome_status)) { chr_status <- run$chromosome_status; chr_messages <- as.character(run$chromosome_messages %||% character(0)) }
  else { an <- as.character(bed$chromosome); cs <- canope_chromosome_status(an, an, modechrom = "A", removeY = TRUE); chr_status <- cs$status; chr_messages <- cs$messages }
  chr_depth <- run$chrom_depth; chr_depth_json <- NULL
  if (!is.null(chr_depth) && !is.null(dim(chr_depth)) && ncol(chr_depth) > 0) {
    o <- order(canope_chr_order(rownames(chr_depth))); chr_depth <- chr_depth[o, , drop = FALSE]
    chr_depth_json <- list(chromosomes = canope_arr(rownames(chr_depth)), samples = canope_arr(colnames(chr_depth)), z = lapply(seq_len(nrow(chr_depth)), function(i) canope_num(chr_depth[i, ], 3)))
  }
  inv <- run$chromosome_inventory; analysed_chr <- unique(as.character(bed$chromosome))
  if (!is.null(inv)) {
    inv_df <- data.frame(
      chromosome = inv$chromosome, class = inv$class,
      targets_in_bed = inv$n_targets_input,
      targets_analysed = vapply(inv$chromosome, function(ch)
        sum(as.character(bed$chromosome) == ch), numeric(1)),
      analysed = inv$chromosome %in% analysed_chr,
      stringsAsFactors = FALSE)
  } else {
    tb <- table(as.character(bed$chromosome))
    inv_df <- data.frame(
      chromosome = names(tb),
      class = canope_chrom_class(names(tb)),
      targets_in_bed = as.integer(tb),
      targets_analysed = as.integer(tb),
      analysed = TRUE, stringsAsFactors = FALSE)
  }
  inv_df <- inv_df[order(canope_chr_order(inv_df$chromosome)), , drop = FALSE]
  if (!is.null(run$chromosome_status)) {
    ca <- chr_status$class
    both_ok <- all(chr_status$analysed[ca %in% c("autosome", "X")] |
                   !chr_status$in_input[ca %in% c("autosome", "X")]) &&
              any(chr_status$analysed[ca == "autosome"]) &&
              any(chr_status$analysed[ca == "X"])
  } else {
    both_ok <- NA
  }
  ext <- tapply(as.numeric(bed$end), format_chr_label(bed$chromosome), max); ext <- ext[order(canope_chr_order(names(ext)))]

  # ---- QC metrics ----
  qc <- NULL
  if (!is.null(qc_metrics_file) && file.exists(qc_metrics_file)) {
    qc <- tryCatch(data.table::fread(qc_metrics_file, header = TRUE,
                                     fill = TRUE, data.table = FALSE),
                   error = function(e) NULL)
  }
  qc_json <- NULL; qc_flags_df <- NULL
  if (!is.null(qc) && nrow(qc) > 0) {
    type_n <- as.data.frame(table(Type = qc$Type), stringsAsFactors = FALSE)
    keep_first <- qc$Type %in% c("Whole sample", "Missing Chromosome", "Chromosome not analysed")
    exon_rows <- which(!keep_first)
    capped <- length(exon_rows) > max_qc_rows
    rows_idx <- c(which(keep_first), utils::head(exon_rows, max_qc_rows))
    gene_n <- NULL
    ex <- qc[!keep_first, , drop = FALSE]
    if (nrow(ex) > 0) {
      g <- as.data.frame(table(Gene = ex$Gene, Type = ex$Type), stringsAsFactors = FALSE)
      g <- g[g$Freq > 0, , drop = FALSE]
      g <- g[order(-g$Freq), , drop = FALSE]
      gene_n <- rep_table(utils::head(g, 60))
    }
    qc_json <- list(types = rep_table(type_n), rows = rep_table(qc[rows_idx, , drop = FALSE]), n_total = nrow(qc), capped = capped, genes = gene_n)
    # 'Severity' is optional in the QC metrics TSV (run_canope_qc_metrics()
    # does not write it): default it *before* the column select, which
    # otherwise fails with "undefined columns selected" and aborts the report.
    if (!"Severity" %in% colnames(qc)) {
      qc$Severity <- ifelse(qc$Type %in% c("Whole exon", "Exon variability"), "info", "warning")
    }
    qc_flags_df <- qc[, c("Sample", "Exon", "Type", "Details", "Gene", "Severity"), drop = FALSE]
  } else {
    qc_flags_df <- data.frame(Sample = character(), Exon = character(),
                              Type = character(), Details = character(),
                              Gene = character(), Severity = character(),
                              stringsAsFactors = FALSE)
  }

  # ---- PCA ----
  pca_json <- NULL
  if (length(sample_cols) >= 3L) {
    pca_json <- tryCatch({
      lm <- log2(as.matrix(counts[, sample_cols, drop = FALSE]) + 1)
      rv <- apply(lm, 1, stats::var, na.rm = TRUE)
      lm <- lm[!is.na(rv) & rv > 0, , drop = FALSE]
      if (nrow(lm) < 2L) {
        NULL
      } else {
        pc <- stats::prcomp(t(lm), scale. = TRUE, center = TRUE)
        k <- min(5L, ncol(pc$x))
        ve <- summary(pc)$importance[2, seq_len(k)] * 100
        list(
          names = canope_arr(sample_cols),
          pc1 = canope_num(pc$x[, 1], 3),
          pc2 = canope_num(pc$x[, 2], 3),
          var1 = round(ve[1], 1), var2 = round(ve[2], 1),
          group = rep("All samples", length(sample_cols)),
          n_rows = nrow(lm), n_samples = length(sample_cols),
          samples = canope_arr(sample_cols),
          pcs = setNames(lapply(seq_len(k), function(i) canope_num(pc$x[, i], 3)),
                         paste0("PC", seq_len(k))),
          var = canope_num(ve, 1), n_targets = nrow(lm))
      }
    }, error = function(e) {
      message("[WARNING] PCA skipped: ", conditionMessage(e))
      NULL
    })
  }

  # ---- calls (lowercase unified column names) ----
  calls_json <- NULL; plots <- list(); n_flagged_bg <- 0L; n_plot_missing <- 0L
  if (nrow(cnv) > 0) {
    iv <- lapply(as.character(cnv$INTERVAL), parse_canope_interval)
    gene_hit <- vapply(as.character(cnv$GENE), function(g) any(trimws(unlist(strsplit(g, "[;, ]"))) %in% low_genes), logical(1))
    cdf <- data.frame(
      id = seq_len(nrow(cnv)), sample = as.character(cnv$SAMPLE), chr = format_chr_label(cnv$CHR),
      start = vapply(iv, function(p) as.numeric(p$start), numeric(1)), end = vapply(iv, function(p) as.numeric(p$end), numeric(1)),
      interval = as.character(cnv$INTERVAL), gene = as.character(cnv$GENE), type = as.character(cnv$CNV),
      kb = round(as.numeric(cnv$KB), 2), num_targ = as.integer(cnv$NUM_TARG), mlcn = as.integer(cnv$MLCN),
      cn_label = as.character(cnv$CN_label), q_some = as.integer(cnv$Q_SOME), num_refs = as.integer(cnv$NUM_REFS),
      confidence = ifelse(is.na(cnv$Confidence), "NA", toupper(as.character(cnv$Confidence))),
      homology_gene = gene_hit, targets = as.character(cnv$TARGETS),
      bf = NA_real_, ratio = NA_real_, observed = NA_real_, expected = NA_real_, correlation = NA_real_, n_comp = NA_integer_, comparators = "", calibration = "",
      stringsAsFactors = FALSE)
    ord <- order(-canope_conf_rank(cdf$confidence), -ifelse(is.na(cdf$q_some), -1L, cdf$q_some))
    embed <- utils::head(ord, max_call_plots)
    bed_chr <- format_chr_label(bed$chromosome); exon_index <- compute_exon_index(bed)
    cdf$bg_flag <- FALSE; cdf$has_plot <- FALSE
    for (i in embed) {
      w <- tryCatch(
        canope_call_window(cnv[i, , drop = FALSE], bed, bed_chr, models,
                          flank = flank, gene_gap = gene_gap,
                          exon_index = exon_index),
        error = function(e) NULL)
      if (!is.null(w)) {
        plots[[as.character(cdf$id[i])]] <- w
        cdf$has_plot[i] <- TRUE
        cdf$bg_flag[i] <- isTRUE(w$bg$flag)
        cdf$calibration[i] <- if (isTRUE(w$bg$flag)) "check" else "ok"
      }
    }
    n_flagged_bg <- sum(cdf$bg_flag); n_plot_missing <- sum(!cdf$has_plot)
    calls_json <- rep_table(cdf)
  }

  # ---- global warnings (CANOPE style) ----
  warns <- list()
  add_w <- function(level, category, text, sample = NULL) warns[[length(warns) + 1L]] <<- list(level = level, category = category, text = text, sample = sample)
  for (m in chr_messages) add_w("warn", "Chromosomes", m)
  if (isTRUE(both_ok == FALSE) || (!is.na(both_ok) && !both_ok)) add_w("fail", "Chromosomes", "Autosomes and gonosomes were NOT both analysed.")
  fl <- samples_df$sample[samples_df$status == "failed"]
  if (length(fl)) add_w("fail", "Samples", sprintf("%d sample(s) FAILED: %s.", length(fl), paste(fl, collapse = ", ")))
  ex <- samples_df$sample[samples_df$status == "excluded"]
  if (length(ex)) add_w("warn", "Samples", sprintf("%d excluded: %s.", length(ex), paste(ex, collapse = ", ")))
  for (s in requested) for (f in flags[[s]]) if (f$level %in% c("fail", "warn")) add_w(f$level, "Sample metrics", f$text, s)
  if (n_flagged_bg > 0) add_w("warn", "Calls", sprintf("%d call(s) background calibration flag.", n_flagged_bg))
  if (n_plot_missing > 0 && !is.null(calls_json)) add_w("info", "Calls", sprintf("%d call(s) no embedded plot.", n_plot_missing))
  if (identical(engine, "legacy_canoes")) add_w("info", "Method", "Legacy CANOES HMM used.")

  # ---- derive ECHO-style alerts from CANOPE warnings ----
  alert_from_warning <- function(w) {
    sev <- switch(w$level,
                  fail = "error", warn = "warning",
                  info = "info", "info")
    sm <- if (!is.null(w$sample)) I(w$sample) else NULL
    list(kind = w$category, severity = sev, title = w$category,
         text = w$text, tab = "qc", samples = sm)
  }
  alerts <- lapply(warns, alert_from_warning)

  # ---- derive ECHO-style scope from CANOPE chromosomes ----
  scope_df <- data.frame(
    chromosome = inv_df$chromosome, class = inv_df$class, label = inv_df$chromosome,
    targets = inv_df$targets_in_bed, first = NA_integer_, last = NA_integer_,
    present = TRUE, analysed = inv_df$analysed, stringsAsFactors = FALSE)

  # ---- log ----
  log_lines <- character(0)
  if (!is.null(log_file) && is.character(log_file) && file.exists(log_file)) log_lines <- readLines(log_file, warn = FALSE)
  log_warn <- log_lines[startsWith(log_lines, "[WARNING]")]
  si <- tryCatch(utils::capture.output(print(utils::sessionInfo())), error = function(e) character(0))

  prm <- run$params %||% list(); prm <- lapply(prm, function(x) if (length(x) > 1) canope_arr(x) else x)
  conf_rules <- list(
    high_q = settings$high_q_score %||% 50,
    med_q = settings$med_q_score %||% 20,
    high_refs = settings$high_num_refs %||% 10,
    med_refs = settings$med_num_refs %||% 5,
    high_targ = settings$high_num_targ %||% 1,
    med_targ = settings$med_num_targ %||% 1,
    genes = canope_arr(low_genes))

  mode_label <- switch(as.character(run$modechrom %||% "?"), A = "Autosomes (A)", XX = "chrX, female references (XX)", XY = "chrX/chrY, male references (XY)", "unknown")
  scope_incomplete <- isTRUE(both_ok == FALSE) || (!is.na(both_ok) && !both_ok)

  list(
    meta = list(
      source = "CANOPE",
      title = title %||% paste("CANOPE CNV report -", sub("^CANOPE_?", "", prefix) %||% prefix),
      prefix = sub("^CANOPE_?", "", prefix) %||% prefix, generated = format(Sys.time(), "%Y-%m-%d %H:%M"),
      version = run$canope_version %||% tryCatch(as.character(utils::packageVersion("CANOPE")), error = function(e) "dev"),
      engine = engine, engine_label = if (identical(engine, "legacy_canoes")) "Legacy CANOES (original HMM)" else "CANOPE (modern HMM engine)",
      genome_build = run$genome_build %||% "unknown", modechrom = run$modechrom %||% NA_character_, mode_label = mode_label, modechrom_label = mode_label,
      genome_version = run$genome_build %||% "", legacy_workspace = is.null(run$sample_status),
      focus = NULL, n_targets = nrow(bed), n_samples = length(requested), n_analysed = sum(samples_df$status == "called"),
      n_calls = nrow(cnv), scope_incomplete = scope_incomplete, qc_available = !is.null(qc_metrics_file) && file.exists(qc_metrics_file),
      windows_limited = n_plot_missing > 0, windows_limit = max_call_plots),
    thresholds = thr, conf_rules = conf_rules, params = prm,
    summary = list(
      samples_requested = length(requested), samples_called = sum(samples_df$status == "called"),
      samples_with_calls = length(unique(cnv$SAMPLE)), calls = nrow(cnv), targets = nrow(bed),
      both_chromosome_classes = if (is.na(both_ok)) NULL else both_ok),
    samples = rep_table(samples_df), flags = lapply(flags, function(f) unname(f)), sample_attrs = list(),
    warnings = unname(warns), alerts = alerts,
    scope = rep_table(scope_df), chr = rep_table(scope_df), chr_all = rep_table(scope_df),
    chromosomes = list(
      status = rep_table(chr_status),
      messages = canope_arr(chr_messages),
      inventory = rep_table(inv_df),
      depth = chr_depth_json,
      extent = list(chr = canope_arr(names(ext)),
                    len = canope_arr(as.numeric(ext))),
      sex = if (!is.null(sex_check)) rep_table(sex_check) else NULL),
    qc = qc_json, qc_flags = rep_table(qc_flags_df),
    corr = corr, pca = pca_json, depth_hist = NULL,
    calls = calls_json, plots = plots, windows = plots,
    refs = lapply(refs, function(x) I(as.character(x))),
    settings = lapply(seq_along(prm), function(i) I(c(names(prm)[i], prm[[i]]))),
    scoring = NULL,
    sample_table = if (!is.null(sample_table) && file.exists(sample_table)) {
      st_df <- utils::read.table(sample_table, header = TRUE, sep = "\t",
                                 stringsAsFactors = FALSE, check.names = FALSE,
                                 comment.char = "")
      rep_table(st_df)
    } else NULL,
    ref_bams = NULL,
    log = list(
      warnings = canope_arr(utils::head(log_warn, 300L)),
      tail = canope_arr(utils::tail(log_lines[!startsWith(log_lines, "[WARNING]")], 300L)),
      session = canope_arr(utils::head(si, 120L)),
      truncated = length(log_lines) > 600L),
    config_yaml = NULL)
}

# ----------------------------------------------------------------------------
# Per-call plot window (unchanged from original CANOPE report)
# ----------------------------------------------------------------------------

canope_call_window <- function(call_row, bed, bed_chr, models, flank = 5L, gene_gap = 1, exon_index = NULL) {
  s <- as.character(call_row$SAMPLE); m <- models[[s]]
  if (is.null(m) || is.null(m$target) || is.null(m$ref_matrix)) return(NULL)
  tb <- parse_canope_targets(call_row$TARGETS)
  if (length(tb) < 2L) return(NULL)
  chr_call <- format_chr_label(call_row$CHR)
  in_model <- bed$target %in% m$target
  call_rows <- which(bed$target >= tb[1] & bed$target <= tb[2] & bed_chr == chr_call & in_model)
  if (length(call_rows) == 0L) return(NULL)
  cand <- which(in_model & bed_chr == chr_call)
  lo_i <- match(min(call_rows), cand); hi_i <- match(max(call_rows), cand)
  win <- cand[seq(max(1L, lo_i - flank), min(length(cand), hi_i + flank))]
  mi <- match(bed$target[win], m$target)
  mean_e <- m$mean[mi]; var_e <- m$var_estimate[mi]; test <- m$test_counts[mi]
  refmat <- m$ref_matrix[mi, , drop = FALSE]
  pc <- 0.5; expected_safe <- pmax(mean_e, 1) + pc
  ratio <- log2((test + pc) / expected_safe)
  disp <- var_e - mean_e; size <- ifelse(disp > 0, mean_e^2 / disp, Inf)
  lo <- log2((stats::qnbinom(0.025, mu = mean_e, size = size) + pc) / expected_safe)
  hi <- log2((stats::qnbinom(0.975, mu = mean_e, size = size) + pc) / expected_safe)
  z_mean <- rowMeans(refmat, na.rm = TRUE); z_sd <- sqrt(pmax(var_e, 1)); z_sd[!is.finite(z_sd) | z_sd < 1e-6] <- 1e-6
  z <- (test - z_mean) / z_sd; ref_z <- (refmat - z_mean) / z_sd
  gp <- compute_gene_gap_positions(bed, win, gap = gene_gap)
  affected <- win %in% call_rows; bg <- check_background_calibration(ratio, lo, hi, affected)
  rows_by_ref <- function(mat, digits) lapply(seq_len(ncol(mat)), function(j) canope_num(mat[, j], digits))
  list(
    sample = s, exon = canope_arr(if (!is.null(exon_index)) exon_index[win] else seq_along(win)),
    gene = canope_arr(as.character(bed$GENE[win])), pos = canope_arr(paste0(bed_chr[win], ":", bed$start[win], "-", bed$end[win])),
    px = canope_num(gp$px, 2), gene_group = canope_arr(gp$gene_group), affected = canope_arr(as.integer(affected)),
    test_log = canope_num(log2(test + pc), 3), ratio = canope_num(ratio, 3), lo = canope_num(lo, 3), hi = canope_num(hi, 3), z = canope_num(z, 2),
    ref_names = canope_arr(colnames(refmat)), ref_log = rows_by_ref(log2(refmat + pc), 2), ref_z = rows_by_ref(ref_z, 2),
    bg = list(n_background = bg$n_background, n_outside = bg$n_outside, pct_outside = if (is.na(bg$pct_outside)) NULL else round(bg$pct_outside, 1), flag = isTRUE(bg$flag)),
    off = canope_arr(rep(FALSE, length(win))), xlab = canope_arr(as.character(seq_along(win))))
}

`%||%` <- function(a, b) if (is.null(a)) b else a