# ============================================================================
# CANOPE static (PDF) plots.
#
# CANOPE only builds the *window* of a call (its negative-binomial model, in
# log2 units); the four panels themselves (coverage, gene tiles, 95% interval,
# z-scores) are the shared ones of cnv_shared_plots.R, the same code and the
# same colours ECHO and the HTML report use.
# ============================================================================

#' Look up per-target model statistics
#' @noRd
model_lookup <- function(model, field, target_ids) {
  if (is.list(model) && !is.null(model$target)) {
    idx <- match(as.integer(target_ids), as.integer(model$target))
    vals <- model[[field]][idx]
  } else if (!is.null(names(model[[field]]))) {
    vals <- model[[field]][as.character(target_ids)]
  } else {
    vals <- model[[field]][target_ids]
  }
  if (any(is.na(vals))) {
    warning("Some targets missing model statistics; using median imputation")
    vals[is.na(vals)] <- stats::median(model[[field]], na.rm = TRUE)
  }
  vals
}

#' Look up per-target model statistics (matrix form)
#' @noRd
model_matrix_lookup <- function(model, field, target_ids) {
  mat <- model[[field]]
  idx <- match(as.integer(target_ids), as.integer(model$target))
  out <- mat[idx, , drop = FALSE]
  if (any(is.na(idx))) {
    warning("Some targets missing model statistics; using median imputation")
    for (j in seq_len(ncol(out))) {
      miss <- is.na(out[, j])
      if (any(miss)) out[miss, j] <- stats::median(mat[, j], na.rm = TRUE)
    }
  }
  out
}

#' Target annotation of a CANOPE workspace in the shape the shared plots expect
#'
#' Adds the lower-case \code{gene}, the \code{off_target} flag and the exon
#' numbers recorded BEFORE QC (so a removed exon does not renumber its
#' neighbours).
#' @noRd
canope_target_annotation <- function(bed_file, counts = NULL) {
  if (!"target" %in% colnames(bed_file))
    bed_file$target <- if (!is.null(counts) && "target" %in% colnames(counts)) counts$target else seq_len(nrow(bed_file))
  if (!any(tolower(colnames(bed_file)) %in% c("gene", "gene_symbol", "symbol"))) bed_file$GENE <- "TargetRegion"
  bed_file <- normalize_target_annotation(bed_file)
  if (!"exon_number" %in% colnames(bed_file) || all(is.na(bed_file$exon_number)))
    bed_file$exon_number <- compute_exon_index(bed_file)
  bed_file
}

#' Plot window of one CANOPE call (feeds the PDF panels AND the HTML report)
#'
#' All coverage values are log2(count + 0.5); the ratio is
#' log2(observed / expected) with a 95\% negative-binomial predictive interval
#' from the HMM's own per-target mean / variance, and the z-scores use the same
#' variance, so a well-calibrated region shows the reference points clustered
#' around zero.
#'
#' @param call_row One row of the CANOPE call table (\code{SAMPLE, CNV, CHR,
#'   TARGETS}).
#' @param counts,bed_file Count matrix and its target annotation
#'   (\code{canope_target_annotation()}).
#' @param models,refs \code{models} and \code{refs} of the workspace.
#' @param gene_gap Extra x units between genes.
#' @param flank Targets shown on each side of the call.
#' @param labels X-axis label of every target (\code{target_axis_labels()}).
#' @return A window list (see cnv_shared_plots.R) or \code{NULL} when the call
#'   cannot be drawn (no model, no references, empty window).
#' @noRd
canope_plot_data <- function(call_row, counts, bed_file, models, refs, gene_gap = 1, flank = 5L,
                             labels = target_axis_labels(bed_file, counts)) {
  sample <- as.character(call_row$SAMPLE)
  model <- models[[sample]]
  ref_samples <- refs[[sample]]
  if (!length(ref_samples) || is.null(model) || !sample %in% colnames(counts)) return(NULL)
  tb <- parse_canope_targets(call_row$TARGETS)
  if (length(tb) < 2L) return(NULL)
  idx_start <- tb[1L]; idx_end <- tb[2L]
  target_chr <- normalize_chromosome_vec(call_row$CHR, bed_file$chromosome)
  bed_chr <- normalize_chromosome_vec(bed_file$chromosome, bed_file$chromosome)

  # rows of the call = every retained target whose id lies in the call's range
  # (ids dropped by QC are simply absent, not "unmatched")
  cnv_rows <- which(bed_file$target >= idx_start & bed_file$target <= idx_end & bed_chr == target_chr)
  if (!length(cnv_rows)) return(NULL)
  full <- seq(max(1L, min(cnv_rows) - flank), min(nrow(bed_file), max(cnv_rows) + flank))
  single_chr <- length(unique(bed_chr[full])) == 1L
  prev <- FALSE; new_chr <- ""; er <- full
  if (!single_chr) {
    prev    <- bed_chr[full[1L]] != target_chr
    new_chr <- bed_file$chromosome[if (prev) full[1L] else utils::tail(full, 1L)]
    er      <- full[bed_chr[full] == target_chr]
  }
  # Only the targets the HMM evaluated for this sample: a target dropped for this
  # sample (low reference mean, near-zero reads) has no model statistics and
  # would otherwise be drawn with imputed values.
  if (!is.null(model$target)) er <- er[bed_file$target[er] %in% model$target]
  if (!length(er)) return(NULL)

  gp <- compute_gene_gap_positions(bed_file, er, gap = gene_gap)
  tid <- bed_file$target[er]
  affected <- er %in% cnv_rows
  pc <- 0.5

  test_counts <- as.numeric(model_lookup(model, "test_counts", tid))
  mu   <- as.numeric(model_lookup(model, "mean", tid))
  vr   <- as.numeric(model_lookup(model, "var_estimate", tid))
  expected  <- pmax(mu, 1)
  disp <- vr - mu
  size <- ifelse(disp > 0, mu^2 / disp, Inf)
  ratio <- log2((test_counts + pc) / (expected + pc))
  lo <- log2((stats::qnbinom(0.025, mu = mu, size = size) + pc) / (expected + pc))
  hi <- log2((stats::qnbinom(0.975, mu = mu, size = size) + pc) / (expected + pc))

  # coverage of the test sample and of the references, the latter rescaled to the test depth
  cnt <- function(s) as.numeric(counts[er, s])
  tcnt <- cnt(sample)
  have <- ref_samples %in% colnames(counts)
  ref_mat <- model_matrix_lookup(model, "ref_matrix", tid)
  if (ncol(ref_mat) != length(ref_samples)) {          # older workspace: keep what lines up
    k <- min(ncol(ref_mat), length(ref_samples)); ref_samples <- ref_samples[seq_len(k)]; have <- have[seq_len(k)]
    ref_mat <- ref_mat[, seq_len(k), drop = FALSE]
  }
  ref_samples <- ref_samples[have]; ref_mat <- ref_mat[, have, drop = FALSE]
  if (!length(ref_samples)) return(NULL)
  ref_log <- stats::setNames(lapply(ref_samples, function(r) {
    v <- cnt(r)
    log2(v * stats::median(tcnt / pmax(v, 1), na.rm = TRUE) + pc)
  }), ref_samples)

  sd_z <- sqrt(pmax(vr, 1)); sd_z[!is.finite(sd_z) | sd_z < 1e-6] <- 1e-6
  rmean <- rowMeans(ref_mat, na.rm = TRUE)
  ref_z <- stats::setNames(lapply(seq_along(ref_samples), function(i) (ref_mat[, i] - rmean) / sd_z), ref_samples)

  bg <- check_background_calibration(ratio, lo, hi, affected)
  subtitle <- if (isTRUE(bg$flag))
    sprintf("%d%% of background exons outside CI (%d/%d) — check region/reference match",
            round(bg$pct_outside), bg$n_outside, bg$n_background) else NULL

  list(sample = sample, exon_range = er, px = gp$px, gene_group = gp$gene_group,
       gene = as.character(bed_file$gene[er]),
       off = as.logical(bed_file$off_target[er] %||% rep(FALSE, length(er))),
       xlab = labels[er], affected = affected,
       test_log = log2(tcnt + pc), ref_log = ref_log,
       ratio = ratio, lo = lo, hi = hi, observed = test_counts, expected = expected,
       z = (test_counts - rmean) / sd_z, ref_z = ref_z,
       bg_calibration = bg, ci_subtitle = subtitle,
       call_type = toupper(substr(as.character(call_row$CNV), 1, 3)),
       baseline = 0,
       labels = list(cov = "log2 coverage", ratio = "log2(observed / expected)",
                     z = "Z-score vs references (NB SD)"),
       single_chr = single_chr, prev = prev, new_chr = new_chr,
       idx_start = idx_start, idx_end = idx_end)
}

#' Generate CANOPE CNV detection plots (one PDF per call)
#'
#' Each PDF has the four panels ECHO draws, from the same code: log2 coverage
#' of the test sample against its references, the gene tiles, the
#' log2(observed / expected) with the 95\% negative-binomial predictive
#' interval, and the z-scores relative to the references.
#'
#' @param rdata_file Workspace \code{.RData} written by \code{run_canope()}.
#' @param output_dir Directory for the PDFs (one sub-folder per sample).
#' @param modechrom \code{"A"}, \code{"XX"} or \code{"XY"} (restricts the calls plotted).
#' @param prefix Filename prefix (files are named \code{CANOPE_<prefix>_<sample>_<gene>_<id>.pdf}).
#' @param flank Targets shown on each side of the call.
#' @param gene_gap Extra x-axis units between genes (\code{0} = none).
#' @param log_file Optional pipeline log (the plot messages are appended).
#' @param debug Also log each PDF written.
#' @return Invisibly, the number of PDFs written.
#' @importFrom stats qnbinom median
#' @export
generate_plots <- function(rdata_file, output_dir = "./plots", modechrom = "A", prefix = NULL,
                           flank = 5L, gene_gap = 1, log_file = NULL, debug = FALSE) {
  lg <- cnv_logger(log_file, tag = "CANOPE", append = TRUE)
  lg$info("BEGIN plot generation")
  objs <- load_rdata(rdata_file, required = c("cnv_calls", "counts", "models", "refs"))
  cnv_calls <- objs$cnv_calls; counts <- objs$counts; models <- objs$models; refs <- objs$refs
  if (length(models) && any(vapply(models, function(m) !all(c("target", "test_counts", "ref_matrix") %in% names(m)), logical(1))))
    stop("[ERROR] 'models' in this RData lack 'test_counts'/'ref_matrix' (older run). Re-run run_canope() to regenerate the workspace.")
  if (is.null(cnv_calls) || nrow(cnv_calls) == 0) { lg$info("No CNV calls to plot."); return(invisible(0L)) }

  bed_file <- canope_target_annotation(objs$bed_file %||% counts[, intersect(c("chromosome", "start", "end", "target", "gc", "GENE"), colnames(counts))], counts)
  labels <- target_axis_labels(bed_file)
  # "CANOPE_<prefix>" whether or not the caller already put the CANOPE_ in front
  base <- sub("^CANOPE_?", "", as.character(prefix %||% ""))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  keep <- canope_mode_mask(cnv_calls$CHR, modechrom, isTRUE(objs$removeY %||% TRUE))
  cnv_plot <- cnv_calls[keep, , drop = FALSE]
  lg$info(sprintf("Number of CNV calls to plot: %d", nrow(cnv_plot)))
  if (!nrow(cnv_plot)) { lg$info("No CNV calls remain after chromosome filtering."); return(invisible(0L)) }

  n_written <- 0L
  for (i in seq_len(nrow(cnv_plot))) {
    tryCatch({
      row <- cnv_plot[i, , drop = FALSE]
      pd <- canope_plot_data(row, counts, bed_file, models, refs, gene_gap = gene_gap, flank = flank, labels = labels)
      if (is.null(pd)) {
        lg$warn(sprintf("Skipping plot for %s (call %s): no model, no references or empty window", row$SAMPLE, row$CNV.ID))
        next
      }
      gene_str   <- substr(sanitize_filename(as.character(row$GENE)), 1, 60)
      sample_dir <- file.path(output_dir, sanitize_filename(pd$sample))
      dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
      file_path <- file.path(sample_dir, paste0(paste(c("CANOPE", if (nzchar(base)) base, pd$sample, gene_str, row$CNV.ID), collapse = "_"), ".pdf"))
      save_cnv_pdf(list(create_coverage_plot(pd), create_gene_tile_plot(pd),
                        create_ci_plot(pd), create_zscore_plot(pd)), file_path)
      n_written <- n_written + 1L
      if (debug) lg$info(paste("Saved plot:", file_path))
    }, error = function(e) lg$error(sprintf("Failed to plot call %s: %s", cnv_plot$CNV.ID[i], conditionMessage(e))))
  }
  lg$info(sprintf("END plot generation (%d/%d PDF(s) written)", n_written, nrow(cnv_plot)))
  invisible(n_written)
}
