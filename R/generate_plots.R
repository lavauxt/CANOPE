#' Generate CANOPE CNV Detection Plots
#'
#' Produces per-call PDFs with reference/test log-coverage, gene tiles,
#' log2(observed/expected) with NB predictive intervals, and a z-score panel
#' relative to reference samples.
#'
#' @param rdata_file  Path to workspace \code{.RData} from \code{run_canope}.
#' @param output_dir  Directory for PDF output.
#' @param modechrom   \code{"A"}, \code{"XX"}, or \code{"XY"} chromosome filter.
#' @param prefix      Filename prefix (default: date stamp).
#' @param flank       Number of flanking targets to show on each side.
#' @param debug       Print extra messages.
#'
#' @importFrom ggplot2 ggplot geom_line geom_point geom_tile geom_ribbon geom_hline aes scale_color_manual theme theme_bw labs scale_x_continuous element_text element_blank
#' @importFrom grDevices pdf dev.off dev.list
#' @importFrom grid grid.newpage grid.layout pushViewport viewport
#' @importFrom stats ave qnbinom median sd
#' @export
generate_plots <- function(
    rdata_file,
    output_dir = "./plots",
    modechrom = "A",
    prefix = NULL,
    flank = 5L,
    debug = FALSE
) {
  message("[INFO] BEGIN plot generation")

  if (!file.exists(rdata_file)) stop("[ERROR] Summary RData not found: ", rdata_file)
  load(rdata_file)

  if (!exists("bed_file")) {
    message("[INFO] 'bed_file' not found in RData; reconstructing from 'counts'.")
    if (exists("counts") && all(c("chromosome", "start", "end") %in% colnames(counts))) {
      bed_file <- counts[, intersect(c("chromosome", "start", "end", "target", "gc", "GENE"),
                                     colnames(counts)), drop = FALSE]
    } else {
      stop("[ERROR] Cannot reconstruct 'bed_file' from counts.")
    }
  }

  required_vars <- c("cnv_calls", "counts", "bed_file", "models", "refs")
  missing_vars <- setdiff(required_vars, ls())
  if (length(missing_vars) > 0)
    stop("[ERROR] Missing variables in RData: ", paste(missing_vars, collapse = ", "))

  model_fields_missing <- vapply(models, function(m)
    !all(c("target", "test_counts", "ref_matrix") %in% names(m)), logical(1))
  if (length(models) > 0 && any(model_fields_missing))
    stop("[ERROR] 'models' in this RData is missing 'test_counts'/'ref_matrix' ",
         "(written by an older run_canope()/call_cnvs()). Re-run run_canope() ",
         "to regenerate rdata_output before plotting.")

  if (!"target" %in% colnames(bed_file) && "target" %in% colnames(counts))
    bed_file$target <- counts$target
  if (!"target" %in% colnames(bed_file))
    bed_file$target <- seq_len(nrow(bed_file))
  if (!"GENE" %in% colnames(bed_file))
    bed_file$GENE <- "TargetRegion"

  bed_file$chromosome <- format_chr_label(bed_file$chromosome)
  if ("chromosome" %in% colnames(counts))
    counts$chromosome <- format_chr_label(counts$chromosome)
  if (nrow(cnv_calls) > 0)
    cnv_calls$CHR <- format_chr_label(cnv_calls$CHR)

  prefix_str <- if (is.null(prefix) || prefix == "") format(Sys.time(), "%Y%m%d") else prefix
  if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  inc <- if (modechrom %in% c("XX", "XY")) "chrX" else NULL
  exc <- if (modechrom == "A") c("chrX", "chrY") else NULL
  cnv_plot <- canope_filter_chromosomes(cnv_calls, include = inc, exclude = exc)

  message(sprintf("[INFO] Number of CNV calls to plot: %d", nrow(cnv_plot)))

  exon_index <- if (ncol(bed_file) >= 5 && identical(colnames(bed_file)[5], "exon")) {
    bed_file[[5]]
  } else {
    compute_exon_index(bed_file)
  }

  if (nrow(cnv_plot) == 0) {
    message("[INFO] No CNV calls match chromosome filter.")
    return(invisible(NULL))
  }

  plot_count <- 0L
  for (call_idx in seq_len(nrow(cnv_plot))) {
    sample_name <- as.character(cnv_plot$SAMPLE[call_idx])
    target_chr  <- as.character(cnv_plot$CHR[call_idx])
    target_bounds <- parse_canope_targets(cnv_plot$TARGETS[call_idx])
    if (length(target_bounds) < 2L) {
      message("[WARN] Malformed TARGETS for call ", call_idx)
      next
    }
    idx_start <- target_bounds[1L]
    idx_end   <- target_bounds[2L]

    row_start <- targets_to_rows(idx_start, bed_file$target)
    row_end   <- targets_to_rows(idx_end, bed_file$target)
    if (is.na(row_start) || is.na(row_end)) {
      message("[WARN] Targets ", idx_start, "..", idx_end, " not in bed for call ", call_idx)
      next
    }

    gene_raw <- if ("GENE" %in% colnames(cnv_plot)) as.character(cnv_plot$GENE[call_idx]) else "Region"
    gene_str <- sanitize_filename(gene_raw)
    exon_range <- seq(max(1L, row_start - flank), min(nrow(bed_file), row_end + flank))

    single_chr <- length(unique(bed_file$chromosome[exon_range])) == 1L
    prev <- FALSE
    new_chr <- ""

    if (!single_chr) {
      prev    <- bed_file$chromosome[exon_range[1L]] != target_chr
      new_chr <- if (prev) bed_file$chromosome[exon_range[1L]]
      else bed_file$chromosome[utils::tail(exon_range, 1L)]
      exon_range <- exon_range[bed_file$chromosome[exon_range] == target_chr]
    }

    if (length(exon_range) == 0L) {
      message(sprintf("[WARN] No targets in range for call %d (%s)", call_idx, sample_name))
      next
    }

    cnv_rows <- targets_to_rows(seq(idx_start, idx_end), bed_file$target)
    cnv_rows <- cnv_rows[!is.na(cnv_rows)]

    tryCatch({
      if (!sample_name %in% names(refs))
        stop(sprintf("Sample '%s' not found in 'refs'", sample_name))
      if (!sample_name %in% names(models))
        stop(sprintf("Sample '%s' not found in 'models'", sample_name))
      if (!sample_name %in% colnames(counts))
        stop(sprintf("Sample '%s' not found in counts", sample_name))

      ref_samples <- refs[[sample_name]]
      if (length(ref_samples) == 0L)
        stop(sprintf("No reference samples for '%s'", sample_name))

      model_mean <- model_lookup(models[[sample_name]], "mean", bed_file$target[exon_range])
      model_var  <- model_lookup(models[[sample_name]], "var_estimate", bed_file$target[exon_range])

      pseudocount <- 0.5
      test_log <- log2(counts[exon_range, sample_name] + pseudocount)

      cov_list <- list()
      for (r in ref_samples) {
        if (!r %in% colnames(counts)) {
          warning(sprintf("Reference '%s' missing from counts, skipping", r))
          next
        }
        scaling <- stats::median(
          counts[exon_range, sample_name] / pmax(counts[exon_range, r], 1),
          na.rm = TRUE
        )
        r_log <- log2((counts[exon_range, r] * scaling) + pseudocount)
        cov_list[[length(cov_list) + 1L]] <- data.frame(
          exon_idx = exon_range, coverage = r_log, group = r,
          color_group = "Reference Sample", stringsAsFactors = FALSE
        )
      }
      if (length(cov_list) == 0L) stop("No valid reference samples")

      cov_list[[length(cov_list) + 1L]] <- data.frame(
        exon_idx = exon_range, coverage = test_log, group = "Test Sample",
        color_group = "Test Sample", stringsAsFactors = FALSE
      )
      cov_data <- do.call(rbind, cov_list)

      pt_data <- data.frame(
        exon_idx = exon_range, coverage = test_log,
        color_group = ifelse(exon_range %in% cnv_rows, "Affected exon(s)", "Test Sample"),
        stringsAsFactors = FALSE
      )
      pt_data <- pt_data[pt_data$color_group == "Affected exon(s)", , drop = FALSE]

      p_cov   <- create_coverage_plot(cov_data, pt_data, single_chr, prev, exon_range, exon_index)
      p_genes <- create_gene_tile_plot(bed_file, exon_range, single_chr, prev, new_chr)

      test_counts   <- model_lookup(models[[sample_name]], "test_counts", bed_file$target[exon_range])
      target_means  <- model_mean
      target_vars   <- model_var
      expected_safe <- pmax(target_means, 1) + pseudocount

      ratio <- log2((test_counts + pseudocount) / expected_safe)
      dispersion <- target_vars - target_means
      size_param <- ifelse(dispersion > 0, target_means^2 / dispersion, Inf)

      alpha <- 0.05
      mins <- stats::qnbinom(alpha / 2, mu = target_means, size = size_param)
      maxs <- stats::qnbinom(1 - alpha / 2, mu = target_means, size = size_param)
      lo <- log2((mins + pseudocount) / expected_safe)
      hi <- log2((maxs + pseudocount) / expected_safe)

      ci_data <- data.frame(
        exon = exon_range, ratio = ratio, lo = lo, hi = hi,
        is_affected = factor(
          exon_range %in% cnv_rows,
          levels = c(FALSE, TRUE),
          labels = c("Observed", "Affected")
        )
      )


      bg_calib <- check_background_calibration(
        ci_data$ratio, ci_data$lo, ci_data$hi, ci_data$is_affected == "Affected"
      )
      ci_subtitle <- "95% NB predictive interval"
      if (isTRUE(bg_calib$flag)) {
        ci_subtitle <- sprintf(
          "95%% NB predictive interval  [%d%% of background exons outside CI (%d/%d) \u2014 check region/reference match]",
          round(bg_calib$pct_outside), bg_calib$n_outside, bg_calib$n_background
        )
      }

      p_ci <- ggplot2::ggplot(ci_data, ggplot2::aes(x = exon)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey90", colour = NA) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
        ggplot2::geom_line(ggplot2::aes(y = ratio), colour = "steelblue", linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes(y = ratio, color = is_affected), size = 3) +
        ggplot2::scale_color_manual(values = c("Observed" = "blue", "Affected" = "red")) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = "", y = "log2(Observed / Expected)",
                      subtitle = ci_subtitle)
      p_ci <- apply_xaxis_formatting(p_ci, single_chr, prev, exon_range, exon_index)

      p_z <- create_zscore_plot(
        models[[sample_name]], exon_range, bed_file$target[exon_range], sample_name, ref_samples,
        cnv_rows, single_chr, prev, exon_index
      )

      interval_label <- as.character(cnv_plot$INTERVAL[call_idx])
      file_path <- file.path(
        output_dir,
        paste0(prefix_str, ".", sample_name, ".", target_chr, "_",
               gene_str, "_", call_idx, ".pdf")
      )

      grDevices::pdf(file_path, useDingbats = FALSE, width = 8, height = 12)
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(8, 1)))
      print(p_cov,   vp = grid::viewport(layout.pos.row = 1:3, layout.pos.col = 1))
      print(p_genes, vp = grid::viewport(layout.pos.row = 4,   layout.pos.col = 1))
      print(p_ci,    vp = grid::viewport(layout.pos.row = 5:6, layout.pos.col = 1))
      print(p_z,     vp = grid::viewport(layout.pos.row = 7:8, layout.pos.col = 1))
      grDevices::dev.off()

      plot_count <- plot_count + 1L
      if (debug)
        message(sprintf("[INFO] Saved %s (%s)", file_path, interval_label))
      else
        message(sprintf("[INFO] Plot %d/%d saved: %s", plot_count, nrow(cnv_plot), basename(file_path)))

    }, error = function(e) {
      if (length(grDevices::dev.list()) > 0L) grDevices::dev.off()
      message("[ERROR] Call ", call_idx, " (", sample_name, "): ", conditionMessage(e))
    })
  }

  message(sprintf("[INFO] END plot generation – %d/%d plots created",
                  plot_count, nrow(cnv_plot)))
  invisible(plot_count)
}


#' Look Up Per-Target Model Statistics
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
  if (any(is.na(vals)))
    warning("Some targets missing model statistics; using median imputation")
  vals[is.na(vals)] <- stats::median(model[[field]], na.rm = TRUE)
  vals
}


#' Look Up Per-Target Model Statistics (matrix form)
#'
#' @noRd
model_matrix_lookup <- function(model, field, target_ids) {
  mat <- model[[field]]
  idx <- match(as.integer(target_ids), as.integer(model$target))
  out <- mat[idx, , drop = FALSE]
  if (any(is.na(idx)))
    warning("Some targets missing model statistics; using median imputation")
  for (j in seq_len(ncol(out))) {
    miss <- is.na(out[, j])
    if (any(miss)) out[miss, j] <- stats::median(mat[, j], na.rm = TRUE)
  }
  out
}


#' Z-Score Panel vs Reference Samples
#'
#' @param model      \code{models[[sample_name]]} — must carry
#'   \code{target}, \code{test_counts}, \code{var_estimate}, and
#'   \code{ref_matrix} (columns named by reference sample, in the same
#'   order as \code{ref_samples}).
#' @param exon_range Row positions (into \code{bed_file}) of the exon window
#'   being plotted — used only for x-axis placement.
#' @param target_ids \code{bed_file$target[exon_range]} — the target IDs
#'   used to align model data, since model rows are not necessarily in the
#'   same order/subset as \code{bed_file}.
#' @noRd
create_zscore_plot <- function(model, exon_range, target_ids, sample_name, ref_samples,
                               cnv_rows, single_chr, prev, exon_index) {
  test_vals <- model_lookup(model, "test_counts", target_ids)
  ref_mat   <- model_matrix_lookup(model, "ref_matrix", target_ids)
  model_var <- model_lookup(model, "var_estimate", target_ids)

  means <- rowMeans(ref_mat, na.rm = TRUE)
  sds   <- sqrt(pmax(model_var, 1))
  sds[!is.finite(sds) | sds < 1e-6] <- 1e-6

  sample_z <- (test_vals - means) / sds
  ref_z <- sapply(seq_along(ref_samples), function(i) (ref_mat[, i] - means) / sds)

  z_lim <- suppressWarnings(max(abs(c(ref_z, sample_z)), na.rm = TRUE) * 1.1)
  if (!is.finite(z_lim)) z_lim <- 1
  z_lim <- max(z_lim, 1)

  z_df <- data.frame(
    exon = rep(exon_range, times = length(ref_samples) + 1L),
    z = c(as.vector(ref_z), sample_z),
    series = rep(c(ref_samples, sample_name), each = length(exon_range)),
    highlight = rep(exon_range %in% cnv_rows, times = length(ref_samples) + 1L)
  )
  z_df$series <- factor(z_df$series, levels = c(ref_samples, sample_name))

  p <- ggplot2::ggplot(z_df, ggplot2::aes(x = exon, y = z, group = series)) +
    ggplot2::geom_line(
      data = subset(z_df, series != sample_name),
      colour = "gray60", linewidth = 0.6, alpha = 0.8
    ) +
    ggplot2::geom_line(
      data = subset(z_df, series == sample_name),
      colour = "red", linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = subset(z_df, series == sample_name & highlight),
      colour = "darkred", size = 2.5
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::coord_cartesian(ylim = c(-z_lim, z_lim)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "Z-score vs references",
                  subtitle = "Test sample (red) with CNV targets highlighted")

  apply_xaxis_formatting(p, single_chr, prev, exon_range, exon_index)
}


#' @noRd
apply_xaxis_formatting <- function(p, single_chr, prev, exon_range, exon_index) {
  if (length(exon_range) == 0L) return(p)
  min_e <- min(exon_range)
  max_e <- max(exon_range)

  if (single_chr) {
    b <- exon_range
    l <- exon_index[exon_range]
    lim <- NULL
  } else if (prev) {
    b <- (min_e - 6L):max_e
    l <- c(rep("", 6L), exon_index[exon_range])
    lim <- c(min_e - 6.75, max_e)
  } else {
    b <- min_e:(max_e + 6L)
    l <- c(exon_index[exon_range], rep("", 6L))
    lim <- c(min_e, max_e + 6.75)
  }

  p + ggplot2::scale_x_continuous(breaks = b, labels = l, limits = lim) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, size = 6, hjust = 1))
}


#' @noRd
create_coverage_plot <- function(cov_data, pt_data, single_chr, prev, exon_range, exon_index) {
  cols <- c("Reference Sample" = "gray", "Test Sample" = "blue", "Affected exon(s)" = "red")
  p_cov <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = subset(cov_data, color_group == "Reference Sample"),
      ggplot2::aes(x = exon_idx, y = coverage, group = group, color = color_group),
      linetype = "dashed", linewidth = 1.1
    ) +
    ggplot2::geom_point(
      data = subset(cov_data, color_group == "Reference Sample"),
      ggplot2::aes(x = exon_idx, y = coverage, color = color_group), size = 2
    ) +
    ggplot2::geom_line(
      data = subset(cov_data, color_group == "Test Sample"),
      ggplot2::aes(x = exon_idx, y = coverage, group = group, color = color_group),
      linetype = "dashed", linewidth = 1.1
    ) +
    ggplot2::geom_point(
      data = subset(cov_data, color_group == "Test Sample"),
      ggplot2::aes(x = exon_idx, y = coverage, color = color_group), size = 2
    ) +
    ggplot2::geom_point(
      data = pt_data,
      ggplot2::aes(x = exon_idx, y = coverage, color = color_group), size = 3
    ) +
    ggplot2::scale_colour_manual(values = cols, guide = "legend") +
    ggplot2::labs(y = "log2(Coverage + 0.5)", x = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())
  apply_xaxis_formatting(p_cov, single_chr, prev, exon_range, exon_index)
}


#' @noRd
create_gene_tile_plot <- function(bed_file, exon_range, single_chr, prev, new_chr) {
  if (length(exon_range) == 0L) return(ggplot2::ggplot())
  temp <- cbind(row = seq_len(nrow(bed_file)), bed_file)[exon_range, , drop = FALSE]

  gene_names <- unique(bed_file$GENE[exon_range])
  gene_tiles <- data.frame(
    gene  = gene_names,
    mid   = vapply(gene_names, function(g) mean(exon_range[temp$GENE == g]), numeric(1)),
    width = vapply(gene_names, function(g) sum(temp$GENE == g), numeric(1)) - 0.5,
    y     = 1,
    stringsAsFactors = FALSE
  )
  if (!single_chr) {
    gene_tiles <- rbind(
      gene_tiles,
      data.frame(
        gene = new_chr,
        mid  = if (prev) min(exon_range) - 5 else max(exon_range) + 5,
        width = 3.5, y = 1, stringsAsFactors = FALSE
      )
    )
  }
  ggplot2::ggplot(gene_tiles, ggplot2::aes(x = mid, y = y, fill = gene, width = width, label = gene)) +
    ggplot2::geom_tile(alpha = 0.7) +
    ggplot2::geom_text(size = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "", y = "")
}
