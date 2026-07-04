#' Stop if a value is NULL or empty
#' @param val Object to validate.
#' @param msg Error message.
#' @export
stop_if_missing <- function(val, msg) {
  if (is.null(val) || length(val) == 0) stop(msg)
}

#' Stop if a file does not exist
#' @param path File path.
#' @param msg Error message.
#' @export
stop_if_not_file <- function(path, msg) {
  if (is.null(path) || !file.exists(path)) stop(msg)
}

#' Load an RData file into an isolated environment
#' @param path Path to the .RData file.
#' @param required Optional character vector of required object names.
#' @return Named list of objects loaded from the file.
#' @export
canope_load_rdata <- function(path, required = NULL) {
  stop_if_not_file(path, paste0("[ERROR] RData not found: ", path))
  env <- new.env()
  load(path, envir = env)
  if (!is.null(required)) {
    missing_vars <- setdiff(required, ls(env))
    if (length(missing_vars) > 0)
      stop("[ERROR] Missing variables in RData: ", paste(missing_vars, collapse = ", "),
           "\nAvailable: ", paste(ls(env), collapse = ", "))
  }
  as.list(env)
}

#' Null-coalescing operator
#' @param a Primary value.
#' @param b Fallback value.
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Sanitize a string for use as a filename
#' @param name Character string.
#' @export
sanitize_filename <- function(name) gsub("[^[:alnum:]._-]", "_", as.character(name))

#' Filter a data frame by a chromosome column
#'
#' Handles both \code{CHR} and \code{chromosome} column names and normalises
#' chr-prefix mismatches automatically.
#'
#' @param df Data frame with a \code{CHR} or \code{chromosome} column.
#' @param include Optional chromosome names to keep.
#' @param exclude Optional chromosome names to remove.
#' @return Filtered data frame.
#' @export
canope_filter_chromosomes <- function(df, include = NULL, exclude = NULL) {
  if (is.null(df) || nrow(df) == 0) return(df)
  chrom_col <- intersect(c("CHR", "chromosome", "Chromosome"), colnames(df))[1]
  if (is.na(chrom_col)) stop("Data frame must have a CHR / chromosome column.")
  norm <- function(x) unique(c(x, sub("^chr", "", x), paste0("chr", sub("^chr", "", x))))
  if (!is.null(include)) df <- df[df[[chrom_col]] %in% norm(include), , drop = FALSE]
  if (!is.null(exclude)) df <- df[!df[[chrom_col]] %in% norm(exclude), , drop = FALSE]
  df
}

#' Assign sequential exon numbers within each gene
#'
#' Sorts by chromosome → start → end then assigns 1..n per gene.
#' Duplicate intervals (same chrom + start + end + gene) are removed.
#'
#' @param bed_df Data frame with chromosome/Chr, start/Start, end/End, gene/Gene/GENE columns.
#' @return The input data frame with an added \code{exon_number} integer column.
#' @importFrom data.table := as.data.table setnames setorder .N
#' @export
assign_exon_numbers_per_gene <- function(bed_df) {
  # Detect column aliases
  chrom_col <- intersect(c("chromosome", "Chr", "CHROM"), names(bed_df))[1]
  start_col <- intersect(c("start", "Start", "START"), names(bed_df))[1]
  end_col   <- intersect(c("end", "End", "END"), names(bed_df))[1]
  gene_col  <- intersect(c("GENE", "gene", "Gene"), names(bed_df))[1]
  stopifnot(!is.na(chrom_col), !is.na(start_col), !is.na(end_col), !is.na(gene_col))

  dt <- data.table::as.data.table(bed_df)
  data.table::setnames(dt, c(chrom_col, start_col, end_col, gene_col),
                       c("._chrom", "._start", "._end", "._gene"))

  # Remove duplicates
  dup_rows <- duplicated(dt, by = c("._chrom", "._start", "._end", "._gene"))
  if (any(dup_rows)) {
    warning(sprintf("Removed %d duplicate BED rows (same chrom, start, end, gene).",
                    sum(dup_rows)), immediate. = TRUE)
    dt <- dt[!dup_rows]
  }

  # Sort genomically
  chrom_levels <- c(paste0("chr", c(1:22, "X", "Y", "M")), c(as.character(1:22), "X", "Y", "M"))
  dt[, .chrom_fac := factor(`._chrom`, levels = unique(c(chrom_levels, unique(`._chrom`))))]
  data.table::setorder(dt, .chrom_fac, `._start`, `._end`)
  dt[, .chrom_fac := NULL]

  # Exon numbering within each gene (genomic order)
  dt[, exon_number := seq_len(.N), by = "._gene"]

  data.table::setnames(dt, c("._chrom", "._start", "._end", "._gene"),
                       c(chrom_col, start_col, end_col, gene_col))
  as.data.frame(dt)
}

#' Compute within-gene exon index from a BED-like data frame
#'
#' Returns the exon_number column if already present, otherwise calls
#' \code{\link{assign_exon_numbers_per_gene}}.
#'
#' @param bed_df Data frame (as above).
#' @return Integer vector, one per row.
#' @export
compute_exon_index <- function(bed_df) {
  if ("exon_number" %in% names(bed_df) && !all(is.na(bed_df$exon_number)))
    return(bed_df$exon_number)
  assign_exon_numbers_per_gene(bed_df)$exon_number
}

#' Normalise chromosome names to match a reference naming style (vectorised)
#'
#' @param chr_vec Character vector of chromosome names to normalise.
#' @param ref_chromosomes Reference vector that defines which style to adopt.
#' @return Character vector with chr-prefix aligned to the reference style.
#' @export
normalize_chromosome_vec <- function(chr_vec, ref_chromosomes) {
  has_chr <- any(grepl("^chr", ref_chromosomes))
  if (has_chr) ifelse(grepl("^chr", chr_vec), chr_vec, paste0("chr", chr_vec))
  else          ifelse(grepl("^chr", chr_vec), sub("^chr", "", chr_vec), chr_vec)
}

#' Check whether a BAM file has an index on disk
#' @param bam Path to BAM file.
#' @return Logical.
#' @noRd
bam_has_index <- function(bam) {
  file.exists(paste0(bam, ".bai")) || file.exists(paste0(bam, ".bam.bai"))
}

#' Parse a CANOPE INTERVAL string into components
#'
#' @param interval Character string of the form \code{"chr1:100-200"}.
#' @return Named list with \code{chrom}, \code{start}, \code{end}.
#' @export
parse_canope_interval <- function(interval) {
  interval <- as.character(interval)
  colon    <- regexpr(":", interval, fixed = TRUE)
  chrom    <- substr(interval, 1L, colon - 1L)
  rest     <- substr(interval, colon + 1L, nchar(interval))
  dash     <- regexpr("-", rest, fixed = TRUE)
  list(
    chrom = chrom,
    start = as.integer(substr(rest, 1L, dash - 1L)),
    end   = as.integer(substr(rest, dash + 1L, nchar(rest)))
  )
}

#' Parse a CANOPE TARGETS string into target ID range
#'
#' @param targets Character string of the form \code{"42..47"}.
#' @return Integer vector \code{c(start_target, end_target)}.
#' @export
parse_canope_targets <- function(targets) {
  parts <- as.integer(strsplit(as.character(targets), "..", fixed = TRUE)[[1]])
  parts[!is.na(parts)]
}

# ─── PCA of coverage profiles ──────────────────────────────────────────────

#' Plot PCA of Sample Coverage Profiles (CANOPE version)
#'
#' Mirrors ECHO's \code{plot_coverage_pca} but uses CANOPE's column layout
#' where sample columns begin after the metadata columns.
#'
#' @param counts      Data frame (rows = targets, columns = metadata + samples).
#' @param sample_names Character vector of sample column names.
#' @param output_pdf  Optional path for PDF output (NULL for on-screen).
#' @param color_by    Optional factor/character vector of group labels per sample.
#' @param scale       Logical; scale variables before PCA. Default TRUE.
#' @return Invisibly returns list(pca, var_exp).
#' @export
plot_coverage_pca <- function(counts, sample_names,
                              output_pdf = NULL, color_by = NULL, scale = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for PCA plot.")

  count_mat <- as.matrix(counts[, sample_names, drop = FALSE])
  log_mat   <- log2(count_mat + 1)
  row_var   <- apply(log_mat, 1, var, na.rm = TRUE)
  keep      <- row_var > 0 & !is.na(row_var)
  if (sum(keep) < 2) {
    warning("Insufficient variation for PCA (< 2 informative rows).")
    return(invisible(NULL))
  }
  log_mat <- log_mat[keep, , drop = FALSE]

  pca     <- stats::prcomp(t(log_mat), scale. = scale, center = TRUE)
  var_exp <- summary(pca)$importance[2, ] * 100
  pca_df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         Sample = sample_names, stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = Sample))

  if (!is.null(color_by) && length(color_by) == nrow(pca_df)) {
    pca_df$Group <- as.factor(color_by)
    p <- p + ggplot2::aes(colour = Group)
    if (length(unique(color_by)) >= 2)
      p <- p + ggplot2::stat_ellipse(ggplot2::aes(colour = Group),
                                      type = "norm", linetype = "dashed")
  }

  p <- p +
    ggplot2::geom_point(size = 3,
                        colour = if (is.null(color_by)) "steelblue" else NULL) +
    ggplot2::labs(
      x     = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y     = paste0("PC2 (", round(var_exp[2], 1), "%)"),
      title = "PCA – CANOPE Coverage Profiles"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.title  = ggplot2::element_text(face = "bold")
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(size = 3, max.overlaps = 15)
  } else {
    p <- p + ggplot2::geom_text(size = 3, hjust = -0.2, check_overlap = TRUE)
  }

  if (!is.null(output_pdf)) {
    ggplot2::ggsave(output_pdf, p, width = 8, height = 6)
    message("[INFO] PCA plot saved: ", output_pdf)
  } else {
    print(p)
  }
  invisible(list(pca = pca, var_exp = var_exp))
}


#' Flag Background-Exon Calibration Issues for a CNV Call Window
#'
#' New in response to a real investigation (see README "Round 5"): checking
#' a real report's data directly showed the test sample's own ratio at
#' *non-called* ("background") exons in a plotted window falling outside
#' the modelled 95% predictive interval far more often than the ~5% a
#' well-calibrated interval implies — but concentrated entirely in specific
#' calls, in a pattern (near-uniform, one-sided elevation across most of the
#' window) much more consistent with real signal extending beyond the
#' called boundary, an atypical reference match, or a technical/batch
#' difference for that sample than with a generically miscalibrated
#' interval. A blanket statistical correction was tested (Monte Carlo, see
#' README) and rejected — it overcorrected the normal case. This makes that
#' same per-call check a permanent, automatic part of the pipeline instead,
#' so it doesn't require manually parsing plot data to notice: every call
#' gets flagged (or not) using a real statistical test against the null.
#'
#' This is a diagnostic flag, not a correction — it doesn't change the
#' interval, the call, or the confidence score. It's meant to prompt a
#' manual look at specific calls, the same way this issue was actually
#' found.
#'
#' @param ratio Numeric vector of log2(observed/expected) for every exon in
#'   the plotted window (background and affected together).
#' @param lo,hi Numeric vectors (same length as \code{ratio}) giving the
#'   95% predictive interval bounds at each exon.
#' @param is_affected Logical vector (same length); \code{TRUE} for exons
#'   already called as part of this CNV — excluded from the check, since
#'   those are expected to sit outside the interval.
#' @param min_n Minimum number of background exons required before
#'   flagging (default 5) — below this the percentage is too noisy on its
#'   own to test meaningfully.
#' @return A list with \code{n_background}, \code{n_outside},
#'   \code{pct_outside}, and \code{flag} — \code{flag} is \code{TRUE} when a
#'   one-sided binomial test of \code{n_outside} against a 5% null rate is
#'   significant at p < 0.05.
#' @export
check_background_calibration <- function(ratio, lo, hi, is_affected, min_n = 5) {
  bg <- !is_affected
  n_bg <- sum(bg, na.rm = TRUE)
  if (n_bg == 0) {
    return(list(n_background = 0L, n_outside = 0L, pct_outside = NA_real_, flag = FALSE))
  }
  outside <- (ratio[bg] < lo[bg]) | (ratio[bg] > hi[bg])
  outside[is.na(outside)] <- FALSE
  n_outside <- sum(outside)
  pct_outside <- 100 * n_outside / n_bg
  flag <- n_bg >= min_n &&
    stats::pbinom(n_outside - 1L, size = n_bg, prob = 0.05, lower.tail = FALSE) < 0.05
  list(n_background = n_bg, n_outside = n_outside, pct_outside = pct_outside, flag = flag)
}
