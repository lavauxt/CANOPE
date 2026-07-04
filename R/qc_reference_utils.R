#' Automatic Reference Sample Inference
#'
#' Identifies clean (low-noise) samples and selects top correlated
#' references for each, excluding CNV-heavy samples.
#'
#' @param counts                Numeric matrix (targets × samples).
#' @param top_n                 Integer. Max references per sample.
#' @param cor_method            Character. Correlation method ("spearman").
#' @param contamination_quantile Numeric. MAD quantile above which a sample
#'                              is considered noisy.
#' @param pseudocount           Numeric. Added before log2 transform.
#'
#' @return Named list: refs, clean_samples, anomaly_score, normalized_counts.
#' @export
infer_reference_samples <- function(
    counts,
    top_n                    = 15L,
    cor_method               = "spearman",
    contamination_quantile   = 0.90,
    pseudocount              = 0.5
) {
  counts <- as.matrix(counts)
  if (ncol(counts) < 2L) {
    stop("infer_reference_samples requires at least 2 samples; got ", ncol(counts), ".")
  }

  log_counts     <- log2(counts + pseudocount)
  sample_medians <- apply(log_counts, 2, median, na.rm = TRUE)
  norm_counts    <- sweep(log_counts, 2, sample_medians, "-")

  mad_scores <- apply(norm_counts, 2, function(x) median(abs(x), na.rm = TRUE))
  threshold  <- stats::quantile(mad_scores, contamination_quantile, na.rm = TRUE)

  valid_idx     <- which(mad_scores <= threshold & is.finite(mad_scores))
  clean_samples <- if (length(valid_idx) >= 3)
    colnames(norm_counts)[valid_idx] else colnames(norm_counts)

  cor_mat <- stats::cor(norm_counts[, clean_samples, drop = FALSE],
                        method = cor_method,
                        use    = "pairwise.complete.obs")

  refs <- lapply(clean_samples, function(s) {
    cors      <- sort(cor_mat[, s], decreasing = TRUE)
    cors      <- cors[!is.na(cors)]
    valid_refs <- names(cors)[names(cors) != s]
    head(valid_refs, top_n)
  })
  names(refs) <- clean_samples

  list(
    refs             = refs,
    clean_samples    = clean_samples,
    anomaly_score    = mad_scores,
    normalized_counts = norm_counts
  )
}


#' Detect Outlier Samples  (robust z-score on MAD noise)
#'
#' @param counts      Numeric matrix or data frame.
#' @param pseudocount Numeric.
#' @param z_threshold Numeric. Robust z-score above which a sample is an outlier.
#'
#' @return Data frame: sample, noise_score, robust_z, is_outlier.
#' @export
detect_outlier_samples <- function(counts, pseudocount = 0.5, z_threshold = 3) {
  counts <- as.matrix(counts)
  log_counts   <- log2(counts + pseudocount)
  sample_noise <- apply(log_counts, 2, stats::mad, na.rm = TRUE)

  med  <- median(sample_noise, na.rm = TRUE)
  mad0 <- stats::mad(sample_noise, na.rm = TRUE)

  # Guard against mad0 == 0 (all samples identical noise) or a single sample
  z_scores  <- if (length(sample_noise) > 1L && is.finite(mad0) && mad0 > 0)
    (sample_noise - med) / mad0 else rep(0, length(sample_noise))

  is_outlier <- abs(z_scores) > z_threshold

  data.frame(
    sample       = names(sample_noise),
    noise_score  = sample_noise,
    robust_z     = z_scores,
    is_outlier   = is_outlier,
    stringsAsFactors = FALSE
  )
}


#' Detect Problematic Exons
#'
#' Flags exons with high MAD (cross-sample noise), low mean coverage,
#' extreme GC content, or non-finite values.  Ensures at least one
#' exon per chromosome is retained to prevent empty chromosomes.
#'
#' @param count_matrix Numeric matrix (targets × samples).
#' @param chromosomes  Optional factor/character vector of chromosome labels.
#' @param mad_quantile Numeric. Top quantile of exon MAD to flag.
#' @param min_mean     Numeric. Minimum mean coverage to retain.
#' @param gc           Optional numeric vector of GC content (0–1 or 0–100).
#' @param gc_min       Numeric. Minimum GC fraction.
#' @param gc_max       Numeric. Maximum GC fraction.
#'
#' @return Data frame: exon (row index), mean, mad, problematic.
#' @export
detect_problematic_exons <- function(
    count_matrix,
    chromosomes  = NULL,
    mad_quantile = 0.90,
    min_mean     = 20,
    gc           = NULL,
    gc_min       = 0.10,
    gc_max       = 0.90
) {
  count_matrix <- as.matrix(count_matrix)
  exon_mean    <- rowMeans(count_matrix, na.rm = TRUE)
  exon_mad     <- apply(count_matrix, 1, stats::mad, na.rm = TRUE)
  mad_thresh   <- stats::quantile(exon_mad, probs = mad_quantile, na.rm = TRUE)

  problematic  <- (exon_mad > mad_thresh | exon_mean < min_mean |
                     !is.finite(exon_mean))

  if (!is.null(gc)) {
    if (length(gc) != nrow(count_matrix))
      stop("'gc' must have one value per row of count_matrix (", nrow(count_matrix),
           " rows, got ", length(gc), ").")
    gc_val       <- if (max(gc, na.rm = TRUE) > 1) gc / 100 else gc
    problematic  <- problematic | gc_val < gc_min | gc_val > gc_max |
      !is.finite(gc_val)
  }

  # Guarantee at least one exon per chromosome
  if (!is.null(chromosomes)) {
    for (chr in unique(chromosomes)) {
      idx <- which(chromosomes == chr)
      if (length(idx) > 0 && all(problematic[idx]))
        problematic[idx[which.min(exon_mad[idx])]] <- FALSE
    }
  }

  data.frame(
    exon        = seq_len(nrow(count_matrix)),
    mean        = exon_mean,
    mad         = exon_mad,
    problematic = problematic
  )
}


#' GC-Content Correction via LOESS
#'
#' Fits a LOESS regression of log2(count) on GC fraction per sample and
#' subtracts the fitted trend, preserving the median level.
#'
#' @param counts Numeric matrix (targets × samples).
#' @param gc     Numeric vector of GC fractions (0–1 or 0–100).
#' @param method Character. Only "loess" is currently supported.
#'
#' @return Integer matrix of corrected counts (non-negative).
#' @export
gc_correct_counts <- function(counts, gc, method = "loess") {
  counts  <- as.matrix(counts)
  if (length(gc) != nrow(counts))
    stop("'gc' must have one value per row of counts (", nrow(counts),
         " rows, got ", length(gc), ").")
  gc_val  <- if (max(gc, na.rm = TRUE) > 1) gc / 100 else gc
  corrected <- counts

  if (method == "loess") {
    for (i in seq_len(ncol(counts))) {
      y      <- counts[, i]
      valid  <- is.finite(y) & is.finite(gc_val) & y > 0
      if (sum(valid) >= 50) {
        fit <- tryCatch(
          stats::loess(log2(y[valid] + 1) ~ gc_val[valid],
                       span = 0.4, degree = 2),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          pred <- stats::predict(fit, gc_val)
          med_pred <- median(pred, na.rm = TRUE)
          pred[!is.finite(pred)] <- med_pred
          corrected[, i] <- round(2^(log2(y + 1) - pred + med_pred) - 1)
        }
      }
    }
  }

  corrected[corrected < 0] <- 0L
  storage.mode(corrected) <- "integer"
  return(corrected)
}
