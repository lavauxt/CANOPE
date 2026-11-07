#' Compute CANOPE QC Metrics
#'
#' Evaluates coverage quality across all samples and exons.  Writes flagged rows
#' to a TSV and returns them invisibly.
#'
#' @param counts       Data frame (rows = targets, metadata + sample columns).
#' @param sample_names Character vector of sample column names.
#' @param bed_df       BED data frame aligned with \code{counts} rows.
#'   Must have chromosome/Chr, start/Start, end/End, GENE/gene columns.
#' @param output_file  Path for the TSV output. Default \code{"./CANOPE_QC_metrics.tsv"}.
#' @param min_corr     Minimum Spearman correlation between any two samples. Default 0.98.
#' @param min_cov      Minimum median read depth per sample / exon. Default 100.
#' @param min_total_reads Minimum total reads per sample. Default 300 000.
#' @param max_exon_cv  Maximum allowed coefficient of variation (CV) per exon. Default 0.5.
#' @param all_chromosomes Optional character vector of every chromosome present
#'   in the input BED / count table \emph{before} \code{modechrom} filtering.
#'   When supplied, a chromosome that is in the BED but was left out of the
#'   analysis (by \code{modechrom}/\code{removeY}) is reported as
#'   \code{"Chromosome not analysed"}, instead of being wrongly reported as
#'   missing from the BED. When \code{NULL} the chromosomes of \code{bed_df}
#'   are used (legacy behaviour).
#' @param modechrom,removeY The \code{run_canope()} chromosome settings, used
#'   to explain why a chromosome class was not analysed.
#' @param excluded_samples Optional data frame (\code{sample}, \code{reason}) of
#'   samples removed before calling (e.g. noisy outliers); they are listed as
#'   whole-sample flags so they do not silently vanish from the QC table.
#'
#' @return Invisibly returns the metrics data frame.
#' @export
run_canope_qc_metrics <- function(counts,
                                   sample_names,
                                   bed_df,
                                   output_file      = "./CANOPE_QC_metrics.tsv",
                                   min_corr         = 0.98,
                                   min_cov          = 100,
                                   min_total_reads  = 300000,
                                   max_exon_cv      = 0.5,
                                   all_chromosomes  = NULL,
                                   modechrom        = "A",
                                   removeY          = TRUE,
                                   excluded_samples = NULL) {
  message("[INFO] BEGIN CANOPE QC Metrics")

  gene_col  <- intersect(c("GENE", "gene", "Gene"), names(bed_df))[1]
  chrom_col <- intersect(c("chromosome", "Chr", "CHROM"), names(bed_df))[1]
  if (is.na(gene_col))  gene_col  <- names(bed_df)[4]
  if (is.na(chrom_col)) chrom_col <- names(bed_df)[1]


  gene_vals  <- bed_df[[gene_col]]
  chrom_vals <- bed_df[[chrom_col]]

  bed_numbered  <- assign_exon_numbers_per_gene(bed_df)
  exon_in_gene  <- bed_numbered$exon_number

  dt_counts     <- data.table::as.data.table(counts[, sample_names, drop = FALSE])
  sample_median <- sapply(dt_counts, stats::median, na.rm = TRUE)
  total_reads   <- colSums(dt_counts, na.rm = TRUE)

  if (length(sample_names) < 2) {
    max_corr <- setNames(rep(NA_real_, length(sample_names)), sample_names)
    message("[WARNING] Correlation QC skipped: fewer than 2 samples.")
  } else {
    corr_matrix <- stats::cor(dt_counts, method = "spearman",
                              use = "pairwise.complete.obs")
    # Exclude the diagonal by position, not by value: filtering on "!= 1" would
    # also hide a second sample that correlates perfectly (a duplicate).
    diag(corr_matrix) <- NA_real_
    max_corr <- apply(corr_matrix, 1, function(x) {
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    })
    names(max_corr) <- colnames(corr_matrix)
  }

  m <- list(Sample  = character(), Exon    = character(),
            Type    = character(), Details = character(),
            Gene    = character())

  add_metric <- function(samples, exons, type, details, genes) {
    m$Sample  <<- c(m$Sample,  samples)
    m$Exon    <<- c(m$Exon,    exons)
    m$Type    <<- c(m$Type,    rep(type, length(samples)))
    m$Details <<- c(m$Details, details)
    m$Gene    <<- c(m$Gene,    genes)
  }

  if (!is.null(excluded_samples) && nrow(excluded_samples) > 0)
    add_metric(as.character(excluded_samples$sample), rep("All", nrow(excluded_samples)),
               "Whole sample",
               paste("Excluded from analysis:", excluded_samples$reason),
               rep("All", nrow(excluded_samples)))

  low_corr <- which(!is.na(max_corr) & max_corr < min_corr)
  if (length(low_corr))
    add_metric(sample_names[low_corr], rep("All", length(low_corr)),
               "Whole sample",
               paste("Low correlation:", round(max_corr[low_corr], 3)),
               rep("All", length(low_corr)))

  low_med <- which(sample_median < min_cov)
  if (length(low_med))
    add_metric(sample_names[low_med], rep("All", length(low_med)),
               "Whole sample",
               paste("Low median depth:", round(sample_median[low_med], 1)),
               rep("All", length(low_med)))

  low_reads <- which(total_reads < min_total_reads)
  if (length(low_reads))
    add_metric(sample_names[low_reads], rep("All", length(low_reads)),
               "Whole sample",
               paste("Low total reads:", format(total_reads[low_reads], big.mark = ",")),
               rep("All", length(low_reads)))

  exon_median <- apply(dt_counts, 1, stats::median, na.rm = TRUE)
  fail_exon   <- which(exon_median < min_cov)
  if (length(fail_exon)) {
    elabels <- paste0(gene_vals[fail_exon], ":", exon_in_gene[fail_exon])
    add_metric(rep("All", length(fail_exon)), elabels,
               "Whole exon",
               paste("Low median depth:", round(exon_median[fail_exon], 1)),
               gene_vals[fail_exon])
  }

  if (requireNamespace("matrixStats", quietly = TRUE)) {
    exon_sd  <- matrixStats::rowSds(as.matrix(dt_counts), na.rm = TRUE)
    exon_mean <- rowMeans(dt_counts, na.rm = TRUE)
  } else {
    exon_sd   <- apply(dt_counts, 1, stats::sd,   na.rm = TRUE)
    exon_mean <- apply(dt_counts, 1, base::mean, na.rm = TRUE)
  }
  exon_cv  <- exon_sd / pmax(exon_mean, 1)
  high_cv  <- which(exon_cv > max_exon_cv & !is.na(exon_cv))
  if (length(high_cv)) {
    elabels <- paste0(gene_vals[high_cv], ":", exon_in_gene[high_cv])
    add_metric(rep("All", length(high_cv)), elabels,
               "Exon variability",
               paste0("High CV (>", max_exon_cv, "): ", round(exon_cv[high_cv], 3)),
               gene_vals[high_cv])
  }

  input_chr <- if (!is.null(all_chromosomes)) as.character(all_chromosomes) else as.character(chrom_vals)
  cs <- canope_chromosome_status(input_chr, as.character(chrom_vals),
                                 modechrom = modechrom, removeY = removeY)
  for (k in seq_len(nrow(cs$status))) {
    st <- cs$status[k, ]
    if (!st$in_input) {
      # Only the sex chromosomes were ever reported as missing; keep that
      # (an autosome-free BED is reported too, as it makes the run meaningless).
      message("[WARNING] ", st$label, " not found in BED/counts - ", st$label, " CNVs cannot be called.")
      add_metric("All", "All", "Missing Chromosome",
                 sprintf("%s not present in BED", st$label), "All")
    } else if (!st$analysed) {
      message("[WARNING] ", st$label, " present in the BED but not analysed: ", st$reason)
      add_metric("All", "All", "Chromosome not analysed",
                 sprintf("%s present in BED but not analysed: %s", st$label, st$reason), "All")
    }
  }

  final_metrics <- data.frame(m, stringsAsFactors = FALSE)
  if (nrow(final_metrics) == 0)
    final_metrics <- data.frame(Sample = character(), Exon = character(),
                                Type = character(), Details = character(),
                                Gene = character(), stringsAsFactors = FALSE)

  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(final_metrics, file = output_file, sep = "\t",
                     quote = FALSE, row.names = FALSE)
  message("[INFO] QC metrics written to: ", output_file)
  message("[INFO] END CANOPE QC Metrics")
  invisible(final_metrics)
}
