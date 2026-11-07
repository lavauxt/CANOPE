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

#' Score CANOPE CNV Calls (HIGH / MEDIUM / LOW)
#'
#' Uses the phred-like \code{Q_SOME} quality score, the number of reference
#' samples (\code{NUM_REFS}), the most-likely copy number (\code{MLCN}),
#' and a list of genes known to be difficult to score to assign one of three
#' confidence tiers.
#'
#' \strong{Assignment logic:}
#' \enumerate{
#'   \item HIGH: Q_SOME ≥ \code{high_q_score} AND NUM_REFS ≥ \code{high_num_refs}
#'         AND MLCN is consistent with CNV type (≤ 1 for DEL, ≥ 3 for DUP)
#'         AND gene not in \code{low_confidence_genes}.
#'   \item MEDIUM: Q_SOME ≥ \code{med_q_score} AND NUM_REFS ≥ \code{med_num_refs}
#'         AND gene not in \code{low_confidence_genes}.
#'   \item LOW: all other calls, including those in \code{low_confidence_genes}.
#' }
#'
#' @param cnv_calls         Data frame of CANOPE CNV calls.
#' @param high_q_score      Minimum Q_SOME for HIGH confidence. Default 50.
#' @param med_q_score       Minimum Q_SOME for MEDIUM confidence. Default 20.
#' @param high_num_refs     Minimum NUM_REFS for HIGH confidence. Default 10.
#' @param med_num_refs      Minimum NUM_REFS for MEDIUM confidence. Default 5.
#' @param high_num_targ     Minimum NUM_TARG for HIGH confidence. Default 1.
#' @param med_num_targ      Minimum NUM_TARG for MEDIUM confidence. Default 1.
#' @param low_confidence_genes Character vector of gene symbols always scored LOW.
#'   Default includes genes with known pseudogene / homology issues.
#'
#' @return The input data frame with an added \code{Confidence} character column
#'   and a \code{CN_label} column encoding the copy-number state as a string.
#' @export
score_canope_confidence <- function(
    cnv_calls,
    high_q_score         = 50L,
    med_q_score          = 20L,
    high_num_refs        = 10L,
    med_num_refs         = 5L,
    high_num_targ        = 1L,
    med_num_targ         = 1L,
    low_confidence_genes = c("PMS2", "SMN1", "CYP2D6", "HBA1", "HBA2",
                             "STRC", "CYP21A2", "GBA1", "CFTR")
) {
  if (is.null(cnv_calls) || nrow(cnv_calls) == 0) return(cnv_calls)

  required <- c("CNV", "Q_SOME", "NUM_REFS", "MLCN", "NUM_TARG", "GENE")
  missing_cols <- setdiff(required, colnames(cnv_calls))
  if (length(missing_cols) > 0) {
    warning("[WARNING] score_canope_confidence: missing columns: ",
            paste(missing_cols, collapse = ", "), ". Scoring skipped.")
    cnv_calls$Confidence <- "LOW"
    cnv_calls$CN_label   <- NA_character_
    return(cnv_calls)
  }

  q_score  <- as.integer(cnv_calls$Q_SOME)
  num_refs <- as.integer(cnv_calls$NUM_REFS)
  num_targ <- as.integer(cnv_calls$NUM_TARG)
  mlcn     <- as.integer(cnv_calls$MLCN)
  cnv_type <- toupper(as.character(cnv_calls$CNV)) 
  genes    <- as.character(cnv_calls$GENE)

  # MLCN must agree with the call direction. (A homozygous deletion, MLCN 0,
  # is already covered by the DEL branch; it must NOT count as consistent for a
  # DUP call.) A missing MLCN is never "consistent".
  mlcn_consistent <- (cnv_type == "DEL" & mlcn <= 1L) |
                     (cnv_type == "DUP" & mlcn >= 3L)
  mlcn_consistent[is.na(mlcn_consistent)] <- FALSE

  gene_flagged <- vapply(genes, function(g) {
    syms <- trimws(unlist(strsplit(g, "[;, ]")))
    any(syms %in% low_confidence_genes)
  }, logical(1))

  cnv_calls$Confidence <- "LOW"

  high_cond <- !gene_flagged &
    !is.na(q_score)  & q_score  >= high_q_score  &
    !is.na(num_refs) & num_refs >= high_num_refs  &
    !is.na(num_targ) & num_targ >= high_num_targ  &
    mlcn_consistent

  med_cond  <- !gene_flagged & !high_cond &
    !is.na(q_score)  & q_score  >= med_q_score   &
    !is.na(num_refs) & num_refs >= med_num_refs   &
    !is.na(num_targ) & num_targ >= med_num_targ

  cnv_calls$Confidence[high_cond] <- "HIGH"
  cnv_calls$Confidence[med_cond]  <- "MEDIUM"

  cnv_calls$CN_label <- dplyr::case_when(
    is.na(mlcn) ~ NA_character_,
    mlcn == 0L ~ "Hom. Del (CN=0)",
    mlcn == 1L ~ "Hem. Del (CN=1)",
    mlcn == 2L ~ "Normal (CN=2)",
    mlcn == 3L ~ "Dup (CN=3)",
    mlcn  > 3L ~ paste0("High-level Dup (CN=", mlcn, ")"),
    TRUE        ~ as.character(mlcn)
  )

  cnv_calls
}
