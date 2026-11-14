#' Compute CANOPE QC Metrics
#'
#' Thin wrapper around the QC code shared with ECHO (\code{cnv_qc_metrics()},
#' \code{cnv_apply_qc()}): the SAME rules, thresholds and table schema
#' (\code{Sample, Exon, Type, Details, Gene, Severity}), computed on the RAW
#' counts of every candidate sample, so a sample or target that was excluded
#' still shows the numbers that excluded it.
#'
#' @param counts       Data frame / matrix of RAW counts (rows = targets,
#'   columns include \code{sample_names}).
#' @param sample_names Character vector of sample columns (test + references).
#' @param bed_df       Data frame aligned with the rows of \code{counts}
#'   (chromosome, start, end, GENE).
#' @param output_file  Path for the TSV output.
#' @param min_corr,min_cov,min_total_reads,max_exon_cv QC thresholds.
#' @param all_chromosomes Every chromosome of the INPUT targets (before
#'   \code{modechrom} filtering), used to explain which chromosome classes were
#'   not analysed.
#' @param modechrom,removeY The \code{run_canope()} chromosome settings.
#' @param qc Result of \code{cnv_apply_qc()} for this run. When \code{NULL} it
#'   is computed here from \code{qc_zscore}, \code{sample_qc} and the default
#'   target rules, so the function can also be used on its own.
#' @param excluded_samples Optional data frame (\code{sample}, \code{reason}) of
#'   samples removed before calling; used when \code{qc} is not supplied.
#' @param removed_targets Optional data frame (\code{label}, \code{gene},
#'   \code{reason}) of targets dropped before \code{counts} was built.
#' @param qc_zscore,sample_qc See \code{\link{run_canope}}.
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
                                  qc               = NULL,
                                  excluded_samples = NULL,
                                  removed_targets  = NULL,
                                  qc_zscore        = 3,
                                  sample_qc        = TRUE) {
  message("[INFO] BEGIN CANOPE QC Metrics")
  mat <- as.matrix(counts[, sample_names, drop = FALSE])
  storage.mode(mat) <- "double"
  bc  <- bed_cols(bed_df)
  chr_here <- as.character(bed_df[[bc$chrom]])

  if (is.null(qc)) {
    if (!is.null(excluded_samples)) {
      qc <- list(excluded_samples = excluded_samples, sample_qc = NULL, target_qc = NULL)
    } else {
      qc <- cnv_apply_qc(mat, chromosomes = chr_here, sample_qc = sample_qc, exon_qc = FALSE,
                         qc_zscore = qc_zscore)
    }
  }

  scope <- canope_scope(all_chromosomes %||% chr_here, modechrom, removeY, analysed_chromosomes = chr_here)
  scope_rows <- chromosome_scope_warnings(scope, "modechrom", CANOPE_SCOPE_HINTS)
  for (k in seq_len(nrow(scope_rows)))
    message(sprintf("[%s] %s", if (scope_rows$Severity[k] == "info") "INFO" else "WARNING", scope_rows$Details[k]))

  metrics <- cnv_qc_metrics(mat, bed_df, min_corr = min_corr, min_cov = min_cov,
                            min_total_reads = min_total_reads, max_exon_cv = max_exon_cv,
                            qc = qc, qc_zscore = qc_zscore, sample_qc = sample_qc,
                            removed_targets = removed_targets, scope_rows = scope_rows)
  cnv_write_qc(metrics, output_file)
  message("[INFO] QC metrics written to: ", output_file, " (", nrow(metrics), " flag(s))")
  message("[INFO] END CANOPE QC Metrics")
  invisible(metrics)
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
