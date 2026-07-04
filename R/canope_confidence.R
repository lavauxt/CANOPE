# =============================================================================
# canope_confidence.R — CNV confidence scoring for CANOPE output
# Adapted from ECHO/cnv_calls.R::score_cnv_confidence.
#
# CANOPE CNV columns used here:
#   SAMPLE, CNV (DEL|DUP), INTERVAL, KB, CHR, MID_BP,
#   TARGETS, NUM_TARG, GENE, MLCN, Q_SOME, NUM_REFS, REF_SAMPLES
# =============================================================================

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

  # ── Guard: ensure required columns exist ──────────────────────────────────
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
  cnv_type <- toupper(as.character(cnv_calls$CNV))   # "DEL" or "DUP"
  genes    <- as.character(cnv_calls$GENE)

  # Is the MLCN consistent with the CNV call?
  mlcn_consistent <- (cnv_type == "DEL" & mlcn <= 1L) |
                     (cnv_type == "DUP" & mlcn >= 3L) |
                     mlcn == 0L                          # homozygous deletion always consistent

  # Gene flagged as difficult?
  gene_flagged <- vapply(genes, function(g) {
    syms <- trimws(unlist(strsplit(g, "[;, ]")))
    any(syms %in% low_confidence_genes)
  }, logical(1))

  # ── Score ──────────────────────────────────────────────────────────────────
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

  # ── Human-readable copy-number label ──────────────────────────────────────
  cnv_calls$CN_label <- dplyr::case_when(
    mlcn == 0L ~ "Hom. Del (CN=0)",
    mlcn == 1L ~ "Hem. Del (CN=1)",
    mlcn == 2L ~ "Normal (CN=2)",
    mlcn == 3L ~ "Dup (CN=3)",
    mlcn  > 3L ~ paste0("High-level Dup (CN=", mlcn, ")"),
    TRUE        ~ as.character(mlcn)
  )

  cnv_calls
}
