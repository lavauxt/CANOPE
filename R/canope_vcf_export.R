#' Export CANOPE CNV calls to VCF (plain-text writer)
#'
#' Converts a CANOPE CNV-call data frame into a VCF 4.2 file with the writer
#' shared with ECHO (\code{vcf_write()}): UTF-8, LF line endings, unpadded
#' numbers, and a valid header-only file when nothing was called. A sample
#' without a call in a record is written as \code{./.} (no call is not evidence
#' of a confident reference genotype).
#'
#' @param cnv_calls   Data frame of CANOPE CNV calls (zero rows allowed).
#' @param output_vcf  Output VCF file path.
#' @param sample_name Optional single sample name (single-sample VCF).
#' @param source      Source description for the VCF header.
#' @param reference   Reference genome assembly label.
#' @param sample_names Optional analysed samples to list in a combined VCF even
#'   when they have no call (default: the samples that have calls).
#' @return Invisibly returns the path written.
#' @export
export_canope_to_vcf <- function(cnv_calls, output_vcf, sample_name = NULL,
                                 source = "CANOPE", reference = "hg19",
                                 sample_names = NULL) {
  if (is.null(cnv_calls)) cnv_calls <- .empty_cnv_calls()
  cnv_calls <- as.data.frame(cnv_calls, stringsAsFactors = FALSE)

  if (!is.null(sample_name)) {
    cnv_calls <- cnv_calls[cnv_calls$SAMPLE == sample_name, , drop = FALSE]
    vcf_samples <- as.character(sample_name)
  } else if (length(sample_names)) {
    vcf_samples <- unique(as.character(sample_names))
  } else {
    vcf_samples <- sort(unique(as.character(cnv_calls$SAMPLE)))
  }
  if (anyNA(vcf_samples) || any(!nzchar(vcf_samples)) || any(grepl("[[:space:]]", vcf_samples)))
    stop("[ERROR] VCF sample names must be non-empty and cannot contain whitespace.")
  absent <- setdiff(as.character(cnv_calls$SAMPLE), vcf_samples)
  if (length(absent))
    stop("[ERROR] CNV calls contain sample(s) absent from the VCF header: ", paste(absent, collapse = ", "))

  if (nrow(cnv_calls) > 0L && !all(c("Start", "End") %in% colnames(cnv_calls))) {
    parsed <- lapply(cnv_calls$INTERVAL, parse_canope_interval)
    cnv_calls$Start <- vapply(parsed, function(p) as.numeric(p$start), numeric(1))
    cnv_calls$End   <- vapply(parsed, function(p) as.numeric(p$end),   numeric(1))
  }

  safe_char <- function(x, default = ".") {
    if (length(x) == 0L || is.na(x[1]) || as.character(x[1]) == "NA") default else as.character(x[1])
  }
  meta <- c(
    "##fileformat=VCFv4.2",
    paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
    paste0("##source=", gsub("[\"\\\\]", "_", source)),
    paste0("##reference=", gsub("[\"\\\\]", "_", reference)),
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant (DUP/DEL)">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=CONFIDENCE,Number=1,Type=String,Description="Best confidence level (HIGH/MEDIUM/LOW) among the carrier samples; see FORMAT/CF per sample">',
    '##INFO=<ID=NUMTARG,Number=1,Type=Integer,Description="Number of targets (exons) spanned">',
    '##INFO=<ID=REFS,Number=.,Type=String,Description="Reference samples used (comma-separated); only written when a single sample carries the call">',
    '##INFO=<ID=NCOMP,Number=1,Type=Integer,Description="Largest number of reference samples among the carrier samples">',
    '##INFO=<ID=QSOME,Number=1,Type=Integer,Description="Highest phred-scaled probability some CNV exists in this interval (SQ) among the carrier samples; see FORMAT/QS per sample">',
    '##INFO=<ID=CARRIERS,Number=1,Type=Integer,Description="Number of samples of the cohort carrying an overlapping call of the same type">',
    '##INFO=<ID=RECURRENT,Number=0,Type=Flag,Description="Call recurs in many samples of the cohort (likely a technical artefact or a common CNV)">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Most-likely copy number (MLCN)">',
    '##FORMAT=<ID=QS,Number=1,Type=Integer,Description="Per-sample Q_SOME phred score">',
    '##FORMAT=<ID=CF,Number=1,Type=String,Description="Per-sample confidence level (HIGH/MEDIUM/LOW)">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=DEL,Description="Deletion">')

  key_df <- if (nrow(cnv_calls)) unique(data.frame(
    Chromosome = as.character(cnv_calls$CHR), Start = as.numeric(cnv_calls$Start),
    End = as.numeric(cnv_calls$End), Type = toupper(as.character(cnv_calls$CNV)),
    stringsAsFactors = FALSE)) else
    data.frame(Chromosome = character(), Start = numeric(), End = numeric(), Type = character(),
               stringsAsFactors = FALSE)
  if (nrow(key_df))
    key_df <- key_df[order(vcf_chrom_order(key_df$Chromosome), key_df$Chromosome,
                           key_df$Start, key_df$End, key_df$Type), , drop = FALSE]

  sample_field <- function(row_df) {
    if (nrow(row_df) == 0L) return("./.:.:.:.")
    mlcn <- suppressWarnings(as.integer(row_df$MLCN[1]))
    qs   <- suppressWarnings(as.integer(row_df$Q_SOME[1]))
    qs_str <- if (is.na(qs)) "." else as.character(qs)
    cf <- if ("Confidence" %in% colnames(row_df)) safe_char(row_df$Confidence[1]) else "."
    if (is.na(mlcn)) return(paste0("./.:.:", qs_str, ":", cf))
    gt <- if (mlcn == 0L) "1/1" else if (mlcn == 2L) "0/0" else "0/1"
    paste0(gt, ":", mlcn, ":", qs_str, ":", cf)
  }

  n <- nrow(key_df)
  chrom_v <- as.character(key_df$Chromosome); pos_v <- key_df$Start
  id_v <- alt_v <- info_v <- character(n)
  fields <- matrix("", nrow = n, ncol = length(vcf_samples), dimnames = list(NULL, vcf_samples))
  rank <- c(HIGH = 3L, MEDIUM = 2L, LOW = 1L)

  for (i in seq_len(n)) {
    k <- key_df[i, ]
    hits <- cnv_calls[as.character(cnv_calls$CHR) == k$Chromosome &
                        as.numeric(cnv_calls$Start) == k$Start & as.numeric(cnv_calls$End) == k$End &
                        toupper(as.character(cnv_calls$CNV)) == k$Type, , drop = FALSE]
    rep1 <- hits[1, , drop = FALSE]
    len  <- k$End - k$Start + 1
    sv_len <- if (k$Type == "DUP") len else -len
    id_v[i]  <- paste0(gsub("^chr", "", k$Chromosome, ignore.case = TRUE), "_", format(k$Start, scientific = FALSE, trim = TRUE),
                       "_", format(k$End, scientific = FALSE, trim = TRUE), "_", k$Type)
    alt_v[i] <- paste0("<", k$Type, ">")
    info <- c(paste0("SVTYPE=", k$Type), paste0("SVLEN=", format(sv_len, scientific = FALSE, trim = TRUE)),
              paste0("END=", format(k$End, scientific = FALSE, trim = TRUE)))
    if ("Confidence" %in% colnames(hits)) {
      cf <- toupper(as.character(hits$Confidence)); cf <- cf[cf %in% names(rank)]
      if (length(cf)) info <- c(info, paste0("CONFIDENCE=", cf[which.max(rank[cf])]))
    }
    if ("NUM_TARG" %in% colnames(rep1) && safe_char(rep1$NUM_TARG) != ".")
      info <- c(info, paste0("NUMTARG=", safe_char(rep1$NUM_TARG)))
    if ("REF_SAMPLES" %in% colnames(rep1) && nrow(hits) == 1L) {
      refs <- safe_char(rep1$REF_SAMPLES, "")
      if (nzchar(refs))
        info <- c(info, paste0("REFS=", paste(vapply(strsplit(refs, ";", fixed = TRUE)[[1]], vcf_info_escape, ""), collapse = ",")))
    }
    if ("NUM_REFS" %in% colnames(hits)) {
      nr <- suppressWarnings(max(as.integer(hits$NUM_REFS), na.rm = TRUE))
      if (is.finite(nr)) info <- c(info, paste0("NCOMP=", nr))
    }
    if ("Q_SOME" %in% colnames(hits)) {
      qm <- suppressWarnings(max(as.integer(hits$Q_SOME), na.rm = TRUE))
      if (is.finite(qm)) info <- c(info, paste0("QSOME=", qm))
    }
    if ("Carriers" %in% colnames(rep1) && !is.na(rep1$Carriers[1]))
      info <- c(info, paste0("CARRIERS=", as.integer(rep1$Carriers[1])))
    if ("Recurrent" %in% colnames(rep1) && isTRUE(rep1$Recurrent[1])) info <- c(info, "RECURRENT")
    info_v[i] <- paste(info, collapse = ";")
    for (s in vcf_samples) fields[i, s] <- sample_field(hits[hits$SAMPLE == s, , drop = FALSE])
  }

  records <- data.frame(CHROM = chrom_v, POS = pos_v, ID = id_v, ALT = alt_v, INFO = info_v,
                        stringsAsFactors = FALSE)
  vcf_write(output_vcf, meta, vcf_samples, records, format = "GT:CN:QS:CF", sample_fields = fields)
  invisible(output_vcf)
}
