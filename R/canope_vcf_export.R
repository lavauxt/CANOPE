#' Export CANOPE CNV calls to VCF
#'
#' Converts a CANOPE CNV-call data frame (as returned by \code{\link{run_canope}}
#' and optionally annotated by \code{\link{score_canope_confidence}}) into a
#' standard VCF file.  When \code{sample_name} is \code{NULL}, a single
#' multi-sample VCF is written; otherwise a single-sample VCF is produced.
#'
#' @param cnv_calls   Data frame of CANOPE CNV calls (the \code{xcnv}-style
#'   table produced by \code{print_cnvs()} / \code{run_canope()}).
#' @param output_vcf  Output VCF file path.
#' @param sample_name Optional character string. If supplied, write a
#'   single-sample VCF for that sample only.
#' @param source      Source description for the VCF header. Default \code{"CANOPE"}.
#' @param reference   Reference genome assembly label (e.g. \code{"hg19"}).
#'
#' @return Invisibly returns the path written.
#' @export
export_canope_to_vcf <- function(cnv_calls, output_vcf, sample_name = NULL,
                                 source = "CANOPE", reference = "hg19") {
  if (is.null(cnv_calls) || nrow(cnv_calls) == 0) {
    message("[INFO] No CNVs to export.")
    return(invisible(NULL))
  }

  if (!is.null(sample_name)) {
    cnv_calls <- cnv_calls[cnv_calls$SAMPLE == sample_name, , drop = FALSE]
    if (nrow(cnv_calls) == 0) {
      warning("No CNVs found for sample ", sample_name)
      return(invisible(NULL))
    }
  }

  # ── Parse INTERVAL into Start/End if not already present ──────────────────
  if (!all(c("Start", "End") %in% colnames(cnv_calls))) {
    parsed <- lapply(cnv_calls$INTERVAL, parse_canope_interval)
    cnv_calls$Start <- vapply(parsed, function(p) p$start, integer(1))
    cnv_calls$End   <- vapply(parsed, function(p) p$end,   integer(1))
  }

  safe_char <- function(x, default = ".") {
    if (length(x) == 0 || is.na(x) || x == "NA") default else as.character(x)
  }
  samples <- sort(unique(as.character(cnv_calls$SAMPLE)))

  header <- c(
    "##fileformat=VCFv4.2",
    paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
    paste0("##source=", source),
    paste0("##reference=", reference),
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant (DUP/DEL)">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=CONFIDENCE,Number=1,Type=String,Description="Confidence level (HIGH/MEDIUM/LOW)">',
    '##INFO=<ID=NUMTARG,Number=1,Type=Integer,Description="Number of targets (exons) spanned">',
    '##INFO=<ID=REFS,Number=.,Type=String,Description="Reference samples used (comma-separated)">',
    '##INFO=<ID=NCOMP,Number=1,Type=Integer,Description="Number of reference samples">',
    '##INFO=<ID=QSOME,Number=1,Type=Integer,Description="Phred-scaled probability some CNV exists in this interval (SQ)">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Most-likely copy number (MLCN)">',
    '##FORMAT=<ID=QS,Number=1,Type=Integer,Description="Per-sample Q_SOME phred score">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=DEL,Description="Deletion">',
    paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
           paste(samples, collapse = "\t"))
  )

  key_df <- unique(data.frame(
    Chromosome = as.character(cnv_calls$CHR),
    Start      = as.integer(cnv_calls$Start),
    End        = as.integer(cnv_calls$End),
    Type       = as.character(cnv_calls$CNV),
    stringsAsFactors = FALSE
  ))

  chrom_order <- function(chr) {
    stripped <- sub("^chr", "", chr)
    num      <- suppressWarnings(as.integer(stripped))
    ifelse(is.na(num), Inf, num)
  }
  key_df <- key_df[order(chrom_order(key_df$Chromosome),
                         key_df$Chromosome,
                         key_df$Start, key_df$End, key_df$Type), , drop = FALSE]

  make_sample_field <- function(row_df) {
    if (nrow(row_df) == 0) return("0/0:2:.")
    mlcn <- suppressWarnings(as.integer(row_df$MLCN[1]))
    if (is.na(mlcn)) mlcn <- 2L
    qs <- suppressWarnings(as.integer(row_df$Q_SOME[1]))
    qs_str <- if (is.na(qs)) "." else as.character(qs)
    gt <- if (mlcn == 0L) "1/1" else if (mlcn < 2L) "0/1" else if (mlcn > 2L) "0/1" else "0/0"
    paste0(gt, ":", mlcn, ":", qs_str)
  }

  records <- character(nrow(key_df))
  for (i in seq_len(nrow(key_df))) {
    key_row  <- key_df[i, ]
    chrom    <- as.character(key_row[["Chromosome"]])
    start    <- as.integer(key_row[["Start"]])
    end      <- as.integer(key_row[["End"]])
    cnv_type <- toupper(as.character(key_row[["Type"]]))
    is_dup   <- cnv_type == "DUP"
    alt      <- if (is_dup) "<DUP>" else "<DEL>"
    sv_type  <- if (is_dup) "DUP" else "DEL"
    sv_len   <- if (is_dup) end - start + 1L else -(end - start + 1L)

    row_hits <- cnv_calls[
      as.character(cnv_calls$CHR) == chrom &
      as.integer(cnv_calls$Start) == start &
      as.integer(cnv_calls$End)   == end   &
      toupper(as.character(cnv_calls$CNV)) == cnv_type, , drop = FALSE]
    representative <- row_hits[1, , drop = FALSE]

    id <- paste0(gsub("^chr", "", chrom), "_", start, "_", end, "_", sv_type)

    info_items <- c(
      paste0("SVTYPE=", sv_type),
      paste0("SVLEN=", sv_len),
      paste0("END=", end)
    )

    conf <- if ("Confidence" %in% colnames(representative)) safe_char(representative$Confidence) else "."
    if (conf != ".") info_items <- c(info_items, paste0("CONFIDENCE=", conf))

    if ("NUM_TARG" %in% colnames(representative)) {
      nt <- safe_char(representative$NUM_TARG)
      if (nt != ".") info_items <- c(info_items, paste0("NUMTARG=", nt))
    }
    if ("REF_SAMPLES" %in% colnames(representative)) {
      refs_val <- safe_char(representative$REF_SAMPLES, "")
      if (nzchar(refs_val)) {
        # REF_SAMPLES is stored ';'-joined elsewhere in the pipeline (CSV
        # output, HTML report). But INFO fields use ';' to separate distinct
        # keys and ',' to separate multiple values *within* one Number=.
        # field (per the VCF spec) -- embedding literal ';' here corrupted
        # this record's INFO field for any downstream VCF parser. Convert
        # to ','-joined just for this VCF value.
        refs_vcf <- gsub(";", ",", refs_val, fixed = TRUE)
        info_items <- c(info_items, paste0("REFS=", refs_vcf))
      }
    }
    if ("NUM_REFS" %in% colnames(representative)) {
      ncomp <- safe_char(representative$NUM_REFS)
      if (ncomp != ".") info_items <- c(info_items, paste0("NCOMP=", ncomp))
    }
    if ("Q_SOME" %in% colnames(representative)) {
      qs <- safe_char(representative$Q_SOME)
      if (qs != ".") info_items <- c(info_items, paste0("QSOME=", qs))
    }

    format_field  <- "GT:CN:QS"
    sample_fields <- vapply(samples, function(s) {
      sample_row <- row_hits[row_hits$SAMPLE == s, , drop = FALSE]
      make_sample_field(sample_row)
    }, character(1))

    records[i] <- paste(chrom, start, id, "N", alt, ".", "PASS",
                        paste(info_items, collapse = ";"),
                        format_field,
                        paste(sample_fields, collapse = "\t"),
                        sep = "\t")
  }

  dir.create(dirname(output_vcf), showWarnings = FALSE, recursive = TRUE)
  con <- file(output_vcf, "w")
  on.exit(close(con), add = TRUE)
  writeLines(header,  con)
  writeLines(records, con)
  message("[INFO] Wrote ", length(records), " CNV records to ", output_vcf)
  invisible(output_vcf)
}
