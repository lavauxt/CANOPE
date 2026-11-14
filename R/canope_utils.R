

#' Pad a BED \emph{file}'s gene-terminal exons (file-level wrapper)
#'
#' Internally converts the (0-based, half-open) BED coordinates to 1-based
#' inclusive before padding, and back to 0-based on write-out, so the
#' padding math is identical to ECHO's (which operates on 1-based BED
#' output from \code{process_bed_file()}).
#'
#' @param input_bed Path to a BED file (>= 3 columns; a 4th "name"/gene
#'   column is required for per-gene padding to be meaningful -- without
#'   one, every row is treated as its own single-exon "gene").
#' @param output_bed Path to write the padded BED to.
#' @param padding Integer >= 0 bp of padding (see \code{\link{pad_gene_terminal_exons}}).
#'   \code{0} just copies \code{input_bed} to \code{output_bed} unchanged.
#' @param chr_lengths Optional named numeric vector of contig lengths (see
#'   \code{\link{pad_gene_terminal_exons}}).
#' @param verbose Logical. Default \code{TRUE}.
#' @return Invisibly returns \code{output_bed}.
#' @export
pad_bed_file <- function(input_bed, output_bed, padding = 0, chr_lengths = NULL, verbose = TRUE) {
  if (is.null(padding) || is.na(padding) || padding <= 0) {
    if (!identical(normalizePath(input_bed, mustWork = FALSE),
                   normalizePath(output_bed, mustWork = FALSE))) {
      file.copy(input_bed, output_bed, overwrite = TRUE)
    }
    return(invisible(output_bed))
  }

  raw <- utils::read.table(input_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (ncol(raw) < 3) stop("[ERROR] pad_bed_file: input BED must have at least 3 columns: ", input_bed)
  colnames(raw)[1:3] <- c("chromosome", "start", "end")
  if (ncol(raw) >= 4) colnames(raw)[4] <- "GENE" else raw$GENE <- paste0("Target_", seq_len(nrow(raw)))

  raw$start <- raw$start + 1L
  padded <- pad_gene_terminal_exons(raw, padding = padding, chr_lengths = chr_lengths, verbose = verbose)
  padded$start <- padded$start - 1L
  keep   <- if (ncol(raw) >= 4) colnames(raw) else setdiff(colnames(raw), "GENE")
  out_df <- padded[, keep, drop = FALSE]

  dir.create(dirname(output_bed), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(out_df, file = output_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  if (verbose) message("[INFO] pad_bed_file: padded BED written: ", output_bed)
  invisible(output_bed)
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


#' Infer each sample's sex from sex-chromosome coverage
#'
#' Heuristic on the depth of chrX (and chrY when present) relative to the
#' autosomal median: a female has ~1x chrX and ~0x chrY, a male ~0.5x chrX and
#' a clearly non-zero chrY. Only meant to flag an annotation / \code{modechrom}
#' mismatch, not to replace a proper sex check.
#'
#' @param chrom_depth Numeric matrix, rows = chromosomes (named), columns =
#'   samples: median target depth on the chromosome divided by the sample's
#'   autosomal median.
#' @param x_threshold chrX ratio below which a sample looks male. Default 0.75.
#' @param y_threshold chrY ratio at/above which a sample looks male. Default 0.1.
#' @return Data frame \code{sample, x_ratio, y_ratio, inferred_sex}
#'   (\code{"F"}, \code{"M"} or \code{"unclear"}), or \code{NULL} when the
#'   table has no chrX.
#' @export
canope_infer_sex <- function(chrom_depth, x_threshold = 0.75, y_threshold = 0.1)
  cnv_infer_sex(chrom_depth, x_threshold, y_threshold)


#' Guess the genome build label from a BSgenome package / FASTA name
#' @param ... Character strings to inspect (BSgenome package, FASTA path,
#'   \code{genome_version}).
#' @return \code{"hg38"}, \code{"hg19"}, \code{"T2T-CHM13"} or \code{"unknown"}.
#' @export
canope_infer_genome_build <- function(...) {
  x <- tolower(paste(unlist(list(...)), collapse = " "))
  if (grepl("hg38|grch38|b38", x)) "hg38"
  else if (grepl("hg19|grch37|b37", x)) "hg19"
  else if (grepl("t2t|chm13", x)) "T2T-CHM13"
  else "unknown"
}



#' Chromosome classes a CANOPE run analyses
#'
#' \code{"A"} analyses the autosomes, \code{"XX"} chrX only and \code{"XY"}
#' chrX (plus chrY when \code{removeY = FALSE}).
#' @param modechrom,removeY The \code{run_canope()} settings.
#' @return Character vector of classes (\code{"autosome"}, \code{"X"}, \code{"Y"}).
#' @noRd
canope_analysed_classes <- function(modechrom = "A", removeY = TRUE) {
  switch(as.character(modechrom),
         XX = "X", XY = if (isTRUE(removeY)) "X" else c("X", "Y"), "autosome")
}

#' Which targets a CANOPE run analyses (logical mask over chromosome labels)
#' @noRd
canope_mode_mask <- function(chromosome, modechrom = "A", removeY = TRUE) {
  chromosome_class(chromosome) %in% canope_analysed_classes(modechrom, removeY)
}

#' Chromosome scope of a CANOPE run (shared layout, see \code{chromosome_scope()})
#' @noRd
canope_scope <- function(chromosomes, modechrom = "A", removeY = TRUE, analysed_chromosomes = NULL) {
  why <- sprintf("not analysed with modechrom = '%s'", modechrom)
  chromosome_scope(chromosomes, analysed_classes = canope_analysed_classes(modechrom, removeY),
                   analysed_chromosomes = analysed_chromosomes,
                   reasons = c(autosome = why, X = why,
                               Y = if (identical(modechrom, "XY") && isTRUE(removeY)) "excluded by removeY = TRUE" else why))
}

# Hints shown when a chromosome class was skipped by the selected modechrom.
CANOPE_SCOPE_HINTS <- c(autosome = "modechrom = 'A'",
                        X = "modechrom = 'XX' (females) or 'XY' (males)",
                        Y = "modechrom = 'XY' with removeY = FALSE")
