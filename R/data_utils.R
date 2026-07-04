#' @importFrom dplyr %>%

#' @noRd
clean_name <- function(x) {
  return(strsplit(tools::file_path_sans_ext(basename(x)), "\\.")[[1]][1])
}

#' Format Chromosome Labels for Output and Plotting
#'
#' Maps integer / bare labels (1-22, 23, 24, X, Y) to UCSC-style \code{chr*} names.
#'
#' @param chr Character or numeric vector of chromosome identifiers.
#' @return Character vector of \code{chr}-prefixed labels (except \code{chrM}).
#' @export
format_chr_label <- function(chr) {
  chr <- as.character(chr)
  dplyr::case_when(
    chr %in% c("23", "X", "chrX")  ~ "chrX",
    chr %in% c("24", "Y", "chrY")  ~ "chrY",
    chr %in% c("M", "MT", "chrM")   ~ "chrM",
    grepl("^chr", chr)              ~ chr,
    TRUE                            ~ paste0("chr", chr)
  )
}

#' Map Target IDs to Row Indices in a Counts Table
#'
#' @param target_ids    Integer vector of target identifiers.
#' @param target_vector \code{target} column from the counts / BED table.
#' @return Integer row indices (\code{NA} if a target ID is not found).
#' @export
targets_to_rows <- function(target_ids, target_vector) {
  rows <- match(as.integer(target_ids), as.integer(target_vector))
  if (any(is.na(rows))) {
    warning(sprintf(
      "[WARNING] %d target id(s) could not be matched to a row in the counts table (returning NA).",
      sum(is.na(rows))
    ), call. = FALSE)
  }
  rows
}

#' Compute GC Content Specifically from BED Regions
#'
#' @param fasta_file Path to an indexed (or indexable) reference FASTA.
#' @param bed_input  Path to a BED file, or a data frame with chromosome/start/end columns.
#' @return Data frame: chromosome, start, end, GENE, GC_CONTENT.
#' @export
compute_gc_from_bed <- function(fasta_file, bed_input) {
  if (is.character(bed_input) && length(bed_input) == 1 && file.exists(bed_input)) {
    bed_df <- utils::read.table(bed_input, header = FALSE, stringsAsFactors = FALSE)
  } else if (is.data.frame(bed_input)) {
    bed_df <- bed_input
  } else {
    stop("ERROR: 'bed_input' must be a valid file path string or a data.frame.")
  }

  if (all(c("V1", "V2", "V3") %in% colnames(bed_df))) {
    colnames(bed_df)[1:3] <- c("chromosome", "start", "end")
  }

  bed <- GenomicRanges::makeGRangesFromDataFrame(
    bed_df,
    start.field = "start",
    end.field = "end",
    seqnames.field = "chromosome",
    keep.extra.columns = TRUE
  )

  if (!file.exists(paste0(fasta_file, ".fai"))) {
    Rsamtools::indexFa(fasta_file)
  }

  genome_fa <- Rsamtools::FaFile(fasta_file)
  genome_seqinfo <- GenomeInfoDb::seqinfo(genome_fa)
  valid_seqnames <- GenomeInfoDb::seqnames(genome_seqinfo)

  bed <- GenomeInfoDb::keepSeqlevels(bed, intersect(GenomeInfoDb::seqlevels(bed), valid_seqnames), pruning.mode = "coarse")

  region_seqs <- Biostrings::getSeq(genome_fa, bed)
  gc_freq <- Biostrings::letterFrequency(region_seqs, letters = "GC", as.prob = TRUE)

  gene_names <- if (!is.null(bed$name)) {
    as.character(bed$name)
  } else {
    rep(NA, length(bed))
  }

  return(data.frame(
    chromosome = as.character(GenomeInfoDb::seqnames(bed)),
    start      = BiocGenerics::start(bed),
    end        = BiocGenerics::end(bed),
    GENE       = gene_names,
    GC_CONTENT = as.numeric(gc_freq)
  ))
}

#' Get Read Coverage from BAM files
#'
#' @param bam_files Character vector of BAM file paths.
#' @param bed_input Path to a BED file, or a \code{GRanges} object of targets.
#' @param single_end Logical or \code{NULL} (default). If \code{NULL}, auto-detect
#'   from the first BAM file via \code{Rsamtools::testPairedEndBam}.
#' @return Data frame of read counts (rows = targets, columns = BAM basenames).
#' @export
get_coverage_from_bams <- function(bam_files, bed_input, single_end = NULL) {
  if (is.character(bed_input)) {
    targets <- rtracklayer::import(bed_input, format = "BED")
  } else if (inherits(bed_input, "GRanges")) {
    targets <- bed_input
  } else {
    stop("ERROR: bed_input must be a path to a BED file or a GRanges object.")
  }

  if (is.null(single_end)) {
    is_paired <- tryCatch(
      Rsamtools::testPairedEndBam(bam_files[1]),
      error = function(e) {
        message("[WARNING] Could not auto-detect BAM pairing (", conditionMessage(e),
                "); assuming paired-end.")
        TRUE
      }
    )
    single_end <- !is_paired
    message("[INFO] Auto-detected BAM type: ", if (single_end) "single-end" else "paired-end")
  }

  bam_list <- Rsamtools::BamFileList(bam_files, yieldSize = 100000)

  se <- GenomicAlignments::summarizeOverlaps(
    features = targets, reads = bam_list, mode = "Union",
    singleEnd = single_end, ignore.strand = TRUE
  )

  counts_matrix <- as.data.frame(SummarizedExperiment::assay(se))
  colnames(counts_matrix) <- basename(bam_files)

  return(counts_matrix)
}


#' Parse a Single Megadepth \code{*.annotation.tsv} File
#'
#' @param ann_file Path to a \code{*.annotation.tsv} file written by megadepth.
#' @param n_expected Expected row count (the number of BED targets); used only
#'   to produce a clear error if the file doesn't match, since a silent length
#'   mismatch would misalign every downstream sample column.
#' @return Numeric vector, one value per BED region, in file order.
#' @noRd
.parse_megadepth_annotation <- function(ann_file, n_expected = NULL) {
  if (!file.exists(ann_file))
    stop("[ERROR] Expected megadepth output not found: ", ann_file,
         ". If you're on a recent megadepth version the output filename or ",
         "column layout may have changed -- run megadepth once by hand and ",
         "check Sys.glob(paste0(dirname(ann_file), '/*')) before relying on this.")

  ann <- utils::read.table(ann_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if (ncol(ann) < 4)
    stop("[ERROR] Unexpected megadepth annotation format in ", ann_file,
         " (expected >= 4 columns: chrom, start, end, ..., score; got ", ncol(ann), ").")

  if (!is.null(n_expected) && nrow(ann) != n_expected)
    stop(sprintf(
      "[ERROR] megadepth returned %d region(s) but %d were expected from the BED file. ",
      nrow(ann), n_expected),
      "Pass keep_order = TRUE, and confirm the BED file's chromosomes all exist in the BAM."
    )

  as.numeric(ann[[ncol(ann)]])
}


#' Get Per-Target Coverage from BAM Files via Megadepth (fast; native Windows binary)
#'
#' @param bam_files    Character vector of BAM file paths.
#' @param bed_input    Path to a BED file of target regions (chrom, start, end, ...).
#' @param op            \code{"sum"} (default; total covered bases per region --
#'   the closer analogue to a read count) or \code{"mean"} (mean depth).
#' @param threads      Integer threads per megadepth invocation (BAM decompression).
#' @param work_dir     Directory for megadepth's intermediate output files.
#'   Default a fresh temp directory.
#' @param keep_order   Logical. Passes \code{--keep-order} so output rows are
#'   guaranteed to match \code{bed_input}'s row order exactly, regardless of
#'   chromosome order in the BAM. Strongly recommended; default \code{TRUE}.
#' @param install_if_missing Logical. If \code{TRUE} (default), calls
#'   \code{megadepth::install_megadepth()} first (a no-op if already installed).
#'
#' @return Data frame of per-target coverage values (rows = BED targets in
#'   \code{bed_input} order, columns = BAM basenames without extension).
#' @export
get_coverage_from_bams_megadepth <- function(bam_files, bed_input,
                                             op = c("sum", "mean"),
                                             threads = 1L,
                                             work_dir = NULL,
                                             keep_order = TRUE,
                                             install_if_missing = TRUE) {
  op <- match.arg(op)
  if (!requireNamespace("megadepth", quietly = TRUE))
    stop("[ERROR] Install the Bioconductor 'megadepth' package first:\n",
         "  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
         "  BiocManager::install('megadepth')")

  if (install_if_missing) {
    tryCatch(megadepth::install_megadepth(),
            error = function(e) message(
              "[WARNING] install_megadepth() failed (", conditionMessage(e),
              "); continuing in case megadepth is already on PATH."
            ))
  }

  if (!file.exists(bed_input)) stop("[ERROR] bed_input not found: ", bed_input)
  bed_df <- utils::read.table(bed_input, header = FALSE, stringsAsFactors = FALSE)
  n_targets <- nrow(bed_df)

  if (is.null(work_dir)) work_dir <- file.path(tempdir(), "megadepth_canope")
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

  sample_names <- vapply(bam_files, clean_name, character(1), USE.NAMES = FALSE)
  result_list  <- vector("list", length(bam_files))

  for (i in seq_along(bam_files)) {
    bam <- bam_files[i]
    if (!file.exists(bam)) stop("[ERROR] BAM not found: ", bam)

    prefix <- file.path(work_dir, paste0("md_", sample_names[i]))
    args <- c(
      "--annotation", shQuote(bed_input),
      "--op", op,
      "--no-annotation-stdout",
      "--prefix", shQuote(prefix),
      "--threads", as.character(threads)
    )
    if (keep_order) args <- c(args, "--keep-order")

    message(sprintf("[INFO] megadepth: %s (%d/%d)", sample_names[i], i, length(bam_files)))
    tryCatch(
      megadepth::megadepth_cmd(bam, paste(args, collapse = " ")),
      error = function(e) stop("[ERROR] megadepth failed on ", bam, ": ", conditionMessage(e))
    )

    ann_file <- paste0(prefix, ".annotation.tsv")
    result_list[[i]] <- .parse_megadepth_annotation(ann_file, n_expected = n_targets)
  }

  out_df <- as.data.frame(result_list)
  colnames(out_df) <- sample_names
  out_df
}
