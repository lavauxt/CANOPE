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

#' Compute GC Content from a FASTA File and a BED File
#'
#' Extracts per‑region GC fraction using a FASTA genome file (must be indexed).
#'
#' @param fasta_file Path to a FASTA file (with a `.fai` index).
#' @param bed_input  Path to a BED file, or a data frame with columns
#'   \code{chromosome, start, end} (and optionally \code{GENE}).
#' @return Data frame: chromosome, start, end, GENE, GC_CONTENT.
#' @export
compute_gc_from_fasta <- function(fasta_file, bed_input) {
  if (!requireNamespace("Rsamtools", quietly = TRUE))
    stop("Package 'Rsamtools' is required.")
  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Package 'Biostrings' is required.")

  # ---- Read BED ----
  if (is.character(bed_input) && length(bed_input) == 1 && file.exists(bed_input)) {
    bed_df <- utils::read.table(bed_input, header = FALSE, stringsAsFactors = FALSE)
  } else if (is.data.frame(bed_input)) {
    bed_df <- bed_input
  } else {
    stop("'bed_input' must be a BED file path or a data frame.")
  }

  # Ensure required columns
  if (all(c("V1", "V2", "V3") %in% colnames(bed_df))) {
    colnames(bed_df)[1:3] <- c("chromosome", "start", "end")
  }
  if (!all(c("chromosome", "start", "end") %in% colnames(bed_df)))
    stop("BED must contain chromosome, start, end columns.")

  # Sort
  if (requireNamespace("gtools", quietly = TRUE)) {
    bed_df <- bed_df[gtools::mixedorder(bed_df$chromosome), ]
  } else {
    bed_df <- bed_df[order(bed_df$chromosome, bed_df$start), ]
  }

  # ---- Open FASTA ----
  fa <- Rsamtools::FaFile(fasta_file)
  if (!file.exists(Rsamtools::index(fa)))
    stop("FASTA index (.fai) not found. Please run 'samtools faidx' on ", fasta_file)

  fa_seqnames <- as.character(GenomeInfoDb::seqnames(GenomicRanges::seqinfo(fa)))

  # ---- Map BED chromosome names to FASTA names ----
  bed_chroms <- unique(bed_df$chromosome)
  chrom_map <- setNames(rep(NA_character_, length(bed_chroms)), bed_chroms)

  for (bc in bed_chroms) {
    if (bc %in% fa_seqnames) {
      chrom_map[bc] <- bc
    } else {
      bc_no_chr <- sub("^chr", "", bc)
      if (bc_no_chr %in% fa_seqnames) {
        chrom_map[bc] <- bc_no_chr
      } else {
        bc_chr <- paste0("chr", bc)
        if (bc_chr %in% fa_seqnames) {
          chrom_map[bc] <- bc_chr
        } else {
          # Try numeric mapping for X/Y/M
          bc_num <- suppressWarnings(as.numeric(bc))
          if (!is.na(bc_num)) {
            if (bc_num == 23 && "X" %in% fa_seqnames) chrom_map[bc] <- "X"
            else if (bc_num == 24 && "Y" %in% fa_seqnames) chrom_map[bc] <- "Y"
            else if (bc_num == 25 && "M" %in% fa_seqnames) chrom_map[bc] <- "M"
          }
        }
      }
    }
  }

  unmapped <- names(chrom_map)[is.na(chrom_map)]
  if (length(unmapped) > 0) {
    message("[WARNING] Dropping chromosomes not found in FASTA: ",
            paste(unmapped, collapse = ", "))
    bed_df <- bed_df[!bed_df$chromosome %in% unmapped, ]
    chrom_map <- chrom_map[setdiff(names(chrom_map), unmapped)]
  }

  if (nrow(bed_df) == 0) stop("No BED regions remain after chromosome filtering.")

  # ---- Compute GC ----
  gc_vals <- numeric(nrow(bed_df))
  gene_names <- if ("GENE" %in% colnames(bed_df)) as.character(bed_df$GENE) else rep(NA, nrow(bed_df))

  for (i in seq_len(nrow(bed_df))) {
    chrom <- bed_df$chromosome[i]
    start <- bed_df$start[i] + 1L   # BED is 0‑based
    end   <- bed_df$end[i]
    fa_chrom <- chrom_map[chrom]

    tryCatch({
      gr <- GenomicRanges::GRanges(seqnames = fa_chrom,
                                   ranges = IRanges::IRanges(start = start, end = end))
      seqs <- Biostrings::getSeq(fa, gr)
      if (length(seqs) > 0 && Biostrings::width(seqs)[1] > 0) {
        gc_vals[i] <- Biostrings::letterFrequency(seqs, letters = "GC", as.prob = TRUE)[1]
      } else {
        gc_vals[i] <- NA_real_
      }
    }, error = function(e) {
      message("[WARNING] Failed to extract sequence for ", chrom, ":", start, "-", end)
      gc_vals[i] <<- NA_real_
    })
  }

  valid <- !is.na(gc_vals)
  if (sum(valid) == 0) stop("All GC computations failed.")
  if (sum(!valid) > 0) {
    message("[WARNING] Dropped ", sum(!valid), " region(s) with failed GC computation.")
    bed_df <- bed_df[valid, , drop = FALSE]
    gc_vals <- gc_vals[valid]
    gene_names <- gene_names[valid]
  }

  data.frame(
    chromosome = bed_df$chromosome,
    start      = bed_df$start,
    end        = bed_df$end,
    GENE       = gene_names,
    GC_CONTENT = as.numeric(gc_vals),
    stringsAsFactors = FALSE
  )
}

#' Compute GC Content Specifically from BED Regions using BSgenome
#'
#' @param bsgenome_pkg Character string naming the BSgenome package to use 
#'   (e.g., "BSgenome.Hsapiens.UCSC.hg38" or "BSgenome.Hsapiens.UCSC.hg19").
#' @param bed_input  Path to a BED file, or a data frame with chromosome/start/end columns.
#' @return Data frame: chromosome, start, end, GENE, GC_CONTENT.
#' @export
compute_gc_from_bed <- function(bsgenome_pkg, bed_input) {
  # ---- 0. Load the BSgenome package ----
  if (!is.character(bsgenome_pkg) || length(bsgenome_pkg) != 1) {
    stop("[ERROR] 'bsgenome_pkg' must be a single character string.")
  }
  
  if (!requireNamespace(bsgenome_pkg, quietly = TRUE)) {
    stop(sprintf(
      "[ERROR] The package '%s' is not installed. Please install it via BiocManager::install('%s').", 
      bsgenome_pkg, bsgenome_pkg
    ))
  }
  
  # Retrieve the actual BSgenome object dynamically
  genome_obj <- getExportedValue(bsgenome_pkg, bsgenome_pkg)
  
  # ---- 1. Read and sort the BED input ----
  if (is.character(bed_input) && length(bed_input) == 1 && file.exists(bed_input)) {
    bed_df <- utils::read.table(bed_input, header = FALSE, stringsAsFactors = FALSE)
  } else if (is.data.frame(bed_input)) {
    bed_df <- bed_input
  } else {
    stop("[ERROR] 'bed_input' must be a valid file path string or a data.frame.")
  }

  if (all(c("V1", "V2", "V3") %in% colnames(bed_df))) {
    colnames(bed_df)[1:3] <- c("chromosome", "start", "end")
  }

  # Sort BED by chromosome then start. NOTE: previously this re-sorted with
  # plain order() right after gtools::mixedorder(), which silently discarded
  # the natural numeric ordering (chr1, chr2, ..., chr10) in favour of plain
  # lexicographic ordering (chr1, chr10, chr11, ..., chr2, ...). Fixed to
  # actually keep the mixedorder() result when gtools is available.
  if (requireNamespace("gtools", quietly = TRUE)) {
    bed_df <- bed_df[gtools::mixedorder(bed_df$chromosome), ]
  } else {
    bed_df <- bed_df[order(bed_df$chromosome, bed_df$start), ]
  }
  message("[INFO] BED regions sorted by chromosome and start.")

  # ---- 2. Create GRanges (CRITICAL: BED is 0-based!) ----
  bed <- GenomicRanges::makeGRangesFromDataFrame(
    bed_df,
    start.field = "start",
    end.field = "end",
    seqnames.field = "chromosome",
    keep.extra.columns = TRUE,
    starts.in.df.are.0based = TRUE
  )

  # ---- 3. BSgenome sequence info and chromosome mapping ----
  genome_seqinfo <- GenomeInfoDb::seqinfo(genome_obj)
  fasta_chroms <- GenomeInfoDb::seqnames(genome_seqinfo)

  # Map BED chromosome names to BSgenome names
  bed_chroms <- GenomeInfoDb::seqlevels(bed)
  chrom_map <- setNames(rep(NA_character_, length(bed_chroms)), bed_chroms)

  for (bc in bed_chroms) {
    if (bc %in% fasta_chroms) {
      chrom_map[bc] <- bc
      next
    }
    bc_no_chr <- sub("^chr", "", bc)
    if (bc_no_chr %in% fasta_chroms) {
      chrom_map[bc] <- bc_no_chr
      next
    }
    bc_chr <- paste0("chr", bc)
    if (bc_chr %in% fasta_chroms) {
      chrom_map[bc] <- bc_chr
      next
    }
    bc_num <- suppressWarnings(as.numeric(bc))
    if (!is.na(bc_num)) {
      if (bc_num == 23 && "X" %in% fasta_chroms) {
        chrom_map[bc] <- "X"
        next
      }
      if (bc_num == 24 && "Y" %in% fasta_chroms) {
        chrom_map[bc] <- "Y"
        next
      }
      if (bc_num == 25 && "M" %in% fasta_chroms) {
        chrom_map[bc] <- "M"
        next
      }
    }
  }

  # Remove unmapped regions and drop unused seqlevels
  unmapped <- names(chrom_map)[is.na(chrom_map)]
  if (length(unmapped) > 0) {
    message("[WARNING] Dropping chromosomes not found in BSgenome: ",
            paste(unmapped, collapse = ", "))
    keep <- !as.character(GenomeInfoDb::seqnames(bed)) %in% unmapped
    bed <- bed[keep, ]
    bed <- GenomeInfoDb::dropSeqlevels(bed, unmapped)
    chrom_map <- chrom_map[setdiff(names(chrom_map), unmapped)]
  }

  if (length(bed) == 0) stop("[ERROR] No BED regions remain after chromosome filtering.")

  # Rename seqlevels to match BSgenome
  new_levels <- chrom_map[GenomeInfoDb::seqlevels(bed)]
  if (any(is.na(new_levels))) {
    stop("[ERROR] Internal error: some seqlevels are NA after mapping. Check mapping logic.")
  }
  GenomeInfoDb::seqlevels(bed) <- new_levels

  # Keep only seqlevels that are actually in the BSgenome
  keep_levels <- intersect(GenomeInfoDb::seqlevels(bed), fasta_chroms)
  if (length(keep_levels) == 0) {
    stop("[ERROR] No BED regions share seqlevels with the BSgenome.")
  }
  bed <- GenomeInfoDb::keepSeqlevels(bed, keep_levels, pruning.mode = "coarse")

  # ---- 4. GC content extraction ----
  gc_freq <- numeric(length(bed))
  gene_names <- if (!is.null(bed$name)) as.character(bed$name) else rep(NA, length(bed))

  for (i in seq_along(bed)) {
    region <- bed[i]
    seq_chrom <- as.character(GenomeInfoDb::seqnames(region))
    seq_start <- BiocGenerics::start(region)
    seq_end   <- BiocGenerics::end(region)

    tryCatch({
      seqs <- Biostrings::getSeq(genome_obj, region)  # returns DNAStringSet from BSgenome
      
      # Check that we have at least one non-empty sequence
      if (length(seqs) > 0 && all(Biostrings::width(seqs) > 0)) {
        # Use the whole DNAStringSet – avoids S4 coercion issues
        gc_freq[i] <- Biostrings::letterFrequency(seqs, letters = "GC", as.prob = TRUE)[1]
      } else {
        message(sprintf("[WARNING] Empty sequence for %s:%d-%d", seq_chrom, seq_start, seq_end))
        gc_freq[i] <- NA_real_
      }
    }, error = function(e) {
      message(sprintf(
        "[WARNING] Failed to read sequence for %s:%d-%d: %s",
        seq_chrom, seq_start, seq_end, conditionMessage(e)
      ))
      gc_freq[i] <<- NA_real_
    })
  }

  # Drop regions with NA GC content
  valid <- !is.na(gc_freq)
  if (sum(valid) == 0) {
    stop("[ERROR] All GC content computations failed. Check that the BSgenome ",
         "contains the positions in your BED file.")
  }
  if (sum(!valid) > 0) {
    message(sprintf("[WARNING] Dropped %d target(s) with failed GC computation.", sum(!valid)))
    bed <- bed[valid]
    gc_freq <- gc_freq[valid]
    gene_names <- gene_names[valid]
  }

  data.frame(
    chromosome = as.character(GenomeInfoDb::seqnames(bed)),
    start      = BiocGenerics::start(bed),
    end        = BiocGenerics::end(bed),
    GENE       = gene_names,
    GC_CONTENT = as.numeric(gc_freq),
    stringsAsFactors = FALSE
  )
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