.align_reference_counts <- function(target_data, reference_data) {
  coordinate_cols <- c("chromosome", "start", "end")
  missing_target <- setdiff(coordinate_cols, names(target_data))
  missing_reference <- setdiff(coordinate_cols, names(reference_data))
  if (length(missing_target) || length(missing_reference)) {
    stop("Both target and reference count tables must contain chromosome, start, and end columns.",
         call. = FALSE)
  }

  make_keys <- function(data) {
    chromosome <- as.character(data$chromosome)
    start <- suppressWarnings(as.numeric(as.character(data$start)))
    end <- suppressWarnings(as.numeric(as.character(data$end)))
    if (anyNA(chromosome) || anyNA(start) || anyNA(end))
      stop("Target coordinates must be non-missing and numeric where applicable.", call. = FALSE)
    paste(chromosome, start, end, sep = "\r")
  }

  target_keys <- make_keys(target_data)
  reference_keys <- make_keys(reference_data)
  if (anyDuplicated(target_keys) || anyDuplicated(reference_keys))
    stop("Duplicate chromosome/start/end coordinates prevent unambiguous reference alignment.",
         call. = FALSE)
  if (length(target_keys) != length(reference_keys) || !setequal(target_keys, reference_keys))
    stop("Target and reference count tables contain different genomic coordinates.", call. = FALSE)

  reference_data[match(target_keys, reference_keys), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Input helpers (each does one job; run_canope() only orchestrates them)
# ---------------------------------------------------------------------------

#' BAM files and sample names from a directory or a list file
#' @return List with \code{bams} and \code{names}.
#' @noRd
canope_resolve_bams <- function(samples) {
  if (length(samples) == 1L && dir.exists(samples)) {
    bams <- list.files(samples, pattern = "\\.bam$", full.names = TRUE)
    if (length(bams) == 0L) stop("[ERROR] No .bam files found in directory: ", samples, call. = FALSE)
  } else if (length(samples) == 1L && file.exists(samples)) {
    bams <- trimws(readLines(samples, warn = FALSE))
    bams <- bams[nzchar(bams)]
  } else {
    stop("[ERROR] 'samples' must be a directory of BAM files or a file listing BAM paths.", call. = FALSE)
  }
  nm <- vapply(bams, clean_name, character(1), USE.NAMES = FALSE)
  if (anyDuplicated(nm)) {
    dups <- unique(nm[duplicated(nm)])
    stop("[ERROR] Sample names derived from the BAM file names are not unique (",
         paste(utils::head(dups, 5), collapse = ", "),
         "): the name is the part of the file name before the first '.', so e.g. ",
         "'P1.A.bam' and 'P1.B.bam' both become 'P1'. Rename the BAM files.", call. = FALSE)
  }
  list(bams = bams, names = nm)
}

#' Read the per-target counts: a pre-computed table or fresh from the BAMs
#'
#' @return List: \code{reads} (chromosome, start, end, GENE, one column per
#'   sample), \code{bed_file} (the BED actually used -- the padded one when
#'   padding was applied) and \code{sample_names}.
#' @noRd
canope_read_input_counts <- function(reads_file, bed_file, bams, sample_names,
                                     pad_terminal_exons, fasta_file, bsgenome_pkg,
                                     coverage_backend, megadepth_op, megadepth_threads,
                                     out_dir, output_prefix, lg) {
  n_meta <- 4L
  if (!is.null(reads_file) && file.exists(reads_file)) {
    lg$info("Reading pre-computed counts: ", basename(reads_file))
    reads <- utils::read.table(reads_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(reads) <= n_meta)
      stop("[ERROR] reads_file has no sample column after the ", n_meta,
           " metadata columns (chromosome, start, end, GENE).", call. = FALSE)
    if (is.null(sample_names)) sample_names <- colnames(reads)[-seq_len(n_meta)]
    if (ncol(reads) - n_meta != length(sample_names)) {
      stop(sprintf("[ERROR] reads_file has %d column(s) after the %d metadata columns (chromosome, start, end, GENE), but %d sample(s) were supplied. ",
                   ncol(reads) - n_meta, n_meta, length(sample_names)),
           "Check that reads_file has exactly one column per sample after those metadata columns.",
           call. = FALSE)
    }
    # the first four columns are the metadata whatever the header spelling,
    # the others are the samples, in the order of `sample_names`
    colnames(reads) <- c("chromosome", "start", "end", "GENE", sample_names)
    return(list(reads = reads, bed_file = bed_file, sample_names = sample_names))
  }

  if (is.null(bed_file) || !file.exists(bed_file))
    stop("[ERROR] A bed_file is required to extract coverage from BAMs.", call. = FALSE)
  if (!length(bams)) stop("[ERROR] No BAM files: supply 'samples' (or a reads_file).", call. = FALSE)

  if (!is.null(pad_terminal_exons) && !is.na(pad_terminal_exons) && pad_terminal_exons > 0) {
    chr_lengths <- NULL
    if (!is.null(fasta_file) && file.exists(fasta_file)) {
      chr_lengths <- tryCatch(GenomeInfoDb::seqlengths(Rsamtools::FaFile(fasta_file)),
        error = function(e) { lg$warn("Could not read contig lengths from fasta_file for padding clamp: ", conditionMessage(e)); NULL })
    } else if (!is.null(bsgenome_pkg) && requireNamespace(bsgenome_pkg, quietly = TRUE)) {
      chr_lengths <- tryCatch(GenomeInfoDb::seqlengths(getExportedValue(bsgenome_pkg, bsgenome_pkg)),
        error = function(e) { lg$warn("Could not read contig lengths from bsgenome_pkg for padding clamp: ", conditionMessage(e)); NULL })
    }
    padded_bed <- file.path(out_dir, paste0(output_prefix, "_padded.bed"))
    lg$info(sprintf("Applying terminal-exon padding (%d bp) to %s", as.integer(pad_terminal_exons), basename(bed_file)))
    pad_bed_file(bed_file, padded_bed, padding = pad_terminal_exons, chr_lengths = chr_lengths)
    bed_file <- padded_bed
  }

  lg$info(sprintf("Extracting coverage from %d BAM file(s) (backend = '%s')", length(bams), coverage_backend))
  counts_df <- switch(coverage_backend,
    bioconductor = get_coverage_from_bams(bam_files = bams, bed_input = bed_file),
    megadepth    = get_coverage_from_bams_megadepth(bam_files = bams, bed_input = bed_file,
                                                    op = megadepth_op, threads = megadepth_threads))
  colnames(counts_df) <- sample_names

  coords_full <- as.data.frame(.read_bed_as_granges(bed_file))
  gene_names  <- if ("name" %in% colnames(coords_full)) coords_full$name
                 else paste0("Target_", seq_len(nrow(coords_full)))
  coords <- data.frame(chromosome = as.character(coords_full$seqnames), start = coords_full$start,
                       end = coords_full$end, GENE = gene_names, stringsAsFactors = FALSE)
  dup <- duplicated(coords[, c("chromosome", "start", "end")])
  if (any(dup)) lg$warn(sprintf("%d duplicated BED interval(s) removed", sum(dup)))
  list(reads = cbind(coords[!dup, , drop = FALSE], counts_df[!dup, , drop = FALSE]),
       bed_file = bed_file, sample_names = sample_names)
}

#' GC fraction of every target (aligned by coordinates)
#' @return Numeric vector, one value per row of \code{targets} (NA when the
#'   interval is not in the GC source).
#' @noRd
canope_target_gc <- function(targets, fasta_file, bsgenome_pkg, gc_file, bed_file, lg) {
  if (!is.null(fasta_file) && !is.null(bed_file)) {
    lg$info("Computing GC content from the FASTA file")
    datagc <- compute_gc_from_fasta(fasta_file = fasta_file, bed_input = bed_file)
  } else if (!is.null(bsgenome_pkg) && !is.null(bed_file)) {
    if (!requireNamespace(bsgenome_pkg, quietly = TRUE))
      stop("[ERROR] BSgenome package '", bsgenome_pkg, "' is not installed.", call. = FALSE)
    lg$info(sprintf("Computing GC content from BSgenome package '%s'", bsgenome_pkg))
    datagc <- compute_gc_from_bed(bsgenome_pkg = bsgenome_pkg, bed_input = bed_file)
  } else if (!is.null(gc_file) && file.exists(gc_file)) {
    lg$info("Using the pre-computed GC file")
    datagc <- utils::read.table(gc_file, header = TRUE, stringsAsFactors = FALSE)
  } else {
    stop("[ERROR] Provide (fasta_file + bed_file), (bsgenome_pkg + bed_file) or a pre-computed gc_file.", call. = FALSE)
  }
  if (!"GC_CONTENT" %in% colnames(datagc)) {
    guess <- grep("^gc", colnames(datagc), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(guess)) colnames(datagc)[colnames(datagc) == guess] <- "GC_CONTENT"
    else if (ncol(datagc) >= 4) colnames(datagc)[4] <- "GC_CONTENT"
    else stop("[ERROR] Could not identify a GC content column.", call. = FALSE)
  }
  gc <- .align_gc_to_targets(datagc, targets)
  if (any(gc > 1, na.rm = TRUE)) {
    lg$info("GC content is on a 0-100 scale; converting to a fraction")
    gc <- gc / 100
  }
  gc
}

#' Add the counts of an external reference panel to the cohort table
#' @return List: \code{reads}, \code{refsample_names}.
#' @noRd
canope_attach_reference_panel <- function(reads, refbams_file, ref_reads, modechrom, lg) {
  if (is.null(refbams_file) || !file.exists(refbams_file))
    return(list(reads = reads, refsample_names = character(0)))
  if (is.null(ref_reads) || !file.exists(ref_reads))
    stop("[ERROR] 'refbams_file' was supplied but 'ref_reads' is missing or does not exist: ",
         "'ref_reads' must be a pre-computed counts table for the reference panel.", call. = FALSE)
  raw <- utils::read.csv(refbams_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if ("gender" %in% colnames(raw)) {
    if (modechrom == "XX")      raw <- raw[raw$gender %in% "F", , drop = FALSE]
    else if (modechrom == "XY") raw <- raw[raw$gender %in% "M", , drop = FALSE]
  }
  ref_names <- tools::file_path_sans_ext(basename(apply(raw, 1, toString)))
  data_ref <- utils::read.table(ref_reads, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  n_meta <- 4L
  if (ncol(data_ref) <= n_meta)
    stop("[ERROR] ref_reads has no sample column after the ", n_meta,
         " metadata columns (chromosome, start, end, GENE).", call. = FALSE)
  missing_refs <- setdiff(ref_names, colnames(data_ref))
  if (length(missing_refs))
    stop(sprintf("[ERROR] %d reference sample(s) from refbams_file not found as columns in ref_reads: %s",
                 length(missing_refs), paste(missing_refs, collapse = ", ")), call. = FALSE)
  data_ref <- .align_reference_counts(reads, data_ref)
  # a panel sample that is also a test sample would be counted twice
  ref_names <- setdiff(ref_names, colnames(reads))
  lg$info(sprintf("External reference panel: %d sample(s)", length(ref_names)))
  list(reads = cbind(reads, data_ref[, ref_names, drop = FALSE]), refsample_names = ref_names)
}

#' Keep the analysed chromosomes, sort, number the targets
#'
#' Target ids increase with genomic position: the HMM, the Q_SOME likelihoods
#' and the plots treat \code{TARGETS = a..b} as a contiguous range. Exon numbers
#' are fixed here, BEFORE QC removes targets, so a label always means the
#' exon's true position in its gene.
#' @return Data frame with \code{target}, \code{chromosome}, \code{start},
#'   \code{end}, \code{GENE}, \code{exon_number} and the sample columns.
#' @noRd
canope_prepare_targets <- function(reads, modechrom, removeY, lg) {
  chr <- as.character(reads$chromosome)
  keep <- canope_mode_mask(chr, modechrom, removeY)
  if (!any(keep))
    stop(sprintf("[ERROR] No targets remain after modechrom = '%s' (removeY = %s). Chromosomes in the input: %s.",
                 modechrom, removeY, paste(utils::head(unique(chr), 8), collapse = ", ")), call. = FALSE)
  if (any(!keep))
    lg$info(sprintf("%d of %d target(s) are outside the chromosomes analysed with modechrom = '%s' (the input chromosome inventory is kept for the report)",
                    sum(!keep), length(keep), modechrom))
  reads <- reads[keep, , drop = FALSE]
  reads$chromosome <- as.character(reads$chromosome)
  reads$start <- as.numeric(reads$start)
  reads$end <- as.numeric(reads$end)
  reads <- reads[order(match(reads$chromosome, chrom_sort_levels(reads$chromosome)), reads$start, reads$end), , drop = FALSE]
  dup <- duplicated(reads[, c("chromosome", "start", "end")])
  if (any(dup)) {
    lg$warn(sprintf("%d duplicated target interval(s) removed", sum(dup)))
    reads <- reads[!dup, , drop = FALSE]
  }
  rownames(reads) <- NULL
  reads$GENE <- as.character(reads$GENE)
  reads$target <- seq_len(nrow(reads))
  reads$exon_number <- assign_exon_numbers_per_gene(reads[, c("chromosome", "start", "end", "GENE")])$exon_number
  reads
}

#' Call every analysed sample (one HMM run each)
#' @return List: \code{calls}, \code{models}, \code{refs}, \code{status},
#'   \code{low_coverage_runs}.
#' @noRd
canope_call_all <- function(reads, samples_to_call, refsample_names, auto_refs, auto_reference,
                            par, excluded_samples, lg) {
  all_cnvs <- list(); models <- list(); refs <- list(); status <- list(); low_cov <- list()
  set_status <- function(sample, status_name, message = "", n_cnvs = 0L, n_targets = NA_integer_,
                         n_low_count = NA_integer_, n_refs = NA_integer_, ref_cor = NA_real_,
                         genotyped = NA) {
    status[[sample]] <<- data.frame(
      sample = sample, status = status_name, message = paste(message, collapse = " | "),
      n_cnvs = as.integer(n_cnvs), n_targets = as.integer(n_targets),
      n_low_count_dropped = as.integer(n_low_count), n_refs = as.integer(n_refs),
      ref_cor_median = as.numeric(ref_cor), genotyped = genotyped, stringsAsFactors = FALSE)
  }
  if (nrow(excluded_samples) > 0)
    for (i in seq_len(nrow(excluded_samples)))
      set_status(excluded_samples$sample[i], "excluded", excluded_samples$reason[i])

  for (k in seq_along(samples_to_call)) {
    samp <- samples_to_call[k]
    lg$info(sprintf("Calling CNVs for sample %s (%d/%d)", samp, k, length(samples_to_call)))
    warn_msgs <- character(0)
    tryCatch({
      refs_to_use <- if (auto_reference && !is.null(auto_refs)) auto_refs[[samp]] else refsample_names
      if (auto_reference && is.null(refs_to_use))
        lg$warn(sprintf("Sample '%s' is not in the auto-inferred reference panel (it is among the noisiest samples); all other samples are candidate references.", samp))
      res <- withCallingHandlers(
        call_cnvs(sample_name = samp, counts = reads, p = par$p_value, Tnum = par$Tnum, D = par$D,
                  numrefs = par$numrefs, get_dfs = FALSE, homdel_mean = par$homdel_mean,
                  refsample_names = refs_to_use, full_output = TRUE, min_cor = par$min_cor,
                  decode_method = par$decode_method, engine = par$engine,
                  min_sample_count = par$min_sample_count),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          lg$warn(samp, ": ", conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      n_called <- if (!is.null(res$cnvs)) nrow(res$cnvs) else 0L
      lg$info(sprintf("  %d CNV(s) found for %s", n_called, samp))
      lg$debug(sprintf("  %s: %s target(s) used (%s dropped for low counts), %d reference(s), median reference correlation %.3f",
                       samp, format(res$n_targets_used), format(res$n_low_count_dropped), length(res$reference_samples),
                       stats::median(res$ref_cor, na.rm = TRUE)))
      if (n_called > 0) all_cnvs[[samp]] <- res$cnvs
      refs[[samp]] <- res$reference_samples
      models[[samp]] <- list(target = res$targets, mean = res$mean, var_estimate = res$var_estimate,
                             test_counts = res$test_counts, ref_matrix = res$ref_matrix)
      if (!is.null(res$low_coverage_runs) && nrow(res$low_coverage_runs) > 0) {
        low_cov[[samp]] <- cbind(SAMPLE = samp, res$low_coverage_runs, stringsAsFactors = FALSE)
        lg$warn(sprintf("%s: %d run(s) of >= 3 consecutive targets with ~no reads were excluded before calling (possible homozygous deletion / capture dropout; not callable with min_sample_count = %d).",
                        samp, nrow(res$low_coverage_runs), as.integer(par$min_sample_count)))
      }
      set_status(samp, "called", warn_msgs, n_cnvs = n_called, n_targets = res$n_targets_used,
                 n_low_count = res$n_low_count_dropped, n_refs = length(res$reference_samples),
                 ref_cor = stats::median(res$ref_cor, na.rm = TRUE), genotyped = isTRUE(res$genotyped))
    }, error = function(e) {
      lg$warn(sprintf("Sample %s FAILED: %s", samp, conditionMessage(e)))
      set_status(samp, "failed", c(conditionMessage(e), warn_msgs))
    })
  }
  list(calls = all_cnvs, models = models, refs = refs, status = status, low_coverage_runs = low_cov)
}

#' Read-level evidence of every call, in the same terms ECHO reports
#'
#' \code{READS_OBSERVED} / \code{READS_EXPECTED} are the summed counts of the
#' call's targets in the test sample and in the (NNLS-weighted) reference
#' model, \code{READS_RATIO} their quotient (ECHO: \code{Reads.ratio}) and
#' \code{REF_CORR} the median correlation of the sample with its references
#' (ECHO: \code{Correlation}).
#' @noRd
canope_call_read_stats <- function(calls, models, status = list()) {
  n <- nrow(calls)
  obs <- exp <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    m <- models[[as.character(calls$SAMPLE[i])]]
    tb <- parse_canope_targets(calls$TARGETS[i])
    if (is.null(m) || is.null(m$target) || length(tb) < 2L) next
    k <- which(m$target >= tb[1L] & m$target <= tb[2L])
    if (!length(k)) next
    obs[i] <- sum(m$test_counts[k], na.rm = TRUE)
    exp[i] <- sum(pmax(m$mean[k], 1), na.rm = TRUE)
  }
  rc <- vapply(as.character(calls$SAMPLE), function(s) {
    v <- status[[s]]$ref_cor_median
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, numeric(1))
  data.frame(READS_OBSERVED = round(obs, 1), READS_EXPECTED = round(exp, 1),
             READS_RATIO = round(obs / exp, 3), REF_CORR = round(unname(rc), 3))
}

#' Post-process the pooled calls: coordinates, scoring, terminal-only rule, recurrence
#'
#' Same sequence (and same options) as ECHO: score, flag the calls that only
#' touch the terminal exons of a gene (\code{remove_terminal_only} /
#' \code{penalize_terminal_only}), then the cohort-level recurrence flag.
#' @return Call table with \code{Start}, \code{End}, \code{Carriers},
#'   \code{Recurrent}, \code{terminal_only}, \code{global_start},
#'   \code{global_end} and \code{CNV.ID}.
#' @noRd
canope_finalise_calls <- function(calls, bed_file, n_samples, scored, confidence_args,
                                  remove_terminal_only, penalize_terminal_only,
                                  recurrence_action, recurrence_min_samples,
                                  recurrence_min_fraction, recurrence_min_overlap, lg,
                                  models = list(), status = list()) {
  if (nrow(calls) == 0L) {
    calls <- .empty_cnv_calls(scored = scored)
    calls$READS_OBSERVED <- numeric(0); calls$READS_EXPECTED <- numeric(0)
    calls$READS_RATIO <- numeric(0); calls$REF_CORR <- numeric(0)
    calls$Start <- integer(0); calls$End <- integer(0)
    calls$Carriers <- integer(0); calls$Recurrent <- logical(0); calls$terminal_only <- logical(0)
    calls$global_start <- integer(0); calls$global_end <- integer(0)
    return(cbind(CNV.ID = integer(0), calls))
  }
  iv <- lapply(as.character(calls$INTERVAL), parse_canope_interval)
  calls$Start <- vapply(iv, function(p) as.integer(p$start), integer(1))
  calls$End   <- vapply(iv, function(p) as.integer(p$end),   integer(1))
  calls <- cbind(calls, canope_call_read_stats(calls, models, status))
  tb <- lapply(as.character(calls$TARGETS), parse_canope_targets)
  # `bed_file` is the table of ALL analysed-chromosome targets (before QC), whose
  # target ids are its row numbers
  calls$global_start <- match(vapply(tb, function(t) as.integer(t[1L]), integer(1)), bed_file$target)
  calls$global_end   <- match(vapply(tb, function(t) as.integer(t[length(t)]), integer(1)), bed_file$target)

  if (scored) {
    lg$info("Scoring CNV confidence")
    calls <- do.call(score_canope_confidence, c(list(cnv_calls = calls), confidence_args))
  }
  calls <- flag_terminal_only_calls(calls, bed_file)
  n_term <- sum(calls$terminal_only, na.rm = TRUE)
  if (n_term > 0L) {
    if (isTRUE(remove_terminal_only)) {
      lg$info(sprintf("Removing %d call(s) that only affect terminal exon(s) of a gene", n_term))
      calls <- calls[!calls$terminal_only %in% TRUE, , drop = FALSE]
    } else if (isTRUE(penalize_terminal_only) && "Confidence" %in% colnames(calls)) {
      lg$info(sprintf("Downgrading %d terminal-exon-only call(s) to LOW confidence", n_term))
      calls$Confidence[calls$terminal_only %in% TRUE] <- "LOW"
    } else {
      lg$info(sprintf("%d call(s) only affect terminal exon(s) of a gene (flagged in the report; not removed)", n_term))
    }
  }
  # cohort-level recurrence (the same call in many samples is a noisy target,
  # not a germline CNV); applied once, on the whole cohort
  calls <- cnv_apply_recurrence(
    calls, list(sample = "SAMPLE", chr = "CHR", start = "Start", end = "End", type = "CNV",
                confidence = "Confidence"),
    n_samples = n_samples, action = recurrence_action, min_samples = recurrence_min_samples,
    min_fraction = recurrence_min_fraction, min_overlap = recurrence_min_overlap)
  rownames(calls) <- NULL
  cbind(CNV.ID = seq_len(nrow(calls)), calls)
}

#' Columns written to the CNV CSV (the others are helpers of the report / plots)
#' @noRd
CANOPE_CSV_EXCLUDE <- c("global_start", "global_end", "terminal_only")

#' Run the CANOPE Copy Number Calling Pipeline
#'
#' End-to-end wrapper: reads count/GC data, decides which samples and targets
#' are removed (with the same rules as ECHO, see \code{\link{cnv_apply_qc}}),
#' corrects the GC bias, infers a reference set per sample, runs the HMM, scores
#' and post-processes the calls, and writes the CSV, the RData workspace, the QC
#' table, the PCA, the VCF, the PDF plots and the interactive HTML report.
#'
#' \strong{Order of the QC steps.} (1) the chromosomes of \code{modechrom} are
#' selected; (2) sample outliers and problematic targets are decided on the RAW
#' counts of every candidate sample (test and reference); (3) the remaining
#' counts are GC-corrected; (4) each sample is called against its own references.
#' The QC table and the report describe exactly those decisions.
#'
#' @param gc_file            Character or NULL. Path to pre-computed GC table.
#' @param fasta_file         Character or NULL. Path to genome FASTA (used with bed_file).
#' @param bsgenome_pkg      Character or NULL. Name of a BSgenome package (e.g. "BSgenome.Hsapiens.UCSC.hg38") – used if \code{fasta_file} is not provided.
#' @param reads_file         Character or NULL. Path to pre-computed counts table
#'   (expected layout: \code{chromosome, start, end, GENE, <sample columns...>}).
#' @param modechrom          Character. "A" (autosomes), "XX", or "XY".
#' @param removeY            Logical. If TRUE (default), exclude chrY targets (so
#'   \code{modechrom = "XY"} analyses chrX only unless \code{removeY = FALSE}).
#' @param samples            Character. Directory of BAMs or path to BAM list file.
#'   May be \code{NULL} when \code{reads_file} is given (the sample names are then
#'   the column names of that file).
#' @param p_value            Numeric. Prior CNV probability per target.
#' @param Tnum               Numeric. Expected number of targets per CNV.
#' @param D                  Numeric. Expected genomic span of a CNV (bp). Default 100 kb.
#' @param numrefs            Integer. Max reference samples.
#' @param min_cor            Numeric or NULL. Minimum Spearman correlation for references.
#' @param homdel_mean        Numeric. Expected read count under homozygous deletion.
#' @param auto_reference     Logical. Auto-infer reference sets per sample.
#' @param reference_noise_quantile Numeric (0-1). With \code{auto_reference}, the
#'   noisiest samples (above this quantile of the per-sample MAD) are not used
#'   as references. Default 0.90.
#' @param sample_qc          Logical. Exclude outlier samples before calling.
#' @param exon_qc            Logical. Exclude problematic targets before calling
#'   (high cross-sample MAD, low mean coverage, extreme or unknown GC).
#' @param qc_zscore          Numeric. Robust z-score threshold for sample QC.
#' @param outlier_sides      \code{"upper"} (default) excludes only samples that
#'   are unusually \emph{noisy}; \code{"both"} also excludes unusually clean ones
#'   (the historical ECHO behaviour). Same option and default in ECHO.
#' @param exon_mad_quantile  Numeric. MAD quantile threshold for target QC.
#' @param gc_extreme_filter  Numeric[2]. Min/max GC (fraction, 0-1) to retain a
#'   target (applied when \code{exon_qc = TRUE}).
#' @param min_exon_mean      Numeric. Minimum mean count to retain a target.
#' @param min_sample_count   Integer. Targets where the sample being called has
#'   fewer reads than this are excluded before HMM calling (default 5). This also
#'   hides homozygous deletions; runs of such targets are listed in the report.
#'   See \code{\link{call_cnvs}}.
#' @param genome_build       Optional genome build label for the VCF header
#'   (e.g. \code{"hg38"}). If \code{NULL} it is inferred from
#'   \code{bsgenome_pkg} / \code{fasta_file}.
#' @param output_file        Character. Path for CSV output.
#' @param rdata_output       Character. Path for RData workspace output.
#' @param output_prefix      Character. Base prefix for all output files (e.g. "CANOPE_TEST").
#'                           If NULL, derived from output_file (strips "_CNVCall.csv").
#' @param refbams_file       Character or NULL. Optional reference BAM list.
#' @param ref_reads          Character or NULL. Pre-computed reference counts table.
#' @param bed_file           Character or NULL. BED file for BAM coverage extraction
#'   (also used for GC computation when \code{fasta_file} is supplied).
#' @param coverage_backend   \code{"bioconductor"} (default) or \code{"megadepth"}.
#'   Only matters when extracting coverage directly from \code{samples}' BAMs
#'   (i.e. \code{reads_file} is NULL). \code{"megadepth"} shells out to the
#'   compiled megadepth tool (much faster; reports base-level coverage rather than
#'   a fragment count, see \code{\link{get_coverage_from_bams_megadepth}}).
#' @param megadepth_op       \code{"sum"} (default) or \code{"mean"}.
#' @param megadepth_threads  Integer. Threads per megadepth invocation.
#' @param decode_method      \code{"distance"} (default) or \code{"stationary"} HMM transitions.
#'   Ignored under \code{engine = "legacy_canoes"}.
#' @param engine              \code{"new"} (default) or \code{"legacy_canoes"}
#'   (a faithful port of the original CANOES HMM; see \code{\link{call_cnvs}}).
#' @param bed_process        \code{"NO"} (default), \code{"STANDARD"}, or \code{"REGEN"}.
#'   If not \code{"NO"}, \code{bed_file} is first run through
#'   \code{\link{process_bed_file}} (see \code{bed_process_args}).
#' @param bed_process_args   Named list of extra arguments forwarded to
#'   \code{\link{process_bed_file}} (e.g. \code{list(genome_version = "hg38")},
#'   \code{off_target_pattern}, \code{off_target_handling}).
#' @param pad_terminal_exons Integer >= 0. Bases to extend the outward-facing
#'   edge of each gene's first and last exon (clamped, never overlapping a
#'   neighbour). Only takes effect when extracting fresh from BAMs.
#' @param plot_gene_gap      Numeric >= 0. Extra x-axis units inserted between
#'   genes in every PDF and HTML-report panel. Default \code{1}.
#' @param remove_terminal_only,penalize_terminal_only Logical. Calls that only
#'   touch the first and/or last exon of a gene are often boundary artefacts:
#'   remove them, or force them to LOW confidence. Same options as ECHO.
#' @param recurrence_action  What to do with calls that recur in many samples of
#'   the cohort (a noisy target rather than a germline CNV): \code{"flag"}
#'   (default; adds \code{Carriers}/\code{Recurrent} and shows them in the
#'   report), \code{"downgrade"} (also LOW confidence), \code{"remove"} or
#'   \code{"off"}. Same options as ECHO.
#' @param recurrence_min_samples,recurrence_min_fraction,recurrence_min_overlap
#'   A call is recurrent when at least \code{recurrence_min_samples} samples AND
#'   at least \code{recurrence_min_fraction} of the analysed samples carry a
#'   same-type call overlapping it by \code{recurrence_min_overlap} (reciprocal).
#' @param run_qc_metrics     Logical. Compute and write the QC metrics table.
#' @param qc_output_file     Path for the QC metrics TSV. Default derived from
#'   \code{output_file}'s directory.
#' @param qc_min_corr,qc_min_cov,qc_min_total_reads,qc_max_exon_cv
#'   Thresholds of the QC metrics table (same names and defaults as ECHO's
#'   \code{min_corr}, \code{min_cov}, \code{min_total_reads}, \code{max_exon_cv}).
#' @param score_confidence   Logical. Annotate \code{Confidence}/\code{CN_label}
#'   via \code{\link{score_canope_confidence}}. Default TRUE.
#' @param confidence_args    Named list of extra arguments forwarded to
#'   \code{\link{score_canope_confidence}}.
#' @param generate_plots_flag Logical. Generate per-call static PDF plots. Default TRUE.
#' @param pca_plot           Logical. Generate a PCA-of-coverage PDF. Default TRUE.
#' @param pca_output_file    Path for the PCA PDF. Default derived from \code{output_file}.
#' @param sample_table       Optional path to a TSV with \code{sample_name} (and
#'   optionally \code{gender}) columns, used to colour the PCA plot and check the sex.
#' @param export_vcf         Logical. Export CNV calls to VCF. Default TRUE.
#' @param vcf_output         Path for the combined multi-sample VCF.
#' @param vcf_per_sample     Logical. Also write one VCF per sample. Default FALSE.
#' @param report             Logical. Render the interactive HTML report
#'   (\code{\link{generate_canope_report}}). Default TRUE.
#' @param report_output_dir  Directory for the rendered HTML.
#' @param log_file           Path of the pipeline log (also shown in the HTML
#'   report). Truncated at the start of every run.
#' @param log_level          Minimum level written to the log and the console:
#'   \code{"DEBUG"}, \code{"INFO"} (default), \code{"WARNING"} or \code{"ERROR"}.
#'
#' @return Invisibly, the data frame of all called CNVs.
#' @export
run_canope <- function(
    gc_file = NULL, fasta_file = NULL, bsgenome_pkg = NULL, reads_file = NULL,
    modechrom = "A", removeY = TRUE, samples = NULL,
    p_value = 1e-08, Tnum = 6, D = 100000,
    numrefs = 30, min_cor = NULL, homdel_mean = 0.2,
    auto_reference = TRUE, reference_noise_quantile = 0.90,
    sample_qc = TRUE, exon_qc = TRUE,
    qc_zscore = 3, outlier_sides = c("upper", "both"), exon_mad_quantile = 0.90,
    gc_extreme_filter = c(0.15, 0.85), min_exon_mean = 20,
    min_sample_count = 5L, genome_build = NULL,
    output_file = "CNVCall.csv", rdata_output = "canope_workspace.RData",
    output_prefix = NULL,
    refbams_file = NULL, ref_reads = NULL, bed_file = NULL,
    coverage_backend = c("bioconductor", "megadepth"),
    megadepth_op = c("sum", "mean"), megadepth_threads = 1L,
    decode_method = c("distance", "stationary"),
    engine = c("new", "legacy_canoes"),
    bed_process = "NO", bed_process_args = list(),
    pad_terminal_exons = 0, plot_gene_gap = 1,
    remove_terminal_only = FALSE, penalize_terminal_only = FALSE,
    recurrence_action = c("flag", "downgrade", "remove", "off"),
    recurrence_min_samples = 3L, recurrence_min_fraction = 0.25, recurrence_min_overlap = 0.5,
    run_qc_metrics = TRUE, qc_output_file = NULL,
    qc_min_corr = 0.98, qc_min_cov = 100, qc_min_total_reads = 300000, qc_max_exon_cv = 0.5,
    score_confidence = TRUE, confidence_args = list(),
    generate_plots_flag = TRUE,
    pca_plot = TRUE, pca_output_file = NULL, sample_table = NULL,
    export_vcf = TRUE, vcf_output = NULL, vcf_per_sample = FALSE,
    report = TRUE, report_output_dir = NULL,
    log_file = NULL, log_level = "INFO"
) {
  decode_method     <- match.arg(decode_method)
  engine            <- match.arg(engine)
  coverage_backend  <- match.arg(coverage_backend)
  megadepth_op      <- match.arg(megadepth_op)
  outlier_sides     <- match.arg(outlier_sides)
  recurrence_action <- match.arg(recurrence_action)
  modechrom         <- toupper(as.character(modechrom)[1])
  if (!modechrom %in% c("A", "XX", "XY"))
    stop("[ERROR] modechrom must be 'A', 'XX' or 'XY' (got '", modechrom, "').", call. = FALSE)

  if (is.null(output_prefix)) {
    base <- basename(output_file)
    output_prefix <- if (grepl("_CNVCall\\.csv$", base)) sub("_CNVCall\\.csv$", "", base)
                     else tools::file_path_sans_ext(base)
  }
  if (!identical(output_prefix, "CANOPE") && !grepl("^CANOPE_", output_prefix))
    output_prefix <- paste0("CANOPE_", output_prefix)
  out_dir <- dirname(output_file)
  if (out_dir == "") out_dir <- "."
  if (is.null(log_file)) log_file <- file.path(out_dir, paste0(output_prefix, "_pipeline.log"))
  for (path in c(output_file, rdata_output))
    if (!is.null(path)) dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  # ONE logger for the whole run (shared with ECHO): timestamped lines in the
  # file and on the console, the file is truncated at the start so the report
  # only shows this run, and lg$capture() routes everything the callees emit
  # with message() / warning() into the same file.
  lg <- cnv_logger(log_file, tag = "CANOPE", level = log_level, total_steps = NA_integer_)
  lg$header(packages = c("nnls", "data.table", "GenomicRanges"))
  ok <- FALSE
  on.exit(lg$footer(ok), add = TRUE)
  run_step <- function(name, expr, on_error = "stop") lg$step(name, expr, on_error = on_error)

  final_cnvs <- NULL
  lg$capture({
    lg$set_total(4L + as.integer(isTRUE(run_qc_metrics)) + as.integer(isTRUE(pca_plot)) +
                   as.integer(isTRUE(generate_plots_flag)) + as.integer(isTRUE(export_vcf)) +
                   as.integer(isTRUE(report)))
    lg$info(sprintf("modechrom = '%s', removeY = %s | engine = '%s' | outlier_sides = '%s' | recurrence = '%s'",
                    modechrom, removeY, engine, outlier_sides, recurrence_action))
    lg$info(sprintf("sample_qc = %s, exon_qc = %s, qc_zscore = %s | terminal-exon padding: %s bp | remove / penalize terminal-only: %s / %s",
                    sample_qc, exon_qc, qc_zscore, pad_terminal_exons %||% 0,
                    remove_terminal_only, penalize_terminal_only))

    run_step("Preparing the BED file and the sample list", {
      if (!is.null(bed_file) && !identical(bed_process, "NO")) {
        processed_bed <- file.path(out_dir, paste0(output_prefix, "_processed.bed"))
        lg$info(sprintf("process_bed_file(mode = '%s') on %s", bed_process, basename(bed_file)))
        do.call(process_bed_file, c(list(input_bed = bed_file, output_bed = processed_bed,
                                         bed_process = bed_process), bed_process_args))
        bed_file <- processed_bed
      }
      bam_info <- if (!is.null(samples)) canope_resolve_bams(samples) else list(bams = character(0), names = NULL)
      if (is.null(samples) && (is.null(reads_file) || !file.exists(reads_file)))
        stop("[ERROR] 'samples' (BAM directory / list) is required when no reads_file is given.", call. = FALSE)
      sample_table_df <- NULL
      if (!is.null(sample_table) && file.exists(sample_table)) {
        sample_table_df <- utils::read.table(sample_table, header = TRUE, sep = "\t",
                                             stringsAsFactors = FALSE, check.names = FALSE)
        if (!"sample_name" %in% colnames(sample_table_df)) sample_table_df <- NULL
      }
    })

    run_step("Reading counts and GC content", {
      inp <- canope_read_input_counts(reads_file, bed_file, bam_info$bams, bam_info$names,
                                      pad_terminal_exons, fasta_file, bsgenome_pkg, coverage_backend,
                                      megadepth_op, megadepth_threads, out_dir, output_prefix, lg)
      bed_file <- inp$bed_file
      samples_to_analyse <- inp$sample_names
      reads <- inp$reads
      for (col in samples_to_analyse) {
        raw_vals <- reads[[col]]
        reads[[col]] <- suppressWarnings(as.numeric(as.character(raw_vals)))
        lost <- is.na(reads[[col]]) & !is.na(raw_vals) & nzchar(as.character(raw_vals))
        if (any(lost)) lg$warn(sprintf("Column '%s' had %d non-numeric value(s) coerced to NA", col, sum(lost)))
      }
      reads$gc <- canope_target_gc(reads, fasta_file, bsgenome_pkg, gc_file, bed_file, lg)
      n_gc_missing <- sum(is.na(reads$gc))
      if (n_gc_missing > 0)
        lg$warn(sprintf("GC content is unknown for %d of %d target(s) (interval not found in the reference sequence)%s.",
                        n_gc_missing, nrow(reads),
                        if (exon_qc) ": they are removed by exon_qc" else ": kept because exon_qc = FALSE"))
      panel <- canope_attach_reference_panel(reads, refbams_file, ref_reads, modechrom, lg)
      reads <- panel$reads
      refsample_names <- panel$refsample_names
      candidates <- c(samples_to_analyse, refsample_names)
      all_sample_names <- samples_to_analyse

      # What the BED / count table holds BEFORE modechrom decides what is
      # analysed: the report states which chromosome classes were left out, and
      # gives a per-chromosome depth / sex check.
      input_chromosomes <- as.character(reads$chromosome)
      chrom_levels <- unique(input_chromosomes)
      chrom_depth <- tryCatch(cnv_chrom_depth(reads[, candidates, drop = FALSE], input_chromosomes), error = function(e) NULL)
    })

    run_step("Selecting targets and removing outlier samples / problematic targets", {
      reads <- canope_prepare_targets(reads, modechrom, removeY, lg)
      n_targets_mode <- nrow(reads)
      raw_mat <- as.matrix(reads[, candidates, drop = FALSE]); storage.mode(raw_mat) <- "double"
      # one row per analysed-chromosome target BEFORE QC; the row number is the
      # target id, so the terminal-exon rule sees the gene's real exon set even
      # when QC removed one of its terminal exons
      raw_bed <- reads[, c("target", "chromosome", "start", "end", "GENE", "exon_number")]

      # ONE QC decision for both pipelines and for the QC-metrics step (see
      # cnv_apply_qc()): outlier samples among ALL candidate samples first, then
      # the problematic targets among the samples that remain -- on RAW counts.
      qc <- cnv_apply_qc(raw_mat, chromosomes = reads$chromosome, gc = reads$gc,
                         sample_qc = sample_qc, exon_qc = exon_qc, qc_zscore = qc_zscore,
                         outlier_sides = outlier_sides, exon_mad_quantile = exon_mad_quantile,
                         min_exon_mean = min_exon_mean, gc_range = gc_extreme_filter)
      excluded_samples <- qc$excluded_samples
      if (nrow(excluded_samples) > 0) {
        lg$warn(sprintf("sample_qc: excluding %d outlier sample(s) from calling (used neither as test nor reference): %s",
                        nrow(excluded_samples),
                        paste(sprintf("%s (z = %.2f)", excluded_samples$sample,
                                      qc$sample_qc$robust_z[match(excluded_samples$sample, qc$sample_qc$sample)]),
                              collapse = ", ")))
      }
      samples_to_analyse <- intersect(samples_to_analyse, qc$kept_samples)
      refsample_names    <- intersect(refsample_names, qc$kept_samples)
      candidates         <- intersect(candidates, qc$kept_samples)
      if (!length(samples_to_analyse))
        stop("[ERROR] sample_qc excluded every sample; relax qc_zscore or set sample_qc = FALSE.", call. = FALSE)
      if (length(candidates) < 4L)
        lg$warn(sprintf("Only %d sample(s) left for calling; CANOPE needs >= 3 reference samples per sample, so calls are unlikely or unreliable.",
                        length(candidates)))

      removed_targets <- NULL
      if (isTRUE(exon_qc) && length(qc$keep_targets) < n_targets_mode) {
        drop_idx <- setdiff(seq_len(n_targets_mode), qc$keep_targets)
        removed_targets <- data.frame(
          chromosome = raw_bed$chromosome[drop_idx], start = raw_bed$start[drop_idx], end = raw_bed$end[drop_idx],
          gene = raw_bed$GENE[drop_idx], label = target_labels(raw_bed)[drop_idx],
          reason = qc$target_qc$reason[drop_idx], stringsAsFactors = FALSE)
        lg$info(sprintf("exon_qc: excluding %d/%d problematic target(s) before calling (%d low mean, %d high variability, %d extreme/unknown GC)",
                        length(drop_idx), n_targets_mode, sum(qc$target_qc$low_mean[drop_idx]),
                        sum(qc$target_qc$high_mad[drop_idx]),
                        sum(qc$target_qc$gc_extreme[drop_idx] | grepl("unavailable", qc$target_qc$reason[drop_idx]))))
      }

      # GC correction on the retained raw counts, then what the HMM sees
      keep <- qc$keep_targets
      reads <- reads[keep, , drop = FALSE]
      rownames(reads) <- NULL
      corrected <- gc_correct_counts(raw_mat[keep, candidates, drop = FALSE], gc = reads$gc, method = "loess")
      corrected <- round(pmax(corrected, 0))
      counts_tab <- cbind(reads[, c("target", "chromosome", "start", "end", "gc", "GENE")],
                          as.data.frame(corrected, check.names = FALSE))
      bed_tab <- reads[, c("chromosome", "start", "end", "target", "gc", "GENE", "exon_number")]
      lg$info(sprintf("%d of %d target(s) and %d sample(s) go to the caller",
                      nrow(counts_tab), n_targets_mode, length(candidates)))
    })

    run_step("Calling CNVs", {
      auto_refs <- NULL
      if (auto_reference) {
        ref_obj <- infer_reference_samples(corrected, top_n = numrefs, contamination_quantile = reference_noise_quantile)
        auto_refs <- ref_obj$refs
      }
      par <- list(p_value = p_value, Tnum = Tnum, D = D, numrefs = numrefs, homdel_mean = homdel_mean,
                  min_cor = min_cor, decode_method = decode_method, engine = engine,
                  min_sample_count = min_sample_count)
      res <- canope_call_all(counts_tab, samples_to_analyse, refsample_names, auto_refs, auto_reference,
                             par, excluded_samples, lg)
      scored <- isTRUE(score_confidence)
      pooled <- if (length(res$calls) > 0) do.call(rbind, res$calls) else .empty_cnv_calls(scored = FALSE)
      rownames(pooled) <- NULL
      final_cnvs <- canope_finalise_calls(
        pooled, raw_bed, n_samples = length(samples_to_analyse), scored = scored,
        confidence_args = confidence_args, remove_terminal_only = remove_terminal_only,
        penalize_terminal_only = penalize_terminal_only, recurrence_action = recurrence_action,
        recurrence_min_samples = recurrence_min_samples, recurrence_min_fraction = recurrence_min_fraction,
        recurrence_min_overlap = recurrence_min_overlap, lg = lg,
        models = res$models, status = res$status)
      failed <- vapply(res$status, function(r) identical(r$status, "failed"), logical(1))
      if (any(failed))
        lg$warn(sprintf("%d sample(s) FAILED and have no calls (this is not the same as 'no CNV found'): %s",
                        sum(failed), paste(names(res$status)[failed], collapse = ", ")))
      # The table is always written (header only if nothing was called), so
      # downstream steps and the HTML report never have to guess.
      csv_cols <- setdiff(colnames(final_cnvs), CANOPE_CSV_EXCLUDE)
      utils::write.csv(final_cnvs[, csv_cols, drop = FALSE], output_file, row.names = FALSE)
      lg$info(sprintf("%d CNV(s) called across %d sample(s); table written to %s",
                      nrow(final_cnvs), length(samples_to_analyse), basename(output_file)))
    })

    if (is.null(genome_build))
      genome_build <- canope_infer_genome_build(bsgenome_pkg, fasta_file, bed_process_args$genome_version)
    sample_status <- if (length(res$status) > 0) do.call(rbind, res$status) else data.frame()
    rownames(sample_status) <- NULL
    sex_check <- if (!is.null(chrom_depth))
      cnv_sex_check(chrom_depth, sample_table_df, expected_sex = switch(modechrom, XX = "F", XY = "M", NA_character_),
                    analysed = samples_to_analyse, mode = modechrom, lg = lg) else NULL
    scope <- canope_scope(input_chromosomes, modechrom, removeY, analysed_chromosomes = counts_tab$chromosome)
    # the QC-metrics step logs the same messages; log them here only when it is skipped
    if (!run_qc_metrics) {
      scope_msgs <- chromosome_scope_warnings(scope, "modechrom", CANOPE_SCOPE_HINTS)
      for (k in seq_len(nrow(scope_msgs)))
        lg$log(scope_msgs$Details[k], if (scope_msgs$Severity[k] == "info") "INFO" else "WARNING")
    }

    qc_info <- list(
      params = list(sample_qc = sample_qc, exon_qc = exon_qc, qc_zscore = qc_zscore,
                    outlier_sides = outlier_sides, exon_mad_quantile = exon_mad_quantile,
                    min_exon_mean = min_exon_mean, gc_extreme_filter = gc_extreme_filter,
                    reference_noise_quantile = reference_noise_quantile),
      sample_qc = qc$sample_qc, excluded_samples = excluded_samples, sample_metrics = qc$sample_metrics,
      removed_targets = removed_targets,
      n_targets = list(input = length(input_chromosomes), mode = n_targets_mode, analysed = nrow(counts_tab),
                       gc_missing = n_gc_missing),
      candidate_samples = c(all_sample_names, refsample_names),
      recurrence = list(action = recurrence_action, min_samples = recurrence_min_samples,
                        min_fraction = recurrence_min_fraction, min_overlap = recurrence_min_overlap),
      terminal = list(remove = remove_terminal_only, penalize = penalize_terminal_only))

    run_info <- list(
      canope_version = tryCatch(as.character(utils::packageVersion("CANOPE")), error = function(e) "dev"),
      generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      prefix = output_prefix, modechrom = modechrom, removeY = removeY, engine = engine,
      decode_method = decode_method, coverage_backend = coverage_backend, genome_build = genome_build,
      params = list(
        p_value = p_value, Tnum = Tnum, D = D, numrefs = numrefs, min_cor = min_cor,
        homdel_mean = homdel_mean, min_sample_count = min_sample_count,
        auto_reference = auto_reference, reference_noise_quantile = reference_noise_quantile,
        sample_qc = sample_qc, exon_qc = exon_qc, qc_zscore = qc_zscore, outlier_sides = outlier_sides,
        exon_mad_quantile = exon_mad_quantile, gc_extreme_filter = gc_extreme_filter,
        min_exon_mean = min_exon_mean, pad_terminal_exons = pad_terminal_exons, bed_process = bed_process,
        remove_terminal_only = remove_terminal_only, penalize_terminal_only = penalize_terminal_only,
        recurrence_action = recurrence_action, recurrence_min_samples = recurrence_min_samples,
        recurrence_min_fraction = recurrence_min_fraction, recurrence_min_overlap = recurrence_min_overlap),
      thresholds = list(qc_min_corr = qc_min_corr, qc_min_cov = qc_min_cov,
                        qc_min_total_reads = qc_min_total_reads, qc_max_exon_cv = qc_max_exon_cv),
      chrom_depth = chrom_depth, sex_check = sex_check,
      samples_requested = all_sample_names, samples_analysed = samples_to_analyse,
      sample_status = sample_status,
      low_coverage_runs = if (length(res$low_coverage_runs) > 0) do.call(rbind, res$low_coverage_runs) else NULL,
      sample_table = sample_table_df)

    if (!is.null(rdata_output)) {
      cnv_calls <- final_cnvs
      counts <- counts_tab; bed_file <- bed_tab
      models <- res$models; refs <- res$refs
      raw_counts <- raw_mat; hmm_engine_used <- engine
      sample_names <- samples_to_analyse; all_samples <- candidates
      save(cnv_calls, counts, bed_file, models, refs, hmm_engine_used, run_info, qc_info,
           raw_counts, raw_bed, input_chromosomes, modechrom, removeY, sample_names, all_samples,
           file = rdata_output)
      lg$info("Workspace written to ", basename(rdata_output))
    }

    plot_dir <- file.path(out_dir, "plots")
    qc_metrics_path <- NULL
    if (run_qc_metrics) {
      if (is.null(qc_output_file)) qc_output_file <- file.path(out_dir, paste0(output_prefix, "_QC_metrics.tsv"))
      run_step("Computing QC metrics", {
        run_canope_qc_metrics(
          counts = raw_mat, sample_names = colnames(raw_mat), bed_df = raw_bed,
          output_file = qc_output_file, min_corr = qc_min_corr, min_cov = qc_min_cov,
          min_total_reads = qc_min_total_reads, max_exon_cv = qc_max_exon_cv,
          all_chromosomes = input_chromosomes, modechrom = modechrom, removeY = removeY,
          qc = qc, removed_targets = removed_targets, qc_zscore = qc_zscore, sample_qc = sample_qc)
        qc_metrics_path <- qc_output_file
      }, on_error = "warn")
    }

    if (pca_plot && ncol(raw_mat) >= 3) {
      if (is.null(pca_output_file)) pca_output_file <- file.path(out_dir, paste0(output_prefix, "_PCA.pdf"))
      run_step("Generating the coverage PCA plot", {
        color_by <- if (!is.null(sample_table_df) && "gender" %in% colnames(sample_table_df))
          sample_table_df$gender[match(colnames(raw_mat), sample_table_df$sample_name)] else NULL
        plot_coverage_pca(raw_mat, colnames(raw_mat), output_pdf = pca_output_file, color_by = color_by)
      }, on_error = "warn")
    }

    if (generate_plots_flag && !is.null(rdata_output)) {
      run_step("Generating the static CNV plots", {
        generate_plots(rdata_file = rdata_output, output_dir = plot_dir, modechrom = modechrom,
                       prefix = output_prefix, gene_gap = plot_gene_gap, log_file = log_file)
      }, on_error = "warn")
    }

    if (export_vcf) {
      if (is.null(vcf_output)) vcf_output <- file.path(out_dir, paste0(output_prefix, "_calls.vcf"))
      run_step("Exporting the calls to VCF", {
        export_canope_to_vcf(final_cnvs, vcf_output, reference = genome_build, sample_names = samples_to_analyse)
        if (vcf_per_sample)
          for (samp in samples_to_analyse)
            export_canope_to_vcf(final_cnvs, file.path(out_dir, paste0(sanitize_filename(samp), ".vcf")),
                                 sample_name = samp, reference = genome_build)
      }, on_error = "warn")
    }

    # The report is rendered even when nothing was called: that is exactly when
    # the QC metrics and the warnings (failed / excluded samples, chromosomes not
    # analysed) matter most.
    if (report) {
      if (is.null(rdata_output)) {
        lg$warn("report = TRUE needs rdata_output; skipping the HTML report.")
      } else {
        if (is.null(report_output_dir)) report_output_dir <- out_dir
        run_step("Rendering the interactive HTML report", {
          generate_canope_report(
            rdata_output = rdata_output, qc_metrics_file = qc_metrics_path, output_dir = report_output_dir,
            prefix = sub("^CANOPE_?", "", output_prefix),
            # user-supplied confidence_args override the defaults shown in the report
            settings = utils::modifyList(
              c(lapply(formals(score_canope_confidence)[c("high_q_score", "med_q_score", "high_num_refs",
                                                          "med_num_refs", "high_num_targ", "med_num_targ",
                                                          "low_confidence_genes")], eval),
                list(qc_min_cov = qc_min_cov, qc_min_total_reads = qc_min_total_reads,
                     qc_min_corr = qc_min_corr, qc_max_exon_cv = qc_max_exon_cv,
                     plot_gene_gap = plot_gene_gap, remove_terminal_only = remove_terminal_only,
                     penalize_terminal_only = penalize_terminal_only)),
              confidence_args),
            sample_table = sample_table, log_file = log_file)
        }, on_error = "warn")
      }
    }
    ok <- TRUE
  })
  invisible(final_cnvs)
}
