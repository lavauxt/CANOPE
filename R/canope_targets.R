# Opt-in {targets} layer for CANOPE.
#
# The expensive, parameter-independent step -- counting reads in the BAMs over
# the BED -- becomes a cached target of its own. Changing p_value, the outlier /
# recurrence options, the confidence rules ... then only re-runs the calling
# and the report; touching a BAM or the BED re-counts.
#
# All the logic sits in plain functions (canope_targets_*), so it can be used
# and tested without {targets}; canope_targets() only wires them into a DAG.

#' CANOPE as a {targets} pipeline (opt-in)
#'
#' Returns a list of targets for a \code{_targets.R} file:
#' \preformatted{
#' library(targets)
#' list(CANOPE::canope_targets("config.yaml"))
#' }
#' Targets: \code{canope_config} (the YAML, tracked), \code{canope_args} (parsed
#' options), \code{canope_bams} (BAM files, tracked by size / mtime),
#' \code{canope_counts} (counts table written once per BAM / BED / counting
#' setting), \code{canope_run} (everything \code{run_canope()} writes) and
#' \code{canope_calls} (the call table as a data frame).
#'
#' Requires the \pkg{targets} package (>= 1.3.0), a suggested dependency.
#'
#' @param config_path Path to the YAML configuration (see \code{\link{canope}}).
#' @param ... Options overriding the config, as in \code{\link{canope}}.
#' @return A list of \code{targets::tar_target()} objects.
#' @export
canope_targets <- function(config_path = "config.yaml", ...) {
  cnv_require_targets()
  dots <- list(...)
  list(
    targets::tar_target(canope_config, config_path, format = "file"),
    targets::tar_target(canope_args, canope_config_args(canope_config, !!!dots)),
    targets::tar_target(canope_bams, canope_targets_bams(canope_args), format = "file_fast"),
    targets::tar_target(canope_counts, canope_targets_counts(canope_args, canope_bams), format = "file"),
    targets::tar_target(canope_run, canope_targets_run(canope_args, canope_counts), format = "file"),
    targets::tar_target(canope_calls, utils::read.csv(canope_run[1L], check.names = FALSE, stringsAsFactors = FALSE))
  )
}

#' BAM files the counts depend on (empty when a pre-computed reads_file is used)
#' @noRd
canope_targets_bams <- function(args) {
  if (!is.null(args$reads_file) && file.exists(args$reads_file)) return(character(0))
  canope_resolve_bams(args$samples)$bams
}

#' The counts table of the run: the config's reads_file, or one counted from the BAMs
#'
#' The first element is the counts table (\code{chromosome, start, end, GENE},
#' one column per sample -- the \code{reads_file} format); when it was counted
#' here, the second is the BED that was actually used (processed / padded).
#' @noRd
canope_targets_counts <- function(args, bams = character(0)) {
  if (!is.null(args$reads_file) && file.exists(args$reads_file)) return(args$reads_file)
  out_dir <- dirname(args$output_file %||% file.path(getwd(), "CANOPE_out", "calls.csv"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- args$output_prefix %||% "CANOPE"
  lg <- cnv_logger(NULL, tag = "CANOPE", level = args$log_level %||% "INFO")
  bed <- args$bed_file
  if (!is.null(bed) && !identical(args$bed_process %||% "NO", "NO")) {
    processed <- file.path(out_dir, paste0(prefix, "_processed.bed"))
    do.call(process_bed_file, c(list(input_bed = bed, output_bed = processed, bed_process = args$bed_process),
                                args$bed_process_args))
    bed <- processed
  }
  info <- canope_resolve_bams(args$samples)
  inp <- canope_read_input_counts(NULL, bed, bams, info$names, args$pad_terminal_exons %||% 0,
                                  args$fasta_file, args$bsgenome_pkg,
                                  args$coverage_backend %||% "bioconductor", args$megadepth_op %||% "sum",
                                  args$megadepth_threads %||% 1L, out_dir, prefix, lg)
  counts_file <- file.path(out_dir, paste0(prefix, "_counts.tsv"))
  utils::write.table(inp$reads, counts_file, sep = "\t", quote = FALSE, row.names = FALSE)
  c(counts_file, inp$bed_file)
}

#' Run the calling on a counts table; returns the files the run wrote (call table first)
#' @noRd
canope_targets_run <- function(args, counts) {
  t0 <- Sys.time()
  args$reads_file <- counts[1L]
  if (length(counts) > 1L) {            # counted here: the BED is already processed / padded
    args$bed_file <- counts[2L]
    args$bed_process <- "NO"
    args$pad_terminal_exons <- 0
    args$samples <- NULL
  }
  attr(args, "config_file") <- NULL
  out_dir <- dirname(args$output_file %||% file.path(getwd(), "CANOPE_out", "calls.csv"))
  do.call(run_canope, args)
  files <- cnv_files_since(out_dir, t0)
  unique(c(args$output_file, files))
}
