#' Directory this file was \code{source()}d from (if it was)
#'
#' The quick-start usage in the README sources every pipeline file —
#' including this one — directly with plain filenames, which only works if
#' the current working directory happens to be the folder all those files
#' (and \code{CANOPE_report.Rmd}) live in. As soon as \code{run_canope()} is
#' called from a different working directory (a wrapper script, an RStudio
#' project rooted elsewhere, a scheduled job, etc.), \code{getwd()} no
#' longer points at that folder and the template can't be found even though
#' it's sitting right next to this very file. This captures that location
#' at source-time as a robust fallback, independent of the caller's
#' \code{getwd()}. Silently evaluates to \code{NULL} when the package is
#' loaded normally (installed and \code{library()}'d) rather than sourced,
#' since \code{system.file()} already covers that case.
#' @noRd
.canope_report_source_dir <- local({
  ofile <- tryCatch({
    of <- NULL
    for (fr in sys.frames()) {
      cand <- get0("ofile", envir = fr, inherits = FALSE)
      if (!is.null(cand)) of <- cand
    }
    of
  }, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile) && file.exists(ofile)) dirname(normalizePath(ofile)) else NULL
})


#' Generate Interactive HTML Report for CANOPE
#'
#' Mirrors ECHO's \code{generate_report()}: loads the CANOPE workspace RData,
#' optionally a QC-metrics TSV and a sample table, then renders
#' \code{CANOPE_report.Rmd} with those objects as parameters.
#'
#' @param rdata_output   Path to the CANOPE workspace \code{.RData}
#'   (as written by \code{\link{run_canope}}; must contain
#'   \code{cnv_calls}, \code{counts}, \code{bed_file}, \code{models}, \code{refs}).
#' @param qc_metrics_file Optional path to a QC-metrics TSV (from
#'   \code{\link{run_canope_qc_metrics}}). If \code{NULL} or missing, the
#'   report's QC section is skipped gracefully.
#' @param output_dir     Directory for the rendered HTML (created if missing).
#' @param settings       Optional named list of confidence-scoring thresholds
#'   (as passed to \code{\link{score_canope_confidence}}) — shown in the
#'   report's "Confidence Score Method" section.
#' @param sample_table   Optional path to a sample table (TSV with at least
#'   a \code{sample_name} column; \code{gender} optional) for PCA coloring.
#' @param log_file       Optional path to a pipeline log file to display.
#' @param prefix         Filename prefix for the output HTML
#'   (\code{CANOPE_<prefix>_report.html}). Default \code{"CANOPE"}.
#' @param template_path  Optional explicit path to \code{CANOPE_report.Rmd}.
#'   If \code{NULL}, looks for it next to the RData, then in the working
#'   directory, then via
#'   \code{system.file("rmd/CANOPE_report.Rmd", package = "CANOPE")} (the
#'   template's actual location, \code{inst/rmd/}, once installed).
#'
#' @return Invisibly returns the path to the rendered HTML file.
#' @export
generate_canope_report <- function(rdata_output,
                                   qc_metrics_file = NULL,
                                   output_dir,
                                   settings = NULL,
                                   sample_table = NULL,
                                   log_file = NULL,
                                   prefix = "CANOPE",
                                   template_path = NULL) {
  # Captured immediately, before anything else in this function could change
  # it, so relative paths the *caller* built (log_file, qc_metrics_file,
  # sample_table, rdata_output itself) keep meaning what the caller meant by
  # them — see the knit_root_dir note below.
  caller_wd <- getwd()

  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is required for report generation.")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("Failed to create output directory: ", output_dir)
  }

  local_env <- new.env()
  load(rdata_output, envir = local_env)
  required <- c("cnv_calls", "counts", "bed_file", "models", "refs")
  missing_vars <- setdiff(required, ls(local_env))
  if (length(missing_vars) > 0)
    stop("[ERROR] Missing variables in RData: ", paste(missing_vars, collapse = ", "))

  qc_metrics <- NULL
  if (!is.null(qc_metrics_file) && file.exists(qc_metrics_file)) {
    qc_metrics <- data.table::fread(qc_metrics_file, header = TRUE, fill = TRUE, data.table = FALSE)
  } else {
    message("[INFO] No QC metrics file supplied/found; report will skip the QC section.")
  }

  sample_df <- NULL
  if (!is.null(sample_table) && file.exists(sample_table)) {
    sample_df <- utils::read.table(sample_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (!"sample_name" %in% colnames(sample_df)) {
      warning("sample_table missing required column 'sample_name'; ignoring.")
      sample_df <- NULL
    }
  }

  if (is.null(template_path)) {
    local_candidate <- file.path(dirname(rdata_output), "CANOPE_report.Rmd")
    # BUG FIX: the template ships in inst/rmd/CANOPE_report.Rmd, so once the
    # package is installed it resolves via system.file() under "rmd/...", not
    # "rmarkdown/...". The old candidate never matched anything post-install
    # (system.file() silently returns "" for a nonexistent subdir), so this
    # fallback always failed unless a copy happened to sit next to the RData
    # or in the working directory — which is why the report wasn't found.
    pkg_candidate   <- system.file("rmd/CANOPE_report.Rmd", package = "CANOPE")
    here_candidate  <- file.path(getwd(), "CANOPE_report.Rmd")
    # BUG FIX: when these files are sourced directly (the README's quick-start
    # usage) rather than installed as a package, the only fallback that ever
    # matched was here_candidate — which depends on getwd() being exactly the
    # folder the scripts were sourced from. Call run_canope() from anywhere
    # else and the template silently couldn't be found, even sitting right
    # next to canope_report.R. Add that folder as a fallback independent of
    # getwd().
    source_candidate <- if (!is.null(.canope_report_source_dir))
      file.path(.canope_report_source_dir, "CANOPE_report.Rmd") else NULL
    candidates <- c(local_candidate, here_candidate, source_candidate, pkg_candidate)
    template_path <- Find(function(p) !is.null(p) && nzchar(p) && file.exists(p), candidates)
    if (is.null(template_path))
      stop("[ERROR] Could not locate CANOPE_report.Rmd. Looked in: ",
           paste(Filter(function(p) !is.null(p) && nzchar(p), candidates), collapse = "; "),
           ". Pass template_path explicitly.")
  }

  hmm_engine_used <- if ("hmm_engine_used" %in% ls(local_env)) local_env$hmm_engine_used else "new"

  # BUG FIX (this is why the report wasn't generating): rmarkdown::render()
  # does NOT resolve a relative `output_file` against the caller's working
  # directory — it resolves it against the directory containing the *input*
  # Rmd file. `template_path` here is typically an absolute path to wherever
  # CANOPE_report.Rmd happens to live (e.g. next to the pipeline's other .R
  # files), which is almost never the same place as `output_dir`. So a
  # perfectly valid, already-created `output_dir` like "results/report"
  # would make render() look for "results/report" *relative to the Rmd's own
  # folder* and fail with "The directory 'results/report' does not exist" —
  # even though that directory unquestionably exists relative to where
  # run_canope() was actually called from. Confirmed by direct reproduction:
  # calling rmarkdown::render() with a relative output_file from a working
  # directory other than the Rmd's own directory throws exactly that error,
  # regardless of whether the target directory exists relative to getwd().
  # Fix: build an absolute output_file. output_dir is already guaranteed to
  # exist at this point (created above), so normalizePath() on it is safe;
  # only the (possibly nonexistent) filename is appended afterward.
  out_file <- file.path(normalizePath(output_dir, mustWork = TRUE),
                        paste0("CANOPE_", prefix, "_report.html"))

  rmarkdown::render(
    template_path,
    params = list(
      cnv_calls    = local_env$cnv_calls,
      qc_metrics   = qc_metrics,
      counts       = local_env$counts,
      bed_file     = local_env$bed_file,
      models       = local_env$models,
      refs         = local_env$refs,
      settings     = settings,
      sample_table = sample_df,
      log_file     = log_file,
      hmm_engine_used = hmm_engine_used
    ),
    output_file = out_file,
    # BUG FIX: by default, rmarkdown::render() evaluates the Rmd's own code
    # chunks with the working directory set to the *Rmd's* directory, not
    # the caller's. That's a second, quieter instance of the same
    # relative-path problem fixed above for `out_file`: any relative path a
    # chunk touches (here, `log_file`, built by run_canope() relative to
    # *its* caller's working directory) would resolve against the wrong
    # base and silently fail its file.exists() check — the report would
    # still render, just missing its pipeline-log tab, with nothing telling
    # you why. Pinning knit_root_dir to the working directory this function
    # was actually called from keeps every relative path inside the Rmd
    # meaning what its caller meant by it.
    knit_root_dir = caller_wd,
    quiet = FALSE
  )
  message("[INFO] CANOPE report written: ", out_file)
  invisible(out_file)
}
