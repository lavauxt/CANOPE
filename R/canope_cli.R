#' Run CANOPE from a YAML Configuration File
#'
#' Reads a YAML configuration file and executes the full CANOPE CNV calling
#' pipeline with the specified parameters.  All arguments in \code{...} override
#' the corresponding values in the config file.
#'
#' @param config_path Path to a YAML configuration file. If the file does not
#'   exist, the function looks in the package's `inst/` directory (for a
#'   default config) and then in the current working directory.
#' @param ... Additional arguments passed to \code{\link{run_canope}}; these
#'   take precedence over values in the config file.
#'
#' @return Data frame of all called CNVs (invisibly).
#' @export
#' @examples
#' \dontrun{
#' canope("config.yaml")
#' canope("config.yaml", p_value = 1e-12, Tnum = 15)
#' }
canope <- function(config_path = "config.yaml", ...) {
  if (!requireNamespace("yaml", quietly = TRUE))
    stop("The 'yaml' package is required to read configuration files. ",
         "Install it with: install.packages('yaml')")

  config_file <- NULL
  if (file.exists(config_path)) {
    config_file <- config_path
  } else {
    pkg_config <- system.file("config.yaml", package = "CANOPE")
    if (file.exists(pkg_config)) {
      config_file <- pkg_config
    } else {
      dev_config <- file.path("inst", config_path)
      if (file.exists(dev_config)) {
        config_file <- dev_config
      }
    }
  }

  if (is.null(config_file))
    stop("Configuration file not found: ", config_path,
         "\nLooked in current directory, system.file('config.yaml', package='CANOPE'), and inst/.")

  cfg <- yaml::read_yaml(config_file)

  if (is.list(cfg$input)) {
    input   <- cfg$input   %||% list()
    output  <- cfg$output  %||% list()
    settings <- cfg$settings %||% list()
  } else {
    input   <- list()
    output  <- list()
    settings <- cfg
  }

  # Look a key up in settings / input / output / top level, in that order, so
  # both the nested layout (README) and the flat layout (inst/config.yaml) work
  # for every option -- e.g. `fasta_file` under `settings:` used to be ignored.
  pick <- function(...) {
    for (nm in c(...)) for (src in list(settings, input, output, cfg)) {
      if (is.list(src) && !is.null(src[[nm]])) return(src[[nm]])
    }
    NULL
  }

  args <- list()

  args$fasta_file <- pick("fasta", "fasta_file")
  args$bed_file   <- pick("bed", "bed_file")
  args$samples    <- pick("bamdir", "samples")
  args$bsgenome_pkg <- pick("bsgenome_pkg")

  # Optional inputs (pre-computed counts / GC / reference panel / sample sheet)
  for (key in c("reads_file", "gc_file", "refbams_file", "ref_reads", "sample_table",
                "genome_build", "vcf_output", "report_output_dir"))
    if (!is.null(v <- pick(key))) args[[key]] <- as.character(v)

  out_dir <- output$dir %||% cfg$output_dir %||% "."
  out_prefix <- output$prefix %||% cfg$output_prefix %||% "CANOPE"
  # "CANOPE_<prefix>"; the default prefix "CANOPE" must not become CANOPE_CANOPE
  base_prefix <- sub("^CANOPE_?", "", as.character(out_prefix))
  full_prefix <- if (nzchar(base_prefix)) paste0("CANOPE_", base_prefix) else "CANOPE"
  args$output_file   <- file.path(out_dir, paste0(full_prefix, "_CNVCall.csv"))
  args$rdata_output  <- file.path(out_dir, paste0(full_prefix, "_workspace.RData"))
  args$output_prefix <- full_prefix

  args$modechrom   <- settings$modechrom %||% cfg$modechrom %||% "A"
  args$min_cor     <- as.numeric(settings$min_cor %||% cfg$min_cor %||% 0.98)
  args$p_value     <- as.numeric(settings$p_value %||% cfg$p_value %||% 1e-08)
  args$Tnum        <- as.integer(settings$Tnum %||% cfg$Tnum %||% 6L)
  args$D           <- as.integer(settings$D %||% cfg$D %||% 100000L)
  args$numrefs     <- as.integer(settings$numrefs %||% cfg$numrefs %||% 30L)
  args$homdel_mean <- as.numeric(settings$homdel_mean %||% cfg$homdel_mean %||% 0.2)
  args$decode_method <- settings$decode_method %||% cfg$decode_method %||% "distance"
  args$engine      <- settings$engine %||% cfg$engine %||% "new"

  args$qc_min_corr         <- settings$qc_min_corr %||% cfg$qc_min_corr %||% 0.98
  args$qc_min_cov          <- settings$qc_min_cov %||% cfg$qc_min_cov %||% 100
  args$qc_min_total_reads  <- settings$qc_min_total_reads %||% cfg$qc_min_total_reads %||% 300000L
  args$qc_max_exon_cv      <- settings$qc_max_exon_cv %||% cfg$qc_max_exon_cv %||% 0.5

  # Extra run_canope() options that used to be settable only from R
  if (!is.null(v <- pick("removeY")))             args$removeY <- as.logical(v)
  if (!is.null(v <- pick("qc_zscore")))           args$qc_zscore <- as.numeric(v)
  if (!is.null(v <- pick("exon_mad_quantile")))   args$exon_mad_quantile <- as.numeric(v)
  if (!is.null(v <- pick("min_exon_mean")))       args$min_exon_mean <- as.numeric(v)
  if (!is.null(v <- pick("min_sample_count")))    args$min_sample_count <- as.integer(v)
  if (!is.null(v <- pick("gc_extreme_filter")))   args$gc_extreme_filter <- as.numeric(unlist(v))

  args$coverage_backend <- settings$coverage_backend %||% cfg$coverage_backend %||% "bioconductor"
  args$megadepth_op     <- settings$megadepth_op %||% cfg$megadepth_op %||% "sum"
  args$megadepth_threads <- settings$megadepth_threads %||% cfg$megadepth_threads %||% 1L


  args$bed_process      <- settings$bed_process %||% cfg$bed_process %||% "NO"
  if (is.list(settings$bed_process_args)) {
    args$bed_process_args <- settings$bed_process_args
  } else if (is.list(cfg$bed_process_args)) {
    args$bed_process_args <- cfg$bed_process_args
  } else {
    args$bed_process_args <- list()
  }

  args$pad_terminal_exons <- as.numeric(settings$pad_terminal_exons %||% cfg$pad_terminal_exons %||% 0)
  args$plot_gene_gap      <- as.numeric(settings$plot_gene_gap %||% cfg$plot_gene_gap %||% 1)

  args$run_qc_metrics      <- settings$run_qc_metrics %||% cfg$run_qc_metrics %||% TRUE
  args$score_confidence    <- settings$score_confidence %||% cfg$score_confidence %||% TRUE
  args$generate_plots_flag <- settings$generate_plots_flag %||% cfg$generate_plots_flag %||% TRUE
  args$pca_plot            <- settings$pca_plot %||% cfg$pca_plot %||% TRUE
  args$export_vcf          <- settings$export_vcf %||% cfg$export_vcf %||% TRUE
  args$vcf_per_sample      <- settings$vcf_per_sample %||% cfg$vcf_per_sample %||% FALSE
  args$report              <- settings$report %||% cfg$report %||% TRUE
  args$sample_qc           <- settings$sample_qc %||% cfg$sample_qc %||% TRUE
  args$exon_qc             <- settings$exon_qc %||% cfg$exon_qc %||% TRUE
  args$auto_reference      <- settings$auto_reference %||% cfg$auto_reference %||% TRUE


  if (!is.null(settings$confidence_args) && is.list(settings$confidence_args)) {
    args$confidence_args <- settings$confidence_args
  } else if (!is.null(cfg$confidence_args) && is.list(cfg$confidence_args)) {
    args$confidence_args <- cfg$confidence_args
  } else {
    args$confidence_args <- list()
  }


  user_args <- list(...)
  args <- modifyList(args, user_args)

  if (is.null(args$fasta_file) && is.null(args$bsgenome_pkg) && is.null(args$gc_file))
    stop("Either fasta_file or bsgenome_pkg (with bed_file), or a pre-computed gc_file, must be provided for GC computation.")
  if (is.null(args$samples))    stop("samples directory is required (input.bamdir in config)")
  if (is.null(args$bed_file))   stop("bed_file is required (input.bed in config)")

  message("[INFO] CANOPE pipeline starting with config: ", config_file)
  do.call(run_canope, args)
}
