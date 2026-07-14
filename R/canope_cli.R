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

  # Try to locate the config file
  config_file <- NULL
  if (file.exists(config_path)) {
    config_file <- config_path
  } else {
    # Look in the installed package's inst/ directory (if any)
    pkg_config <- system.file("config.yaml", package = "CANOPE")
    if (file.exists(pkg_config)) {
      config_file <- pkg_config
    } else {
      # Try inst/ relative to current working directory (development)
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

  # Flatten nested sections for easier access
  # Support both: input/output/settings structure, and flat structure
  if (is.list(cfg$input)) {
    # Nested structure (original ECHO-style)
    input   <- cfg$input   %||% list()
    output  <- cfg$output  %||% list()
    settings <- cfg$settings %||% list()
  } else {
    # Flat structure (simple key-value)
    input   <- list()
    output  <- list()
    settings <- cfg
  }

  # ------------------------------------------------------------------------
  # Map config values to run_canope() arguments
  # ------------------------------------------------------------------------
  args <- list()

  # Input files
  args$fasta_file <- input$fasta %||% cfg$fasta_file
  args$bed_file   <- input$bed   %||% cfg$bed_file
  args$samples    <- input$bamdir %||% cfg$samples

  # BSgenome package (optional, used for GC computation if fasta_file is not provided)
  args$bsgenome_pkg <- settings$bsgenome_pkg %||% cfg$bsgenome_pkg %||% NULL

  # Output files
  out_dir <- output$dir %||% cfg$output_dir %||% "."
  out_prefix <- output$prefix %||% cfg$output_prefix %||% "CANOPE"
  args$output_file   <- file.path(out_dir, paste0("CANOPE_", out_prefix, "_CNVCall.csv"))
  args$rdata_output  <- file.path(out_dir, paste0("CANOPE_", out_prefix, "_workspace.RData"))
  args$output_prefix <- paste0("CANOPE_", out_prefix)  # <-- pass the full prefix to run_canope

  # Settings
  args$modechrom   <- settings$modechrom %||% cfg$modechrom %||% "A"
  args$min_cor     <- as.numeric(settings$min_cor %||% cfg$min_cor %||% 0.98)
  args$p_value     <- as.numeric(settings$p_value %||% cfg$p_value %||% 1e-08)
  args$Tnum        <- as.integer(settings$Tnum %||% cfg$Tnum %||% 6L)
  args$D           <- as.integer(settings$D %||% cfg$D %||% 100000L)
  args$numrefs     <- as.integer(settings$numrefs %||% cfg$numrefs %||% 30L)
  args$homdel_mean <- as.numeric(settings$homdel_mean %||% cfg$homdel_mean %||% 0.2)
  args$decode_method <- settings$decode_method %||% cfg$decode_method %||% "distance"
  args$engine      <- settings$engine %||% cfg$engine %||% "new"

  # QC thresholds
  args$qc_min_corr         <- settings$qc_min_corr %||% cfg$qc_min_corr %||% 0.98
  args$qc_min_cov          <- settings$qc_min_cov %||% cfg$qc_min_cov %||% 100
  args$qc_min_total_reads  <- settings$qc_min_total_reads %||% cfg$qc_min_total_reads %||% 300000L
  args$qc_max_exon_cv      <- settings$qc_max_exon_cv %||% cfg$qc_max_exon_cv %||% 0.5

  # Coverage backend (megadepth / bioconductor)
  args$coverage_backend <- settings$coverage_backend %||% cfg$coverage_backend %||% "bioconductor"
  args$megadepth_op     <- settings$megadepth_op %||% cfg$megadepth_op %||% "sum"
  args$megadepth_threads <- settings$megadepth_threads %||% cfg$megadepth_threads %||% 1L

  # BED preprocessing
  args$bed_process      <- settings$bed_process %||% cfg$bed_process %||% "NO"
  if (is.list(settings$bed_process_args)) {
    args$bed_process_args <- settings$bed_process_args
  } else if (is.list(cfg$bed_process_args)) {
    args$bed_process_args <- cfg$bed_process_args
  } else {
    args$bed_process_args <- list()
  }

  # Terminal-exon padding / plot gene-gap spacing (ported from ECHO)
  args$pad_terminal_exons <- as.numeric(settings$pad_terminal_exons %||% cfg$pad_terminal_exons %||% 0)
  args$plot_gene_gap      <- as.numeric(settings$plot_gene_gap %||% cfg$plot_gene_gap %||% 1)

  # Flags
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

  # Confidence scoring thresholds (passed as a list to score_canope_confidence).
  # The documented config format is a single `confidence_args:` block (see
  # config.yaml / README) with keys matching score_canope_confidence()'s own
  # argument names directly (high_q_score, med_q_score, ...). If that block
  # isn't supplied, fall through to that function's own built-in defaults
  # by passing an empty list, rather than re-declaring the defaults here
  # under key names (score_high_q, score_med_q, ...) that no documented
  # config actually uses.
  if (!is.null(settings$confidence_args) && is.list(settings$confidence_args)) {
    args$confidence_args <- settings$confidence_args
  } else if (!is.null(cfg$confidence_args) && is.list(cfg$confidence_args)) {
    args$confidence_args <- cfg$confidence_args
  } else {
    args$confidence_args <- list()
  }

  # Override with any arguments passed directly to canope()
  user_args <- list(...)
  args <- modifyList(args, user_args)

  # Ensure required arguments are present
  if (is.null(args$fasta_file) && is.null(args$bsgenome_pkg)) 
    stop("Either fasta_file or bsgenome_pkg must be provided (with bed_file) for GC computation.")
  if (is.null(args$samples))    stop("samples directory is required (input.bamdir in config)")
  if (is.null(args$bed_file))   stop("bed_file is required (input.bed in config)")

  message("[INFO] CANOPE pipeline starting with config: ", config_file)
  do.call(run_canope, args)
}


#' Null-coalescing operator for use inside canope()
#' @noRd
`%||%` <- function(x, y) if (!is.null(x)) x else y