# End-to-end CANOPE run on a small simulated count table (no BAM, no FASTA needed):
# planted deletion + noisy sample -> calls, VCF, unified report payload, CLI config.

skip_if_not_installed("jsonlite")
skip_if_not_installed("yaml")
skip_if_not_installed("GenomicRanges")

sim_reads <- function(dir) {
  set.seed(11)
  n_genes <- 40; per <- 4
  tg <- do.call(rbind, lapply(seq_len(n_genes), function(g) {
    chr <- if (g <= 26) "chr1" else "chr2"
    s <- (if (g <= 26) g else g - 26) * 200000 + (seq_len(per) - 1) * 1500
    data.frame(chromosome = chr, start = s, end = s + 150 + sample(0:100, per, TRUE),
               GENE = sprintf("GENE%02d", g), stringsAsFactors = FALSE)
  }))
  n <- nrow(tg); gc <- pmin(0.8, pmax(0.25, rnorm(n, 0.5, 0.08)))
  mu <- exp(rnorm(n, 5.2, 0.35)) * (1 + 1.2 * (gc - 0.5))
  samples <- sprintf("S%02d", 1:14)
  m <- vapply(samples, function(s) rnbinom(n, mu = mu * exp(rnorm(1, 0, 0.1)), size = 60), numeric(n))
  del <- which(tg$GENE == "GENE13")
  m[del, "S03"] <- rnbinom(length(del), mu = mu[del] * 0.5, size = 60)
  m[, "S14"] <- rnbinom(n, mu = mu, size = 3)                       # a noisy outlier sample
  utils::write.table(cbind(tg, as.data.frame(m)), file.path(dir, "reads.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(tg[, 1:3], GC_CONTENT = gc), file.path(dir, "gc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(list(targets = tg, samples = samples))
}

test_that("run_canope calls the planted deletion, excludes the noisy sample and writes a unified report", {
  d <- tempfile("canope_e2e"); dir.create(d)
  sim <- sim_reads(d)
  res <- suppressMessages(run_canope(
    reads_file = file.path(d, "reads.tsv"), gc_file = file.path(d, "gc.tsv"), numrefs = 8,
    output_file = file.path(d, "T_CNVCall.csv"), rdata_output = file.path(d, "T.RData"), output_prefix = "T",
    report = TRUE, report_output_dir = d, generate_plots_flag = FALSE, pca_plot = FALSE,
    outlier_sides = "both", recurrence_action = "downgrade"))
  expect_true(all(c("SAMPLE", "CNV", "CHR", "INTERVAL", "Q_SOME", "Confidence") %in% colnames(res)))
  hit <- res[res$SAMPLE == "S03" & res$CNV == "DEL" & res$GENE == "GENE13", ]
  expect_gte(nrow(hit), 1L)
  expect_false("S14" %in% res$SAMPLE)                                   # excluded as a noisy outlier

  vcf <- list.files(d, "\\.vcf$", full.names = TRUE); expect_length(vcf, 1L)
  expect_true(any(grepl("^#CHROM", readLines(vcf))))

  lg <- list.files(d, "pipeline\\.log$", full.names = TRUE); expect_length(lg, 1L)
  lines <- readLines(lg)
  expect_true(any(grepl("^\\[INFO\\] \\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2} ", lines)))
  expect_true(any(grepl("outlier_sides = 'both'", lines, fixed = TRUE)))

  html <- list.files(d, "report\\.html$", full.names = TRUE); expect_gte(length(html), 1L)
  txt <- paste(readLines(html[1], warn = FALSE), collapse = "\n")
  js <- sub('(?s).*<script id="unified-data" type="application/json">(.*?)</script>.*', "\\1", txt, perl = TRUE)
  pl <- jsonlite::fromJSON(js, simplifyVector = FALSE)
  expect_setequal(names(pl), REPORT_PAYLOAD_KEYS)
  expect_equal(pl$meta$source, "CANOPE")
  expect_true(!is.null(pl$meta$palette$del))                            # the palette the JS paints with
})

test_that("canope() reads the outlier / recurrence / log options from the YAML config", {
  d <- tempfile("canope_cli"); dir.create(d)
  sim_reads(d)
  cfg <- list(input = list(reads_file = file.path(d, "reads.tsv"), gc_file = file.path(d, "gc.tsv")),
              settings = list(numrefs = 8, outlier_sides = "both", recurrence_action = "off", log_level = "DEBUG",
                              report = FALSE, generate_plots_flag = FALSE, pca_plot = FALSE),
              output = list(dir = file.path(d, "out"), prefix = "CLI"))
  yml <- file.path(d, "config.yaml"); yaml::write_yaml(cfg, yml)
  res <- suppressMessages(canope(yml))
  expect_s3_class(res, "data.frame")
  lg <- list.files(file.path(d, "out"), "pipeline\\.log$", full.names = TRUE, recursive = TRUE)
  expect_gte(length(lg), 1L)
  expect_true(any(grepl("outlier_sides = 'both'.*recurrence = 'off'", readLines(lg[1]))))
  expect_true(any(grepl("^\\[DEBUG\\] .* reference\\(s\\)", readLines(lg[1]))))   # log_level DEBUG reaches the log file
  bad <- cfg; bad$input <- NULL; yaml::write_yaml(bad, file.path(d, "bad.yaml"))
  expect_error(suppressMessages(canope(file.path(d, "bad.yaml"))), "required")
})

test_that("the targets layer's plain functions run the pipeline without {targets}", {
  d <- tempfile("canope_tg"); dir.create(d); sim_reads(d)
  cfg <- list(input = list(reads_file = file.path(d, "reads.tsv"), gc_file = file.path(d, "gc.tsv")),
              settings = list(numrefs = 8, report = FALSE, generate_plots_flag = FALSE, pca_plot = FALSE),
              output = list(dir = file.path(d, "out"), prefix = "TG"))
  yml <- file.path(d, "config.yaml"); yaml::write_yaml(cfg, yml)
  args <- canope_config_args(yml, p_value = 1e-10)
  expect_equal(args$p_value, 1e-10)                                   # '...' overrides the config
  expect_equal(canope_targets_bams(args), character(0))               # reads_file given: no BAM to track
  counts <- canope_targets_counts(args, character(0))
  expect_equal(counts, file.path(d, "reads.tsv"))
  files <- suppressMessages(canope_targets_run(args, counts))
  expect_equal(basename(files[1]), "CANOPE_TG_CNVCall.csv")
  expect_true(all(file.exists(files)))
  expect_gt(length(files), 2L)                                        # table + VCF + log + QC metrics ...
})

test_that("canope_targets() builds the DAG, or says what is missing", {
  if (requireNamespace("targets", quietly = TRUE)) {
    tl <- canope_targets("config.yaml")
    nm <- vapply(tl, function(t) t$settings$name, "")
    expect_equal(nm, c("canope_config", "canope_args", "canope_bams", "canope_counts", "canope_run", "canope_calls"))
  } else {
    expect_error(canope_targets("config.yaml"), "targets")
  }
})
