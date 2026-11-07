skip_if_not_installed("jsonlite")

test_that("end-to-end run finds the spiked CNVs and writes VCF, QC table and workspace", {
  d <- withr::local_tempdir()
  make_synthetic_cohort(file.path(d, "in"))
  res <- run_synthetic(file.path(d, "in"), file.path(d, "out"))

  expect_s3_class(res, "data.frame")
  del <- res[res$SAMPLE == "S03" & res$CNV == "DEL", ]
  dup <- res[res$SAMPLE == "S07" & res$CNV == "DUP", ]
  expect_true(any(grepl("G2_2", del$GENE)))
  expect_true(any(grepl("G10_1", dup$GENE)))

  expect_true(file.exists(file.path(d, "out", "CANOPE_T_calls.vcf")))
  expect_true(file.exists(file.path(d, "out", "CANOPE_T_QC_metrics.tsv")))

  e <- new.env(); load(file.path(d, "out", "CANOPE_T_workspace.RData"), envir = e)
  expect_true("run_info" %in% ls(e))
})

test_that("calls do not depend on the row order of the target table", {
  d <- withr::local_tempdir()
  make_synthetic_cohort(file.path(d, "karyo"), order = "karyo")
  make_synthetic_cohort(file.path(d, "lex"), order = "lex")
  a <- run_synthetic(file.path(d, "karyo"), file.path(d, "oa"))
  b <- run_synthetic(file.path(d, "lex"), file.path(d, "ob"))
  key <- function(x) { x <- x[order(x$SAMPLE, x$INTERVAL), c("SAMPLE", "CNV", "INTERVAL", "Q_SOME", "MLCN")]; rownames(x) <- NULL; x }
  expect_equal(key(a), key(b))
})

test_that("an autosome-only run records and reports that gonosomes were not analysed", {
  d <- withr::local_tempdir()
  make_synthetic_cohort(file.path(d, "in"))
  run_synthetic(file.path(d, "in"), file.path(d, "out"), modechrom = "A")

  qc <- utils::read.delim(file.path(d, "out", "CANOPE_T_QC_metrics.tsv"), check.names = FALSE)
  expect_true(any(grepl("not analysed", unlist(qc), ignore.case = TRUE)))

  html <- generate_canope_report(file.path(d, "out", "CANOPE_T_workspace.RData"),
                                 qc_metrics_file = file.path(d, "out", "CANOPE_T_QC_metrics.tsv"),
                                 output_dir = file.path(d, "out"), prefix = "T",
                                 plotly_source = "cdn")
  expect_true(file.exists(html))
  txt <- paste(readLines(html, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(txt, 'id="canope-data"', fixed = TRUE)
  dat <- jsonlite::fromJSON(sub('(?s).*<script id="canope-data" type="application/json">(.*?)</script>.*', "\\1", txt, perl = TRUE),
                            simplifyVector = FALSE)
  w <- vapply(dat$warnings, function(x) x$text, "")
  expect_true(any(grepl("NOT both analysed", w)))
  expect_false(grepl("</script><script>alert", txt, fixed = TRUE))
})

test_that("the report handles several samples and a run without any call", {
  d <- withr::local_tempdir()
  make_synthetic_cohort(file.path(d, "in"), spike = FALSE)
  run_synthetic(file.path(d, "in"), file.path(d, "out"), report = TRUE, p_value = 1e-30)
  expect_true(file.exists(file.path(d, "out", "CANOPE_T_report.html")))
})
