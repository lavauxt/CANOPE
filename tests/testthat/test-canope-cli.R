# Tests for canope_cli.R's YAML config -> run_canope() argument mapping.
#
# canope()'s only job before it hands off to run_canope() is this mapping,
# so we intercept the call to run_canope() with local_mocked_bindings()
# (testthat >= 3.1.4) and inspect what would have been passed, rather than
# actually running the full pipeline (which needs real BAMs/genome
# resources).

test_that("canope() correctly maps nested-YAML settings to run_canope() arguments (regression)", {
  testthat::skip_if_not_installed("testthat", "3.1.4")
  testthat::skip_if_not_installed("yaml")

  # This mirrors the *documented* nested config style from the README
  # (input: / output: / settings:) -- before this fix, p_value/Tnum/D/min_cor
  # supplied here were silently ignored in favour of hardcoded defaults.
  cfg <- list(
    input = list(bed = "dummy.bed", bamdir = "dummy_bam_dir"),
    output = list(dir = tempdir(), prefix = "TESTRUN"),
    settings = list(
      bsgenome_pkg = "BSgenome.Hsapiens.UCSC.hg19",
      p_value = 1e-12,
      Tnum = 15,
      D = 250000,
      min_cor = 0.87
    )
  )
  cfg_path <- tempfile(fileext = ".yaml")
  yaml::write_yaml(cfg, cfg_path)
  on.exit(unlink(cfg_path))

  captured <- NULL
  testthat::local_mocked_bindings(
    run_canope = function(...) { captured <<- list(...); invisible(NULL) },
    .package = "CANOPE"
  )

  canope(cfg_path)

  expect_equal(captured$p_value, 1e-12)
  expect_equal(captured$Tnum, 15L)
  expect_equal(captured$D, 250000L)
  expect_equal(captured$min_cor, 0.87)
})

test_that("canope() still works with the flat config style", {
  testthat::skip_if_not_installed("testthat", "3.1.4")
  testthat::skip_if_not_installed("yaml")

  cfg <- list(
    bsgenome_pkg = "BSgenome.Hsapiens.UCSC.hg19",
    bed_file = "dummy.bed",
    samples = "dummy_bam_dir",
    p_value = 1e-09,
    Tnum = 8,
    D = 50000,
    min_cor = 0.9
  )
  cfg_path <- tempfile(fileext = ".yaml")
  yaml::write_yaml(cfg, cfg_path)
  on.exit(unlink(cfg_path))

  captured <- NULL
  testthat::local_mocked_bindings(
    run_canope = function(...) { captured <<- list(...); invisible(NULL) },
    .package = "CANOPE"
  )

  canope(cfg_path)

  expect_equal(captured$p_value, 1e-09)
  expect_equal(captured$Tnum, 8L)
  expect_equal(captured$D, 50000L)
  expect_equal(captured$min_cor, 0.9)
})

test_that("canope() requires bed_file, samples, and a GC source", {
  testthat::skip_if_not_installed("yaml")

  cfg_path <- tempfile(fileext = ".yaml")
  yaml::write_yaml(list(bed_file = "x.bed"), cfg_path)   # missing samples & GC source
  on.exit(unlink(cfg_path))

  expect_error(canope(cfg_path), "samples directory is required|fasta_file or bsgenome_pkg")
})

test_that("canope() errors clearly when the config file cannot be found", {
  expect_error(canope(tempfile(fileext = ".yaml")), "Configuration file not found")
})
