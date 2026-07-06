# Tests for canope_utils.R

test_that("%||% falls back only on NULL, not on other falsy values", {
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% 2, 2)
  expect_equal(FALSE %||% 2, FALSE)   # FALSE is not NULL -- must not fall through
  expect_equal(0 %||% 2, 0)
  expect_equal(NA %||% 2, NA)
})

test_that("stop_if_missing / stop_if_not_file validate as documented", {
  expect_error(stop_if_missing(NULL, "boom"), "boom")
  expect_error(stop_if_missing(character(0), "boom"), "boom")
  expect_silent(stop_if_missing("ok", "boom"))

  expect_error(stop_if_not_file(NULL, "no file"), "no file")
  expect_error(stop_if_not_file(tempfile(), "no file"), "no file")
  tf <- tempfile()
  file.create(tf)
  on.exit(unlink(tf))
  expect_silent(stop_if_not_file(tf, "no file"))
})

test_that("sanitize_filename replaces unsafe characters only", {
  expect_equal(sanitize_filename("sample 1/2:test"), "sample_1_2_test")
  expect_equal(sanitize_filename("Sample-1.2_ok"), "Sample-1.2_ok")
})

test_that("canope_filter_chromosomes handles chr-prefix mismatches both ways", {
  df <- data.frame(
    CHR = c("chr1", "chr2", "chrX", "3"),
    value = 1:4,
    stringsAsFactors = FALSE
  )

  # include a bare chromosome name, expect it to match the chr-prefixed row
  out <- canope_filter_chromosomes(df, include = "1")
  expect_equal(out$value, 1)

  # exclude using a chr-prefixed name against a bare row
  out2 <- canope_filter_chromosomes(df, exclude = "chr3")
  expect_false(3 %in% out2$value)
  expect_setequal(out2$value, c(1, 2, 4))

  # empty/NULL input passes through unchanged
  expect_null(canope_filter_chromosomes(NULL))
  empty_df <- df[0, ]
  expect_identical(canope_filter_chromosomes(empty_df), empty_df)
})

test_that("normalize_chromosome_vec adopts the reference's chr-prefix style", {
  expect_equal(
    normalize_chromosome_vec(c("1", "chr2", "X"), ref_chromosomes = c("chr1", "chr2")),
    c("chr1", "chr2", "chrX")
  )
  expect_equal(
    normalize_chromosome_vec(c("chr1", "2", "chrX"), ref_chromosomes = c("1", "2")),
    c("1", "2", "X")
  )
})

test_that("parse_canope_interval splits chrom:start-end correctly", {
  p <- parse_canope_interval("chr1:100-200")
  expect_equal(p$chrom, "chr1")
  expect_equal(p$start, 100L)
  expect_equal(p$end, 200L)

  p2 <- parse_canope_interval("22:25713988-25756059")
  expect_equal(p2$chrom, "22")
  expect_equal(p2$start, 25713988L)
  expect_equal(p2$end, 25756059L)
})

test_that("parse_canope_targets splits '..'-joined target ids", {
  expect_equal(parse_canope_targets("42..47"), c(42L, 47L))
  expect_equal(parse_canope_targets("5..5"), c(5L, 5L))
})

## ---- assign_exon_numbers_per_gene(): the core bug-fix regression tests ----

test_that("assign_exon_numbers_per_gene numbers exons in genomic order, independent of input row order", {
  bed_df <- data.frame(
    chromosome = c("chr1", "chr1", "chr1"),
    start      = c(500L, 100L, 300L),
    end        = c(600L, 200L, 400L),
    GENE       = c("GENE3", "GENE3", "GENE3"),
    stringsAsFactors = FALSE
  )

  out <- assign_exon_numbers_per_gene(bed_df)

  expect_equal(nrow(out), 3L)
  # Row order of the *output* must match the *input* order exactly.
  expect_equal(out$start, bed_df$start)
  # Genomic order is 100-200 (exon 1), 300-400 (exon 2), 500-600 (exon 3);
  # the input rows are given in the order 500-600, 100-200, 300-400.
  expect_equal(out$exon_number, c(3L, 1L, 2L))
})

test_that("assign_exon_numbers_per_gene preserves row count and shares exon_number across duplicate rows (regression)", {
  # This is the exact bug class fixed in this pass: previously, duplicate
  # (chrom, start, end, gene) rows were dropped from the *output*, breaking
  # every caller that binds exon_number back onto a same-sized table by
  # position (canope_qc_metrics.R, generate_plots.R, CANOPE_report.Rmd).
  bed_df <- data.frame(
    chromosome = c("chr1", "chr1", "chr1", "chr2"),
    start      = c(300L, 100L, 100L, 50L),
    end        = c(400L, 200L, 200L, 60L),
    GENE       = c("GENE1", "GENE1", "GENE1", "GENE2"),
    stringsAsFactors = FALSE
  )
  # row 3 is an exact duplicate of row 2 (same chrom/start/end/gene)

  expect_warning(
    out <- assign_exon_numbers_per_gene(bed_df),
    "duplicate"
  )

  # Row count and row order must be preserved exactly.
  expect_equal(nrow(out), nrow(bed_df))
  expect_equal(out$start, bed_df$start)
  expect_equal(out$GENE, bed_df$GENE)

  # GENE1: genomic order is 100-200 (exon 1), 300-400 (exon 2).
  # Row 1 (300-400) -> exon 2; rows 2 and 3 (both 100-200, duplicates of
  # each other) -> exon 1, sharing the same number instead of being
  # miscounted as two separate exons.
  expect_equal(out$exon_number, c(2L, 1L, 1L, 1L))
})

test_that("assign_exon_numbers_per_gene works with alternate column name spellings", {
  bed_df <- data.frame(
    Chr   = c("chr1", "chr1"),
    Start = c(200L, 100L),
    End   = c(300L, 200L),
    Gene  = c("G", "G"),
    stringsAsFactors = FALSE
  )
  out <- assign_exon_numbers_per_gene(bed_df)
  expect_equal(out$exon_number, c(2L, 1L))
})

test_that("compute_exon_index reuses an existing exon_number column when present", {
  bed_df <- data.frame(
    chromosome = "chr1", start = 1:3, end = 2:4, GENE = "G",
    exon_number = c(9L, 8L, 7L)
  )
  expect_equal(compute_exon_index(bed_df), c(9L, 8L, 7L))
})

test_that("compute_exon_index computes exon numbers when the column is absent", {
  bed_df <- data.frame(
    chromosome = c("chr1", "chr1"), start = c(200L, 100L), end = c(300L, 200L),
    GENE = c("G", "G"), stringsAsFactors = FALSE
  )
  expect_equal(compute_exon_index(bed_df), c(2L, 1L))
})

## ---- check_background_calibration() ----

test_that("check_background_calibration does not flag a well-calibrated window", {
  ratio       <- rep(0, 7)
  lo          <- rep(-1, 7)
  hi          <- rep(1, 7)
  is_affected <- c(rep(FALSE, 5), TRUE, TRUE)

  res <- check_background_calibration(ratio, lo, hi, is_affected)
  expect_equal(res$n_background, 5L)
  expect_equal(res$n_outside, 0L)
  expect_false(res$flag)
})

test_that("check_background_calibration flags a window where background exons are mostly outside the interval", {
  ratio       <- c(rep(2, 5), 3, 3)     # background ratio (2) far outside [-1, 1]
  lo          <- rep(-1, 7)
  hi          <- rep(1, 7)
  is_affected <- c(rep(FALSE, 5), TRUE, TRUE)

  res <- check_background_calibration(ratio, lo, hi, is_affected)
  expect_equal(res$n_background, 5L)
  expect_equal(res$n_outside, 5L)
  expect_true(res$flag)
})

test_that("check_background_calibration respects min_n and never flags below it", {
  ratio       <- c(2, 2, 2, 3)
  lo          <- rep(-1, 4)
  hi          <- rep(1, 4)
  is_affected <- c(FALSE, FALSE, FALSE, TRUE)   # only 3 background exons

  res <- check_background_calibration(ratio, lo, hi, is_affected, min_n = 5)
  expect_equal(res$n_background, 3L)
  expect_false(res$flag)   # below min_n, regardless of how extreme
})

test_that("check_background_calibration handles zero background exons gracefully", {
  res <- check_background_calibration(1, 0, 2, TRUE)
  expect_equal(res$n_background, 0L)
  expect_true(is.na(res$pct_outside))
  expect_false(res$flag)
})
