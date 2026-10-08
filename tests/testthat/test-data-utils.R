# Tests for the self-contained helpers in data_utils.R.
# (compute_gc_from_fasta / compute_gc_from_bed / get_coverage_from_bams* need
# real genome/BAM resources and are integration-level; not unit tested here.)

test_that("format_chr_label maps every accepted input style to UCSC-style labels", {
  expect_equal(
    format_chr_label(c("1", "chr2", "X", "23", "Y", "24", "M", "MT", "chrM", "chrY")),
    c("chr1", "chr2", "chrX", "chrX", "chrY", "chrY", "chrM", "chrM", "chrM", "chrY")
  )
})

test_that("format_chr_label leaves already-prefixed non-special chromosomes untouched", {
  expect_equal(format_chr_label("chr7"), "chr7")
  expect_equal(format_chr_label("7"), "chr7")
})

test_that("targets_to_rows maps target ids to row positions", {
  target_vector <- c(101L, 102L, 103L, 104L, 105L)
  expect_equal(targets_to_rows(c(102L, 105L), target_vector), c(2L, 5L))
})

test_that("targets_to_rows warns and returns NA for unmatched target ids", {
  target_vector <- c(101L, 102L, 103L)
  expect_warning(rows <- targets_to_rows(c(102L, 999L), target_vector), "could not be matched")
  expect_equal(rows, c(2L, NA_integer_))
})

test_that("clean_name strips path, directory, and (only) the first dot-separated suffix", {
  expect_equal(clean_name("/data/bams/sample1.sorted.bam"), "sample1")
  expect_equal(clean_name("sampleA.bam"), "sampleA")
  expect_equal(clean_name("/x/y/Patient-007.recal.bam"), "Patient-007")
})

test_that("GC values align to target coordinates and report unmatched targets", {
  gc_data <- data.frame(
    chromosome = c("chr2", "chr1"),
    start = c(20L, 10L),
    end = c(30L, 15L),
    GC_CONTENT = c(0.6, 0.4)
  )
  targets <- data.frame(
    chromosome = c("chr1", "chr2", "chr3"),
    start = c(10L, 20L, 5L),
    end = c(15L, 30L, 8L)
  )

  expect_warning(
    gc <- .align_gc_to_targets(gc_data, targets),
    "GC content is missing for 1 target"
  )
  expect_equal(gc, c(0.4, 0.6, NA_real_))
})

test_that("GC values align when rows are reordered but lengths match", {
  gc_data <- data.frame(
    chromosome = c("chr2", "chr1"),
    start = c(20L, 10L),
    end = c(30L, 15L),
    GC_CONTENT = c(0.6, 0.4)
  )
  targets <- gc_data[2:1, c("chromosome", "start", "end")]

  expect_equal(.align_gc_to_targets(gc_data, targets), c(0.4, 0.6))
})

test_that("GC values with mismatched row counts require coordinate columns", {
  targets <- data.frame(chromosome = c("chr1", "chr2"), start = 1:2, end = 3:4)
  expect_error(
    .align_gc_to_targets(data.frame(GC_CONTENT = 0.5), targets),
    "coordinate columns are unavailable for alignment"
  )
})

test_that("compute_gc_from_bed reports missing required BED columns", {
  expect_error(
    compute_gc_from_bed("BSgenome.not.installed", data.frame(other = 1)),
    "missing required column"
  )
})
