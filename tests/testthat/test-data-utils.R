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
