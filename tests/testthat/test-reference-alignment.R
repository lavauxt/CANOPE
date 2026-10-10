test_that("reference counts are aligned to target genomic coordinates", {
  targets <- data.frame(
    chromosome = c("chr1", "chr2"), start = c(10, 20), end = c(15, 25)
  )
  references <- data.frame(
    chromosome = c("chr2", "chr1"), start = c(20, 10), end = c(25, 15),
    REF = c(200, 100)
  )

  aligned <- CANOPE:::.align_reference_counts(targets, references)
  expect_equal(aligned$REF, c(100, 200))
})

test_that("reference alignment rejects ambiguous or mismatched coordinates", {
  targets <- data.frame(
    chromosome = c("chr1", "chr1"), start = c(10, 10), end = c(15, 15)
  )
  references <- data.frame(
    chromosome = c("chr1", "chr1"), start = c(10, 11), end = c(15, 16), REF = 1:2
  )

  expect_error(CANOPE:::.align_reference_counts(targets, references), "Duplicate")
  targets$start <- c(10, 20)
  expect_error(CANOPE:::.align_reference_counts(targets, references), "different genomic coordinates")
})
