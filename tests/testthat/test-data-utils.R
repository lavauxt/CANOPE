test_that(".align_gc_to_targets matches by coordinates whatever the row order", {
  targets <- data.frame(chromosome = c("chr1", "chr1", "chr2"), start = c(100, 500, 100),
                        end = c(200, 600, 200), stringsAsFactors = FALSE)
  gc <- data.frame(chromosome = c("chr2", "chr1", "chr1"), start = c(100, 500, 100),
                   end = c(200, 600, 200), GC_CONTENT = c(0.3, 0.2, 0.1))
  expect_equal(.align_gc_to_targets(gc, targets), c(0.1, 0.2, 0.3))
})

test_that(".align_gc_to_targets applies a uniform +1 start shift (0- vs 1-based) but never repairs partial matches", {
  targets <- data.frame(chromosome = "chr1", start = c(100, 500), end = c(200, 600))
  gc <- data.frame(chromosome = "chr1", start = c(101, 501), end = c(200, 600), GC_CONTENT = c(0.4, 0.6))
  expect_equal(suppressMessages(.align_gc_to_targets(gc, targets)), c(0.4, 0.6))

  gc2 <- data.frame(chromosome = "chr1", start = c(100, 501), end = c(200, 600), GC_CONTENT = c(0.4, 0.6))
  expect_warning(res <- .align_gc_to_targets(gc2, targets), "missing for 1")
  expect_equal(res, c(0.4, NA))
})

test_that(".align_gc_to_targets refuses conflicting duplicate coordinates", {
  targets <- data.frame(chromosome = "chr1", start = 100, end = 200)
  gc <- data.frame(chromosome = "chr1", start = c(100, 100), end = c(200, 200), GC_CONTENT = c(0.4, 0.5))
  expect_error(.align_gc_to_targets(gc, targets), "conflicting")
})

test_that("assign_exon_numbers_per_gene keeps row order and leaves gene-less rows unnumbered", {
  bed <- data.frame(chromosome = c("chr1", "chr1", "chr1", "chr1"),
                    start = c(900, 100, 500, 300), end = c(950, 150, 550, 350),
                    GENE = c("A", "A", NA, "A"), stringsAsFactors = FALSE)
  out <- assign_exon_numbers_per_gene(bed)
  expect_equal(out$start, bed$start)
  expect_equal(out$exon_number, c(3L, 1L, NA, 2L))
})

test_that("check_background_calibration only flags a significant excess of outliers", {
  ok  <- check_background_calibration(ratio = rep(1, 40), lo = rep(0.9, 40), hi = rep(1.1, 40), is_affected = rep(FALSE, 40))
  bad <- check_background_calibration(ratio = c(rep(2, 12), rep(1, 28)), lo = rep(0.9, 40), hi = rep(1.1, 40), is_affected = rep(FALSE, 40))
  expect_false(ok$flag)
  expect_true(bad$flag)
  expect_equal(bad$n_outside, 12)
})
