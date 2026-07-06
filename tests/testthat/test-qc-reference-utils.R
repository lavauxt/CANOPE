# Tests for qc_reference_utils.R

test_that("infer_reference_samples requires at least 2 samples", {
  m <- matrix(1:10, ncol = 1)
  expect_error(infer_reference_samples(m), "at least 2 samples")
})

test_that("infer_reference_samples returns well-formed, self-consistent output", {
  set.seed(42)
  n_targets <- 40L
  base_signal <- round(200 + 40 * sin(seq_len(n_targets) / 3))

  counts <- sapply(1:5, function(i) pmax(base_signal + rnorm(n_targets, sd = 5), 1))
  colnames(counts) <- paste0("S", 1:5)

  res <- infer_reference_samples(counts, top_n = 3)

  expect_named(res, c("refs", "clean_samples", "anomaly_score", "normalized_counts"))
  expect_true(all(names(res$refs) == res$clean_samples))

  for (s in names(res$refs)) {
    # a sample must never be its own reference
    expect_false(s %in% res$refs[[s]])
    # never more than top_n references
    expect_lte(length(res$refs[[s]]), 3L)
    # every listed reference must be a real sample column
    expect_true(all(res$refs[[s]] %in% colnames(counts)))
  }

  expect_equal(length(res$anomaly_score), ncol(counts))
  expect_equal(dim(res$normalized_counts), dim(counts))
})

test_that("detect_outlier_samples flags a sample with much higher cross-target noise", {
  n_targets <- 20L
  small_wave <- rep(c(0, 1, 0, 2), length.out = n_targets)      # low-noise "normal" samples
  big_wave   <- rep(c(0, 90, 0, 95), length.out = n_targets)    # one wildly noisy sample

  counts <- data.frame(
    S1 = 100 + small_wave,
    S2 = 100 + rev(small_wave),
    S3 = 100 + small_wave + 1,
    S4 = 100 + rev(small_wave) + 1,
    S5 = 100 + big_wave
  )

  res <- detect_outlier_samples(counts, z_threshold = 3)

  expect_equal(nrow(res), 5L)
  expect_true(res$is_outlier[res$sample == "S5"])
  expect_false(any(res$is_outlier[res$sample != "S5"]))
  expect_true(res$robust_z[res$sample == "S5"] > 3)
})

test_that("detect_problematic_exons flags low-mean exons but keeps at least one exon per chromosome", {
  # chr1: three exons, all with mean < min_mean (all "bad") -- the least
  # noisy of the three must be rescued so chr1 isn't wiped out entirely.
  # chr2: four exons, only one with mean < min_mean, and it has three
  # perfectly good neighbours -- that one exon should stay flagged.
  count_matrix <- rbind(
    matrix(c(1, 1, 1, 1), nrow = 1),     # chr1, row 1: mean 1
    matrix(c(2, 2, 2, 2), nrow = 1),     # chr1, row 2: mean 2
    matrix(c(3, 3, 3, 3), nrow = 1),     # chr1, row 3: mean 3
    matrix(c(100, 105, 95, 110), nrow = 1),  # chr2, row 4: good
    matrix(c(100, 105, 95, 110), nrow = 1),  # chr2, row 5: good
    matrix(c(100, 105, 95, 110), nrow = 1),  # chr2, row 6: good
    matrix(c(1, 1, 1, 1), nrow = 1)      # chr2, row 7: mean 1 (bad)
  )
  chromosomes <- c("chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2")

  res <- detect_problematic_exons(count_matrix, chromosomes = chromosomes, min_mean = 20)

  expect_equal(res$problematic, c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))
})

test_that("detect_problematic_exons flags exons with extreme GC content", {
  count_matrix <- rbind(
    c(100, 105, 95, 110),
    c(100, 105, 95, 110),
    c(100, 105, 95, 110)
  )
  gc <- c(0.5, 0.96, 0.4)   # row 2 is out of the default [0.10, 0.90] range

  res <- detect_problematic_exons(count_matrix, gc = gc)
  expect_equal(res$problematic, c(FALSE, TRUE, FALSE))
})

test_that("detect_problematic_exons validates gc length", {
  count_matrix <- rbind(c(1, 2), c(3, 4))
  expect_error(detect_problematic_exons(count_matrix, gc = c(0.5)), "one value per row")
})

test_that("gc_correct_counts is a no-op on small inputs (below the 50-point LOESS minimum)", {
  set.seed(1)
  counts <- matrix(sample(50:150, 20, replace = TRUE), nrow = 10, ncol = 2)
  gc <- runif(10, 0.3, 0.6)

  corrected <- gc_correct_counts(counts, gc)
  expect_equal(dim(corrected), dim(counts))
  expect_equal(corrected, matrix(as.integer(counts), nrow = 10, ncol = 2))
})

test_that("gc_correct_counts validates gc length and never returns negative counts", {
  counts <- matrix(1:10, ncol = 2)
  expect_error(gc_correct_counts(counts, gc = c(0.5)), "one value per row")

  corrected <- gc_correct_counts(counts, gc = runif(5, 0.2, 0.8))
  expect_true(all(corrected >= 0))
})
