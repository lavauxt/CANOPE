# Tests for call_cnvs.R

test_that("get_distances: single row is 0, same-chromosome is start-to-end floored at 0, cross-chromosome is the sentinel", {
  d1 <- get_distances(data.frame(target = 1L, chromosome = 1L, start = 100L, end = 200L))
  expect_equal(d1$distance, 0)

  same_chr <- data.frame(target = 1:3, chromosome = c(1L, 1L, 1L),
                          start = c(100L, 300L, 250L), end = c(200L, 400L, 500L))
  d2 <- get_distances(same_chr)
  expect_equal(d2$distance[2], 100L)   # 300 - 200
  expect_equal(d2$distance[3], 0L)     # 250 - 400 is negative -> floored at 0 (overlapping targets)

  cross_chr <- data.frame(target = 1:2, chromosome = c(1L, 2L), start = c(100L, 50L), end = c(200L, 100L))
  d3 <- get_distances(cross_chr)
  expect_equal(d3$distance[2], 1e12)
})

## ---- regression tests for the new n >= 2 guards ----

test_that("call_cnvs errors clearly (rather than mis-behaving) with fewer than 2 targets after chromosome filtering", {
  counts <- data.frame(
    target = 1L, chromosome = "chrM", start = 1L, end = 100L, gc = 0.4,
    S1 = 100, S2 = 100, S3 = 100, S4 = 100
  )
  expect_error(
    call_cnvs("S1", counts, p = 0.01, Tnum = 3, D = 1e5),
    "At least 2 targets"
  )
})

test_that("call_cnvs errors clearly with fewer than 2 targets after coverage filtering", {
  counts <- data.frame(
    target = 1:2, chromosome = "1", start = c(1L, 300L), end = c(100L, 400L), gc = c(0.4, 0.4),
    S1 = c(100, 1),      # second target fails the sample_name >= 5 filter
    S2 = c(100, 150), S3 = c(105, 148), S4 = c(95, 152)
  )
  expect_error(
    call_cnvs("S1", counts, p = 0.01, Tnum = 3, D = 1e5),
    "At least 2 valid targets"
  )
})

## ---- full synthetic end-to-end integration test ----

test_that("call_cnvs runs end-to-end on synthetic data and detects a large, obvious deletion", {
  n_targets <- 30L
  base <- round(200 + 30 * sin(seq_len(n_targets) / 2.5))
  del_region <- 10:19   # 10 contiguous targets

  counts <- data.frame(
    target     = seq_len(n_targets),
    chromosome = rep("1", n_targets),
    start      = seq(1L, by = 200L, length.out = n_targets),
    end        = seq(150L, by = 200L, length.out = n_targets),
    gc         = rep(0.45, n_targets)
  )

  # Test sample: the shared baseline pattern, halved (heterozygous deletion)
  # over a contiguous 10-target run.
  counts$S1 <- base
  counts$S1[del_region] <- round(base[del_region] * 0.5)

  # Seven reference samples: the same underlying shape, each scaled by a
  # different (but purely multiplicative, rank-preserving) factor -- highly
  # correlated with each other and with S1 everywhere except the deletion.
  ref_scale <- c(S2 = 1.00, S3 = 1.05, S4 = 0.95, S5 = 1.02, S6 = 0.98, S7 = 1.03, S8 = 0.97)
  for (nm in names(ref_scale)) counts[[nm]] <- round(base * ref_scale[[nm]])

  result <- call_cnvs("S1", counts, p = 0.01, Tnum = 3, D = 1e5, numrefs = 7)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("SAMPLE", "CNV", "INTERVAL", "KB", "CHR", "MID_BP", "TARGETS",
                    "NUM_TARG", "GENE", "MLCN", "Q_SOME", "NUM_REFS", "REF_SAMPLES") %in%
                    colnames(result)))

  # A ~50% depth reduction over 10 contiguous targets, against a reference
  # panel that is otherwise a near-perfect (purely rescaled) match, is about
  # as unambiguous a deletion signal as synthetic data can produce.
  expect_gt(nrow(result), 0)
  expect_true(all(result$SAMPLE == "S1"))
  expect_true(any(result$CNV == "DEL"))
})
