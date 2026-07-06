# Tests for canope_confidence.R

make_calls <- function() {
  data.frame(
    CNV      = c("DEL", "DUP", "DEL", "DEL",  "DEL", "DEL"),
    Q_SOME   = c(60,     30,    10,    80,     90,    99),
    NUM_REFS = c(12,     6,     2,     15,     20,    99),
    NUM_TARG = c(5,      2,     1,     3,      4,     99),
    MLCN     = c(1,      3,     1,     1,      3,     0),
    GENE     = c("BRCA1","TP53","TP53","STRC", "EGFR","EGFR"),
    stringsAsFactors = FALSE
  )
}

test_that("score_canope_confidence assigns HIGH/MEDIUM/LOW per the documented rules", {
  calls <- make_calls()
  out <- score_canope_confidence(calls)

  expect_equal(
    out$Confidence,
    c(
      "HIGH",   # row 1: DEL, MLCN=1 consistent, all thresholds cleared, gene not flagged
      "MEDIUM", # row 2: DUP, MLCN=3 consistent, but Q_SOME below the HIGH bar
      "LOW",    # row 3: everything below even the MEDIUM bar
      "LOW",    # row 4: otherwise HIGH-worthy, but STRC is in low_confidence_genes -- always LOW
      "MEDIUM", # row 5: MLCN=3 is NOT consistent with a DEL call -> capped at MEDIUM
      "HIGH"    # row 6: MLCN=0 (homozygous deletion) is always consistent
    )
  )
})

test_that("score_canope_confidence maps MLCN to human-readable CN_label", {
  calls <- data.frame(
    CNV = rep("DEL", 5), Q_SOME = rep(99, 5), NUM_REFS = rep(99, 5),
    NUM_TARG = rep(99, 5), MLCN = c(0, 1, 2, 3, 5), GENE = rep("EGFR", 5),
    stringsAsFactors = FALSE
  )
  out <- score_canope_confidence(calls)
  expect_equal(
    out$CN_label,
    c("Hom. Del (CN=0)", "Hem. Del (CN=1)", "Normal (CN=2)",
      "Dup (CN=3)", "High-level Dup (CN=5)")
  )
})

test_that("score_canope_confidence recognises flagged genes joined by ';' or ', '", {
  calls <- data.frame(
    CNV = c("DEL", "DEL"), Q_SOME = c(99, 99), NUM_REFS = c(99, 99),
    NUM_TARG = c(99, 99), MLCN = c(1, 1),
    GENE = c("EGFR;STRC", "EGFR, STRC"),
    stringsAsFactors = FALSE
  )
  out <- score_canope_confidence(calls)
  expect_equal(out$Confidence, c("LOW", "LOW"))
})

test_that("score_canope_confidence handles missing columns and empty input gracefully", {
  expect_null(score_canope_confidence(NULL))

  empty <- data.frame(CNV = character(0))
  expect_identical(score_canope_confidence(empty), empty)

  incomplete <- data.frame(CNV = "DEL", MLCN = 1)
  expect_warning(out <- score_canope_confidence(incomplete), "missing columns")
  expect_equal(out$Confidence, "LOW")
  expect_true(is.na(out$CN_label))
})
