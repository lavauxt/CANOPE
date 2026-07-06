# Tests for canoes_legacy_engine.R -- the deliberately-faithful port of the
# original published CANOES HMM core (see file header: "kept faithful on
# purpose"). Not modified in this review pass; these tests only verify its
# documented behaviour and the same forward-backward/constrained-likelihood
# identities used for the modern engine in test-hmm-engine.R.

test_that("legacy_get_distances: single row returns 0, multi-row is start-to-start with no floor", {
  d1 <- legacy_get_distances(data.frame(target = 1L, chromosome = 1, start = 1000))
  expect_equal(d1$distance, 0)

  same_chr <- data.frame(target = 1:3, chromosome = c(1, 1, 1), start = c(1000, 1500, 3000))
  d2 <- legacy_get_distances(same_chr)
  expect_equal(d2$distance, c(0, 500, 1500))

  cross_chr <- data.frame(target = 1:3, chromosome = c(1, 1, 2), start = c(1000, 1500, 100))
  d3 <- legacy_get_distances(cross_chr)
  expect_equal(d3$distance[1:2], c(0, 500))
  expect_equal(d3$distance[3], 100 - 1500 + 1e12 * (2 - 1))
})

test_that("legacy_transition_logprobs rows are proper probability distributions", {
  for (d in c(0, 1000, 1e6)) {
    m <- legacy_transition_logprobs(d, p = 0.01, Tnum = 6, D = 1e5)
    expect_equal(dim(m), c(3L, 3L))
    expect_equal(rowSums(exp(m)), c(1, 1, 1), tolerance = 1e-8)
  }
})

test_that("legacy_phred matches the original (unguarded) formula, including its documented NaN edge case", {
  expect_equal(legacy_phred(0), 0)
  expect_equal(legacy_phred(0.99), 20)
  expect_equal(legacy_phred(1), 99)
  # Deliberately unguarded: prob slightly > 1 makes 1 - prob negative, and
  # log10() of a negative number is NaN in R -- this is documented, expected
  # behaviour for this engine, not a bug.
  expect_true(is.nan(legacy_phred(1 + 1e-6)))
})

make_synthetic_legacy_case <- function() {
  em <- data.frame(
    target     = 1:6,
    delprob    = c(-20, -20, -0.001, -0.001, -20, -20),
    normalprob = c(-0.001, -0.001, -20, -20, -0.001, -0.001),
    dupprob    = c(-25, -25, -25, -25, -25, -25)
  )
  distances <- data.frame(target = 1:6, distance = c(0, 5000, 5000, 5000, 5000, 5000))
  list(em = em, distances = distances, p = 0.01, Tnum = 6, D = 1e5)
}

test_that("legacy forward-backward marginal likelihood is position-invariant", {
  s <- make_synthetic_legacy_case()
  fwd <- legacy_forward_matrix(s$em, s$distances, s$p, s$Tnum, s$D)
  bwd <- legacy_backward_matrix(s$em, s$distances, s$p, s$Tnum, s$D)

  liks <- vapply(1:6, function(i) legacy_likelihood_at(fwd, bwd, i), numeric(1))
  expect_equal(liks, rep(liks[1], 6), tolerance = 1e-6)
})

test_that("legacy_viterbi recovers the encoded DEL run", {
  s <- make_synthetic_legacy_case()
  decoded <- legacy_viterbi(s$em, s$distances, s$p, s$Tnum, s$D)

  expect_equal(decoded$target, 1:6)
  expect_equal(decoded$viterbi_state,
               c(legacy_state_normal, legacy_state_normal, legacy_state_deletion,
                 legacy_state_deletion, legacy_state_normal, legacy_state_normal))
})

test_that("legacy_modified_likelihood can never exceed the unmodified baseline", {
  s <- make_synthetic_legacy_case()
  fwd <- legacy_forward_matrix(s$em, s$distances, s$p, s$Tnum, s$D)
  bwd <- legacy_backward_matrix(s$em, s$distances, s$p, s$Tnum, s$D)

  res <- legacy_modified_likelihood(
    fwd, bwd, s$em, s$distances,
    start_target = 3, end_target = 4,
    disallowed_states = legacy_state_deletion,
    p = s$p, Tnum = s$Tnum, D = s$D
  )
  expect_length(res, 2)
  expect_lte(res[1], res[2] + 1e-8)   # modified (constrained) <= unmodified baseline
})

test_that("legacy_emission_probs has the same column layout as the modern engine's emission_probs (values may legitimately differ)", {
  target_means <- rep(200, 4)
  var_estimate <- rep(250, 4)
  em <- legacy_emission_probs(rep(200, 4), target_means, var_estimate, targets = 1:4)
  expect_equal(colnames(em), c("target", "delprob", "normalprob", "dupprob"))
  expect_true(all(em[, "normalprob"] > em[, "delprob"]))
})
