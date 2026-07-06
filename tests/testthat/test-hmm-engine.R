# Tests for hmm_engine.R.
#
# hmm_engine.R implements the core HMM math and is intentionally left
# unmodified in this review pass (per earlier direction on this file) -- these
# tests exercise its documented behaviour and well-known HMM identities
# (transition rows are proper probability distributions; forward-backward
# marginal likelihood is position-invariant; constrained likelihood can only
# be <= the unconstrained baseline) rather than hard-coded reference numbers,
# so they remain valid regardless of exactly how the math is implemented
# internally.

test_that("hmm_params validates its inputs and computes q correctly", {
  p <- hmm_params(0.01, expected_targets = 6, span_bp = 1e5)
  expect_equal(p$q, 1 / 6)
  expect_s3_class(p, "cnv_hmm_params")

  expect_error(hmm_params(-0.01, 6, 1e5))
  expect_error(hmm_params(0.01, 0, 1e5))
  expect_error(hmm_params(0.01, 6, -1))
})

test_that("hmm_transition_logprobs always returns valid rows that sum to 1 in probability space", {
  params <- hmm_params(0.01, expected_targets = 6, span_bp = 1e5)
  for (d in c(0, 1, 100, 1e5, 1e7)) {
    m <- hmm_transition_logprobs(d, params)
    expect_equal(dim(m), c(3L, 3L))
    expect_equal(rowSums(exp(m)), c(1, 1, 1), tolerance = 1e-8)
  }
})

test_that("hmm_transition_logprobs is distance-independent when stationary or when span_bp is non-finite", {
  params <- hmm_params(0.01, expected_targets = 6, span_bp = 1e5)
  a <- hmm_transition_logprobs(0,      params, stationary = TRUE)
  b <- hmm_transition_logprobs(999999, params, stationary = TRUE)
  expect_equal(a, b)

  params_inf <- hmm_params(0.01, expected_targets = 6, span_bp = Inf)
  c1 <- hmm_transition_logprobs(0,      params_inf)
  c2 <- hmm_transition_logprobs(999999, params_inf)
  expect_equal(c1, c2)
})

test_that("log_sum_exp matches direct computation and handles -Inf correctly", {
  expect_equal(log_sum_exp(c(0, 0)), log(2))
  expect_equal(log_sum_exp(5), 5)
  expect_equal(log_sum_exp(c(-Inf, -Inf)), -Inf)
  expect_equal(log_sum_exp(c(-Inf, 3)), 3)
})

test_that("add_two_probabilities matches log-space addition of two probabilities", {
  expect_equal(add_two_probabilities(-Inf, -Inf), -Inf)
  expect_equal(add_two_probabilities(log(0.5), log(0.5)), log(1), tolerance = 1e-9)
  expect_equal(add_two_probabilities(log(0.3), -Inf), log(0.3))
})

## ---- forward/backward/viterbi on a small synthetic scenario ----

make_synthetic_hmm_case <- function() {
  # Targets 1,2,5,6 strongly favour NORMAL; targets 3,4 strongly favour DEL;
  # DUP is never favoured anywhere. The gap (~19 nats) is made large enough
  # to dominate any reasonable transition penalty for entering/exiting the
  # abnormal state twice, so the Viterbi path should reliably reflect it.
  em <- data.frame(
    target           = 1:6,
    del_log_prob     = c(-20,  -20, -0.001, -0.001, -20,  -20),
    normal_log_prob  = c(-0.001, -0.001, -20, -20, -0.001, -0.001),
    dup_log_prob     = c(-25, -25, -25, -25, -25, -25)
  )
  distances <- c(0, 5000, 5000, 5000, 5000, 5000)
  params <- hmm_params(0.01, expected_targets = 6, span_bp = 1e5)
  list(em = em, distances = distances, params = params)
}

test_that("forward-backward marginal likelihood is the same at every position (HMM identity)", {
  s <- make_synthetic_hmm_case()
  fwd <- hmm_forward(s$em, s$distances, s$params)
  bwd <- hmm_backward(s$em, s$distances, s$params)

  liks <- vapply(1:6, function(i) hmm_likelihood_at(fwd, bwd, i), numeric(1))
  expect_equal(liks, rep(liks[1], 6), tolerance = 1e-6)
})

test_that("hmm_viterbi recovers the encoded DEL run in the middle of the synthetic scenario", {
  s <- make_synthetic_hmm_case()
  decoded <- hmm_viterbi(s$em, s$distances, s$params)

  expect_equal(decoded$target, 1:6)
  expect_true(all(decoded$state %in% c(state_deletion, state_normal, state_duplication)))
  expect_equal(decoded$state, c(state_normal, state_normal, state_deletion, state_deletion,
                                state_normal, state_normal))
})

test_that("hmm_constrained_likelihood can never exceed the unconstrained baseline", {
  s <- make_synthetic_hmm_case()
  fwd <- hmm_forward(s$em, s$distances, s$params)
  bwd <- hmm_backward(s$em, s$distances, s$params)

  res <- hmm_constrained_likelihood(
    fwd, bwd, s$em, s$distances,
    start_target = 3, end_target = 4,
    forbidden_states = state_deletion,
    params = s$params
  )
  expect_length(res, 2)
  expect_lte(res[1], res[2] + 1e-8)   # constrained <= baseline (+ float tolerance)
})

test_that("decode_hmm_states requires either params or (p, expected_targets)", {
  s <- make_synthetic_hmm_case()
  expect_error(decode_hmm_states(s$em, s$distances), "Supply params")

  out <- decode_hmm_states(s$em, s$distances, p = 0.01, expected_targets = 6, span_bp = 1e5)
  expect_equal(nrow(out), 6L)
})

test_that("emission_probs favours the state whose mean matches the observed count", {
  target_means <- rep(200, 4)
  var_estimate <- rep(250, 4)   # modest overdispersion (dispersion = 50 > 0)
  targets <- 1:4

  em_normal <- emission_probs(rep(200, 4), target_means, var_estimate, targets)
  expect_equal(colnames(em_normal), c("target", "del_log_prob", "normal_log_prob", "dup_log_prob"))
  expect_equal(em_normal[, "target"], targets)
  expect_true(all(em_normal[, "normal_log_prob"] > em_normal[, "del_log_prob"]))
  expect_true(all(em_normal[, "normal_log_prob"] > em_normal[, "dup_log_prob"]))

  em_del <- emission_probs(rep(100, 4), target_means, var_estimate, targets)
  expect_true(all(em_del[, "del_log_prob"] > em_del[, "normal_log_prob"]))
})
