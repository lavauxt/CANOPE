#' @export
num_states <- 3L
#' @export
num_abnormal_states <- 2L
#' @export
state_deletion <- 1L
#' @export
state_normal <- 2L
#' @export
state_duplication <- 3L
#' @export
prior_cnv <- 0.0075

state_names <- c("DEL", "NORMAL", "DUP")


#' @export
hmm_params <- function(
    p,
    expected_targets,
    span_bp = 100000,
    prior_abnormal = prior_cnv
) {
  p <- as.numeric(p)
  expected_targets <- as.numeric(expected_targets)
  span_bp <- as.numeric(span_bp)
  prior_abnormal <- as.numeric(prior_abnormal)

  if (length(p) != 1L || !is.finite(p) || p <= 0 || p >= 0.5)
    stop("p must be a single finite value strictly between 0 and 0.5.", call. = FALSE)
  if (length(expected_targets) != 1L || !is.finite(expected_targets) || expected_targets < 1)
    stop("expected_targets must be a single finite value of at least 1.", call. = FALSE)
  if (length(span_bp) != 1L || is.na(span_bp) || span_bp <= 0)
    stop("span_bp must be a single positive value (Inf is allowed).", call. = FALSE)
  if (length(prior_abnormal) != 1L || !is.finite(prior_abnormal) ||
      prior_abnormal <= 0 || prior_abnormal >= 1)
    stop("prior_abnormal must be a single finite value strictly between 0 and 1.", call. = FALSE)
  if (p + 1 / expected_targets > 1)
    stop("p + 1 / expected_targets must not exceed 1.", call. = FALSE)
  structure(
    list(
      p = p,
      expected_targets = expected_targets,
      span_bp = span_bp,
      prior_abnormal = prior_abnormal,
      q = 1 / expected_targets
    ),
    class = "cnv_hmm_params"
  )
}


#' @noRd
.hmm_initial_logprobs <- function(params) {
  abn <- params$prior_abnormal / num_abnormal_states
  log(c(abn, 1 - params$prior_abnormal, abn))
}


#' @export
hmm_transition_logprobs <- function(distance, params, stationary = FALSE) {
  p <- params$p
  q <- params$q
  f <- if (stationary || !is.finite(params$span_bp) || params$span_bp <= 0) {
    1
  } else {
    exp(-distance / params$span_bp)
  }

  p_stay_abn <- f * (1 - q) + (1 - f) * p
  p_abn_norm <- f * q + (1 - f) * (1 - 2 * p)
  p_switch_abn <- (1 - f) * p
  p_norm_norm <- 1 - 2 * p
  p_norm_abn <- p

  probs <- c(
    p_stay_abn, p_abn_norm, p_switch_abn,
    p_norm_abn, p_norm_norm, p_norm_abn,
    p_switch_abn, p_abn_norm, p_stay_abn
  )
  probs <- pmax(probs, .Machine$double.eps)
  log(matrix(probs, num_states, num_states, byrow = TRUE))
}


#' @noRd
.hmm_emission_matrix <- function(emission_matrix) {
  as.matrix(emission_matrix[, c("del_log_prob", "normal_log_prob", "dup_log_prob"), drop = FALSE])
}


#' @noRd
.hmm_step_transition <- function(distances, step, params, stationary) {
  d <- if (is.data.frame(distances)) distances$distance[step] else distances[step]
  hmm_transition_logprobs(d, params, stationary = stationary)
}


#' @export
hmm_forward <- function(emission_matrix, distances, params, stationary = FALSE) {
  params <- .coerce_hmm_params(params)
  em <- .hmm_emission_matrix(emission_matrix)
  n <- nrow(em)
  fwd <- matrix(NA_real_, nrow = n, ncol = num_states)
  fwd[1, ] <- .hmm_initial_logprobs(params) + em[1, ]

  for (i in 2:n) {
    trans <- .hmm_step_transition(distances, i, params, stationary)
    temp <- fwd[i - 1, ] + trans
    fwd[i, ] <- apply(temp, 2L, log_sum_exp) + em[i, ]
  }
  fwd
}


#' @export
hmm_backward <- function(emission_matrix, distances, params, stationary = FALSE) {
  params <- .coerce_hmm_params(params)
  em <- .hmm_emission_matrix(emission_matrix)
  n <- nrow(em)
  bwd <- matrix(NA_real_, nrow = n, ncol = num_states)
  bwd[n, ] <- 0

  for (i in (n - 1):1) {
    trans <- .hmm_step_transition(distances, i + 1L, params, stationary)
    temp <- trans +
      matrix(bwd[i + 1, ], num_states, num_states, byrow = TRUE) +
      matrix(em[i + 1, ], num_states, num_states, byrow = TRUE)
    bwd[i, ] <- apply(temp, 1L, log_sum_exp)
  }
  bwd
}


#' @export
hmm_likelihood_at <- function(forward, backward, position) {
  log_sum_exp(forward[position, ] + backward[position, ])
}


#' @export
hmm_viterbi <- function(emission_matrix, distances, params, stationary = FALSE) {
  params <- .coerce_hmm_params(params)
  targets <- emission_matrix[, "target"]
  em <- .hmm_emission_matrix(emission_matrix)
  n <- nrow(em)

  score <- matrix(NA_real_, nrow = n, ncol = num_states)
  ptr <- matrix(NA_integer_, nrow = n, ncol = num_states)
  score[1, ] <- .hmm_initial_logprobs(params) + em[1, ]

  for (i in 2:n) {
    trans <- .hmm_step_transition(distances, i, params, stationary)
    temp <- score[i - 1, ] + trans
    score[i, ] <- vapply(seq_len(num_states), function(j) {
      col <- temp[, j]
      if (!any(is.finite(col))) -Inf else max(col)
    }, numeric(1)) + em[i, ]
    ptr[i, ] <- vapply(seq_len(num_states), function(j) {
      col <- temp[, j]
      if (!any(is.finite(col))) 1L else which.max(col)
    }, integer(1))
  }

  states <- integer(n)
  states[n] <- which.max(score[n, ])
  for (i in (n - 1):1) states[i] <- ptr[i + 1, states[i + 1]]

  data.frame(target = targets, state = states, viterbi_state = states)
}


#' @export
hmm_constrained_likelihood <- function(
    forward,
    backward,
    emission_matrix,
    distances,
    start_target,
    end_target,
    forbidden_states,
    params,
    stationary = FALSE
) {
  params <- .coerce_hmm_params(params)
  targets <- emission_matrix[, "target"]
  em <- .hmm_emission_matrix(emission_matrix)
  n <- nrow(em)

  left <- min(which(targets >= start_target))
  right <- max(which(targets <= end_target))
  stopifnot(right >= left)

  anchor <- min(right + 1L, n)
  baseline <- hmm_likelihood_at(forward, backward, anchor)
  fwd <- forward
  mod_em <- em
  mod_em[left:right, forbidden_states] <- -Inf

  if (left == 1L) {
    fwd[1L, ] <- .hmm_initial_logprobs(params) + mod_em[1L, ]
    left <- 2L
  }

  for (i in seq(left, anchor)) {
    trans <- .hmm_step_transition(distances, i, params, stationary)
    temp <- fwd[i - 1L, ] + trans
    sums <- apply(temp, 2L, log_sum_exp)
    use_em <- if (i <= right) mod_em[i, ] else em[i, ]
    fwd[i, ] <- sums + use_em
  }

  c(hmm_likelihood_at(fwd, backward, anchor), baseline)
}


#' @export
decode_hmm_states <- function(
    emission_matrix,
    distances,
    params = NULL,
    transition_model = c("distance", "stationary"),
    p = NULL,
    expected_targets = NULL,
    span_bp = NULL
) {
  transition_model <- match.arg(transition_model)

  if (is.null(params)) {
    if (is.null(p) || is.null(expected_targets))
      stop("Supply params or p, expected_targets, and span_bp.")
    params <- hmm_params(p, expected_targets, span_bp %||% 100000)
  } else {
    params <- .coerce_hmm_params(params, p, expected_targets, span_bp)
  }

  stationary <- transition_model == "stationary"
  hmm_viterbi(emission_matrix, distances, params, stationary = stationary)
}


#' @noRd
.coerce_hmm_params <- function(params, p = NULL, expected_targets = NULL, span_bp = NULL) {
  if (inherits(params, "cnv_hmm_params")) return(params)
  if (is.list(params) && !is.null(params$p)) {
    return(hmm_params(
      params$p,
      params$expected_targets %||% params$Tnum,
      params$span_bp %||% params$D %||% span_bp %||% 100000,
      params$prior_abnormal %||% prior_cnv
    ))
  }
  if (!is.null(p) && !is.null(expected_targets)) {
    return(hmm_params(p, expected_targets, span_bp %||% 100000))
  }
  if (is.numeric(params) && length(params) == 3L) {
    return(hmm_params(params[1], params[2], params[3]))
  }
  stop("Invalid HMM parameters.")
}


#' @export
emission_probs <- function(test_counts, target_means, var_estimate, targets) {
  n <- length(test_counts)
  state_means <- t(apply(data.frame(x = target_means), 1, function(x) c(x * 0.5, x, x * 1.5)))
  dispersion <- var_estimate - target_means
  size <- ifelse(dispersion > 0, target_means^2 / dispersion, Inf)

  em <- matrix(NA_real_, nrow = n, ncol = 4L)
  colnames(em) <- c("target", "del_log_prob", "normal_log_prob", "dup_log_prob")
  em[, "del_log_prob"] <- stats::dnbinom(test_counts, mu = state_means[, 1], size = size / 2, log = TRUE)
  em[, "normal_log_prob"] <- stats::dnbinom(test_counts, mu = state_means[, 2], size = size, log = TRUE)
  em[, "dup_log_prob"] <- stats::dnbinom(test_counts, mu = state_means[, 3], size = size * 3 / 2, log = TRUE)
  em[, "target"] <- targets

  bad <- which(apply(em[, 2:4, drop = FALSE], 1L, function(x) all(is.infinite(x))))
  if (length(bad) > 0L) {
    eps <- log(1e-12)
    for (i in bad) {
      if (test_counts[i] >= state_means[i, 3L]) {
        em[i, 2:4] <- c(eps, eps, 0)
      } else if (test_counts[i] <= state_means[i, 1L]) {
        em[i, 2:4] <- c(0, eps, eps)
      } else {
        em[i, 2:4] <- c(eps, 0, eps)
      }
    }
  }
  em
}


#' @export
log_sum_exp <- function(x) {
  max_x <- max(x)
  if (is.infinite(max_x)) return(max_x)
  max_x + log(sum(exp(x - max_x)))
}


#' @export
legacy_num_states <- 3L
#' @export
legacy_num_abnormal_states <- 2L
#' @export
legacy_state_deletion <- 1L
#' @export
legacy_state_normal <- 2L
#' @export
legacy_state_duplication <- 3L
#' @export
legacy_initial_abnormal_prior <- 0.0075


#' @noRd
.legacy_initial_logprobs <- function() {
  log(c(
    legacy_initial_abnormal_prior / legacy_num_abnormal_states,
    1 - legacy_initial_abnormal_prior,
    legacy_initial_abnormal_prior / legacy_num_abnormal_states
  ))
}


#' Inter-Target Distances, Original CANOES Formula
#'
#' Port of \code{GetDistances()}. Unlike \code{\link{get_distances}()} (the
#' modern engine's version, start-to-end and floored at 0), this measures
#' \strong{start-to-start} and applies no floor — exactly the original.
#'
#' @param counts Data frame with \code{target}, \code{chromosome}, \code{start}
#'   columns, sorted by chromosome then start.
#' @return Data frame: \code{target}, \code{distance}.
#' @export
legacy_get_distances <- function(counts) {
  chromosome <- counts[["chromosome"]]
  startbase  <- counts[["start"]]
  n <- length(startbase)
  if (n <= 1L) return(data.frame(target = counts[["target"]], distance = 0))
  distances <- c(
    0,
    startbase[2:n] - startbase[1:(n - 1L)] +
      1000000000000 * (chromosome[2:n] - chromosome[1:(n - 1L)])
  )
  data.frame(target = counts[["target"]], distance = distances)
}


#' Transition Log-Probabilities, Original CANOES Formula
#'
#' Port of \code{GetTransitionMatrix()}. Always distance-decayed (\code{f =
#' exp(-distance/D)}) — the original never had a "stationary" mode.
#'
#' @param distance Numeric. Distance to the previous target.
#' @param p     Numeric. Per-target transition rate.
#' @param Tnum  Numeric. Expected number of targets per CNV.
#' @param D     Numeric. Expected genomic span of a CNV (bp).
#' @return 3x3 log-probability transition matrix (rows = from-state, cols = to-state).
#' @export
legacy_transition_logprobs <- function(distance, p, Tnum, D) {
  q <- 1 / Tnum
  f <- exp(-distance / D)

  prob_abn_abn      <- f * (1 - q) + (1 - f) * p
  prob_abn_norm      <- f * q + (1 - f) * (1 - 2 * p)
  prob_abn_diff_abn <- (1 - f) * p
  prob_norm_norm     <- 1 - 2 * p
  prob_norm_abn      <- p

  probs <- c(
    prob_abn_abn, prob_abn_norm, prob_abn_diff_abn,
    prob_norm_abn, prob_norm_norm, prob_norm_abn,
    prob_abn_diff_abn, prob_abn_norm, prob_abn_abn
  )
  log(matrix(probs, legacy_num_states, legacy_num_states, byrow = TRUE))
}


#' Emission Log-Probabilities, Original CANOES Formula
#'
#' Port of \code{EmissionProbs()}. Same negative-binomial model as
#' \code{\link{emission_probs}()} in the modern engine (state means at
#' 0.5x/1x/1.5x of the expected count, size halved/doubled accordingly); the
#' only behavioural difference is the all-impossible-states edge case, which
#' the original fills with hard \code{-Inf}/\code{-0.01} sentinels rather than
#' the softer finite values the modern engine uses.
#'
#' @param test_counts Integer vector of observed read counts for the test sample (one per target).
#' @param target_means Numeric vector of expected read counts (mean) for each target (from reference panel).
#' @param var_estimate Numeric vector of variance estimates for each target.
#' @param targets Integer vector of target identifiers.
#' @return Matrix with columns \code{target, delprob, normalprob, dupprob} (log-scale).
#' @export
legacy_emission_probs <- function(test_counts, target_means, var_estimate, targets) {
  n <- length(test_counts)
  state_means <- t(apply(data.frame(x = target_means), 1, function(x) c(x * 0.5, x, x * 1.5)))
  size <- target_means^2 / (var_estimate - target_means)

  em <- matrix(NA_real_, n, 4L)
  colnames(em) <- c("target", "delprob", "normalprob", "dupprob")

  size_del <- size / 2
  size_dup <- size * 3 / 2
  em[, "delprob"]    <- stats::dnbinom(test_counts, mu = state_means[, 1], size = size_del, log = TRUE)
  em[, "normalprob"] <- stats::dnbinom(test_counts, mu = state_means[, 2], size = size,     log = TRUE)
  em[, "dupprob"]    <- stats::dnbinom(test_counts, mu = state_means[, 3], size = size_dup, log = TRUE)
  em[, "target"] <- targets

  row_all_inf <- which(apply(em, 1, function(x) all(is.infinite(x))))
  if (length(row_all_inf) > 0) {
    for (i in row_all_inf) {
      if (test_counts[i] >= state_means[i, 3]) {
        em[i, 2:4] <- c(-Inf, -Inf, -0.01)
      } else if (test_counts[i] <= state_means[i, 1]) {
        em[i, 2:4] <- c(-0.01, -Inf, -Inf)
      } else {
        em[i, 2:4] <- c(-Inf, -0.01, -Inf)
      }
    }
  }
  em
}


#' @noRd
legacy_add_two_probabilities <- function(x, y) {
  if (is.infinite(x)) return(y)
  if (is.infinite(y)) return(x)
  max(x, y) + log1p(exp(-abs(x - y)))
}


#' @noRd
legacy_sum_probabilities <- function(x) {
  s <- x[1]
  for (i in 2:length(x)) s <- legacy_add_two_probabilities(s, x[i])
  s
}


#' Viterbi Decoding, Original CANOES Formula
#'
#' Port of \code{Viterbi()}.
#'
#' @param emission_probs_matrix Output of \code{\link{legacy_emission_probs}}.
#' @param distances Output of \code{\link{legacy_get_distances}}.
#' @inheritParams legacy_transition_logprobs
#' @return Data frame: \code{target, viterbi_state} (1 = DEL, 2 = NORMAL, 3 = DUP).
#' @export
legacy_viterbi <- function(emission_probs_matrix, distances, p, Tnum, D) {
  targets <- emission_probs_matrix[, 1]
  em <- as.matrix(emission_probs_matrix[, 2:4, drop = FALSE])
  n <- nrow(em)

  score <- matrix(NA_real_, n, legacy_num_states)
  ptr   <- matrix(NA_integer_, n, legacy_num_states)
  score[1, ] <- .legacy_initial_logprobs() + em[1, ]

  for (i in 2:n) {
    trans <- score[i - 1, ] + legacy_transition_logprobs(distances$distance[i], p, Tnum, D)
    score[i, ] <- apply(trans, 2, max) + em[i, ]
    ptr[i, ]   <- apply(trans, 2, which.max)
  }

  states <- integer(n)
  states[n] <- which.max(score[n, ])
  for (i in (n - 1):1) states[i] <- ptr[i + 1, states[i + 1]]

  data.frame(target = targets, viterbi_state = states)
}


#' Forward Matrix, Original CANOES Formula
#' @inheritParams legacy_viterbi
#' @export
legacy_forward_matrix <- function(emission_probs_matrix, distances, p, Tnum, D) {
  em <- as.matrix(emission_probs_matrix[, 2:4, drop = FALSE])
  n <- nrow(em)
  fwd <- matrix(NA_real_, n, legacy_num_states)
  fwd[1, ] <- .legacy_initial_logprobs() + em[1, ]

  for (i in 2:n) {
    temp <- fwd[i - 1, ] + legacy_transition_logprobs(distances$distance[i], p, Tnum, D)
    fwd[i, ] <- apply(temp, 2, legacy_sum_probabilities) + em[i, ]
  }
  fwd
}


#' Backward Matrix, Original CANOES Formula
#' @inheritParams legacy_viterbi
#' @export
legacy_backward_matrix <- function(emission_probs_matrix, distances, p, Tnum, D) {
  em <- as.matrix(emission_probs_matrix[, 2:4, drop = FALSE])
  n <- nrow(em)
  bwd <- matrix(NA_real_, n, legacy_num_states)
  bwd[n, ] <- rep(0, legacy_num_states)

  for (i in (n - 1):1) {
    temp <- legacy_transition_logprobs(distances$distance[i + 1], p, Tnum, D) +
      matrix(bwd[i + 1, ], 3, 3, byrow = TRUE) +
      matrix(em[i + 1, ], 3, 3, byrow = TRUE)
    bwd[i, ] <- apply(temp, 1, legacy_sum_probabilities)
  }
  bwd
}


#' @noRd
legacy_likelihood_at <- function(forward_matrix, backward_matrix, x) {
  legacy_sum_probabilities(forward_matrix[x, ] + backward_matrix[x, ])
}


#' Constrained Likelihood, Original CANOES Formula
#'
#' Port of \code{GetModifiedLikelihood()}: the likelihood of the data when
#' certain states are disallowed between \code{start_target} and
#' \code{end_target}, used by \code{\link{genotype_cnvs}} for Q_SOME/NQ
#' scoring under \code{engine = "legacy_canoes"}.
#'
#' @inheritParams legacy_viterbi
#' @param forward_matrix,backward_matrix Output of \code{\link{legacy_forward_matrix}} /
#'   \code{\link{legacy_backward_matrix}}.
#' @param start_target,end_target Integer target IDs bounding the region.
#' @param disallowed_states Integer vector of states to forbid in that region.
#' @return \code{c(modified_likelihood, unmodified_likelihood)}.
#' @export
legacy_modified_likelihood <- function(
    forward_matrix, backward_matrix, emission_probs_matrix, distances,
    start_target, end_target, disallowed_states, p, Tnum, D
) {
  targets <- emission_probs_matrix[, 1]
  em <- as.matrix(emission_probs_matrix[, 2:4, drop = FALSE])

  left  <- min(which(targets >= start_target))
  right <- max(which(targets <= end_target))
  n <- nrow(em)
  unmodified_likelihood <- legacy_likelihood_at(forward_matrix, backward_matrix, min(right + 1, n))
  stopifnot(right >= left)

  mod_em <- em
  mod_em[left:right, disallowed_states] <- -Inf

  if (left == 1L) {
    fwd <- forward_matrix
    fwd[1, ] <- .legacy_initial_logprobs() + mod_em[1, ]
    left <- left + 1L
  } else {
    fwd <- forward_matrix
  }

  for (i in seq(left, min(right + 1, n))) {
    temp <- fwd[i - 1, ] + legacy_transition_logprobs(distances$distance[i], p, Tnum, D)
    sums <- apply(temp, 2, legacy_sum_probabilities)
    if (!i == (right + 1)) {
      fwd[i, ] <- sums + mod_em[i, ]
    } else {
      fwd[i, ] <- sums + em[i, ]
    }
  }

  modified_likelihood <- legacy_likelihood_at(fwd, backward_matrix, min(right + 1, n))
  c(modified_likelihood, unmodified_likelihood)
}


#' Phred Score, Original CANOES Formula (no numerical floor)
#'
#' Port of the inline \code{Phred()} closure in \code{GenotypeCNVs()}.
#' \strong{Deliberately} has no protection against \code{1 - prob <= 0}
#' (which the modern engine's call_cnvs.R adds via
#' \code{max(1 - prob, 1e-10)}) — this can return \code{NaN} for \code{prob}
#' slightly outside \code{[0, 1]} due to floating-point error, exactly as the
#' original could. Kept faithful on purpose for this engine.
#'
#' @param prob Numeric probability.
#' @return Integer Phred score, or \code{NA}/\code{NaN} if \code{prob} pushes
#'   \code{1 - prob} to zero or negative.
#' @export
legacy_phred <- function(prob) {
  suppressWarnings(round(min(99, -10 * log10(1 - prob))))
}
