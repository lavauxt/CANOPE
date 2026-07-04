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
#' @inheritParams emission_probs
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
  round(min(99, -10 * log10(1 - prob)))
}
