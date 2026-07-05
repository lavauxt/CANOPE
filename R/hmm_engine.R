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
  # Force numeric mode (YAML often parses 1e-08 as a string)
  p <- as.numeric(p)
  expected_targets <- as.numeric(expected_targets)
  span_bp <- as.numeric(span_bp)
  prior_abnormal <- as.numeric(prior_abnormal)

  stopifnot(p > 0, expected_targets > 0, span_bp > 0, prior_abnormal > 0)
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
build_cnv_hmm <- function(emission_matrix, distances, p, expected_targets, span_bp) {
  params <- hmm_params(p, expected_targets, span_bp)
  dist_vec <- if (is.data.frame(distances)) distances$distance else distances
  structure(
    list(
      emissions = emission_matrix,
      distances = dist_vec,
      params = params,
      n_targets = nrow(emission_matrix)
    ),
    class = "cnv_hmm"
  )
}


#' @export
decode_hmm_states <- function(
    emission_matrix,
    distances,
    params = NULL,
    transition_model = c("distance", "stationary"),
    backend = c("native", "HMM"),
    p = NULL,
    expected_targets = NULL,
    span_bp = NULL
) {
  transition_model <- match.arg(transition_model)
  backend <- match.arg(backend)

  if (is.null(params)) {
    if (is.null(p) || is.null(expected_targets))
      stop("Supply params or p, expected_targets, and span_bp.")
    params <- hmm_params(p, expected_targets, span_bp %||% 100000)
  } else {
    params <- .coerce_hmm_params(params, p, expected_targets, span_bp)
  }

  stationary <- transition_model == "stationary"
  if (stationary && backend == "HMM") {
    return(.decode_stationary_hmm_package(emission_matrix, params))
  }
  hmm_viterbi(emission_matrix, distances, params, stationary = stationary)
}


#' @noRd
.decode_stationary_hmm_package <- function(emission_matrix, params) {
  hmm_obj <- build_hmm_object(emission_matrix, params$p, params$expected_targets)
  decoded <- HMM::viterbi(hmm_obj, as.character(seq_len(nrow(emission_matrix))))
  states <- rep(state_normal, length(decoded))
  states[decoded == "DEL"] <- state_deletion
  states[decoded == "DUP"] <- state_duplication
  data.frame(
    target = emission_matrix[, "target"],
    state = states,
    viterbi_state = states
  )
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


#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


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
add_two_probabilities <- function(x, y) {
  if (is.infinite(x)) return(y)
  if (is.infinite(y)) return(x)
  max(x, y) + log1p(exp(-abs(x - y)))
}


#' @export
build_hmm_object <- function(emission_matrix, p, expected_targets) {
  params <- hmm_params(p, expected_targets, span_bp = Inf)
  symbols <- as.character(seq_len(nrow(emission_matrix)))
  q <- params$q
  pi_mat <- matrix(
    c(1 - q - p, q, p, p, 1 - 2 * p, p, p, q, 1 - q - p),
    nrow = 3, byrow = TRUE
  )
  pi_mat <- pmax(pi_mat, .Machine$double.eps)
  colnames(pi_mat) <- rownames(pi_mat) <- state_names

  log_em <- .hmm_emission_matrix(emission_matrix)
  shifted <- log_em - apply(log_em, 1, max)
  emission_probs_mat <- t(exp(shifted))
  colnames(emission_probs_mat) <- symbols
  rownames(emission_probs_mat) <- state_names

  HMM::initHMM(
    States = state_names,
    Symbols = symbols,
    startProbs = exp(.hmm_initial_logprobs(params)),
    transProbs = pi_mat,
    emissionProbs = emission_probs_mat
  )
}
