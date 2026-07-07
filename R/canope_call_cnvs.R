#' Compute Inter-Target Distances (for HMM transition modelling)
#'
#' @param counts Data frame with \code{target}, \code{chromosome}, \code{start},
#'   \code{end} columns, already sorted by chromosome then start.
#' @return Data frame: \code{target}, \code{distance} (bp to the previous target,
#'   or a large sentinel value across a chromosome break).
#' @export
get_distances <- function(counts) {
  n <- nrow(counts)
  dst <- numeric(n)
  chr_break <- 1e12
  if (n > 1) {
    for (i in 2:n) {
      if (counts$chromosome[i] == counts$chromosome[i - 1]) {
        dst[i] <- max(counts$start[i] - counts$end[i - 1], 0)
      } else {
        dst[i] <- chr_break
      }
    }
  }
  data.frame(target = counts$target, distance = dst)
}


#' Call CNVs for a Single Sample (CANOPE HMM Engine Wrapper)
#'
#' Builds a weighted reference panel, fits per-target emission probabilities,
#' decodes the most likely CNV state path, and genotypes each resulting call.
#' \strong{The underlying HMM math in \code{hmm_engine.R}
#' (\code{hmm_params}, \code{emission_probs}, \code{decode_hmm_states},
#' \code{hmm_forward}/\code{hmm_backward}, \code{hmm_constrained_likelihood})
#' is untouched here} — this function only prepares inputs for, and
#' post-processes outputs from, that engine (or, if \code{engine =
#' "legacy_canoes"}, the equally-untouched faithful port of the original
#' CANOES HMM core in \code{canoes_legacy_engine.R}).
#'
#' @param sample_name Character. Column name of the test sample in \code{counts}.
#' @param counts Data frame: \code{target, chromosome, start, end, gc}, then one
#'   column per sample.
#' @param p HMM per-target switch probability.
#' @param Tnum Expected number of targets per CNV.
#' @param D Expected genomic span of a CNV (bp). Default 100000.
#' @param numrefs Max number of reference samples to use.
#' @param get_dfs Logical. If TRUE, return emission probs + distances only (used internally).
#' @param homdel_mean Expected read count under homozygous deletion.
#' @param refsample_names Optional character vector restricting the reference pool.
#' @param full_output Logical. If TRUE, return full diagnostic output list.
#' @param min_cor Optional minimum Spearman correlation for reference samples.
#' @param decode_method \code{"distance"} (default) or \code{"stationary"}. Ignored
#'   (forced to \code{"distance"}, with a warning) when \code{engine =
#'   "legacy_canoes"} — the original CANOES HMM never had a stationary mode.
#' @param max_genotype Integer. Skip genotyping if more than this many CNVs are called
#'   (safety valve against runaway computation on pathological samples).
#' @param engine \code{"new"} (default) uses the modern engine in
#'   \code{hmm_engine.R}. \code{"legacy_canoes"} uses a faithful, numerically
#'   verified port of the original published CANOES HMM core instead (see
#'   \code{canoes_legacy_engine.R} for exactly what that does and doesn't
#'   change vs. \code{"new"}). Everything upstream of the HMM — reference
#'   selection, NNLS weighting, variance estimation — is identical either way.
#'
#' @return A data frame of called CNVs (or a list if \code{full_output = TRUE}).
#' @export
call_cnvs <- function(
    sample_name,
    counts,
    p,
    Tnum,
    D = 100000,
    numrefs = 30,
    get_dfs = FALSE,
    homdel_mean = 0.2,
    refsample_names = NULL,
    full_output = FALSE,
    min_cor = NULL,
    decode_method = c("distance", "stationary"),
    max_genotype = 200L,
    engine = c("new", "legacy_canoes")
) {
  decode_method <- match.arg(decode_method)
  engine <- match.arg(engine)
  if (engine == "legacy_canoes" && decode_method == "stationary") {
    warning(
      "engine = 'legacy_canoes' has no 'stationary' transition mode (the ",
      "original CANOES HMM is always distance-decayed); ignoring decode_method = 'stationary'.",
      call. = FALSE
    )
    decode_method <- "distance"
  }
  transition_model <- decode_method
  hmm_par <- hmm_params(p, Tnum, D)

  if (!sample_name %in% colnames(counts))
    stop("No column for sample ", sample_name, " in counts matrix")

  required_meta <- c("target", "chromosome", "start", "end", "gc")
  if (length(setdiff(colnames(counts)[1:5], required_meta)) > 0)
    stop("First five columns must be: target, chromosome, start, end, gc")

  # ============================================================
  # 1. FORCE ALL SAMPLE COLUMNS TO NUMERIC (with proper NA detection)
  # ============================================================
  meta_cols <- c("target", "chromosome", "start", "end", "gc", "GENE", "mean")
  sample_names <- setdiff(colnames(counts), meta_cols)
  for (nm in sample_names) {
    raw_vals <- counts[[nm]]
    counts[[nm]] <- suppressWarnings(as.numeric(as.character(raw_vals)))
    lost <- is.na(counts[[nm]]) & !is.na(raw_vals) & nzchar(as.character(raw_vals))
    if (any(lost)) {
      warning(sprintf(
        "Column '%s' had %d non-numeric value(s) coerced to NA",
        nm, sum(lost)
      ), call. = FALSE)
    }
  }
  if ("mean" %in% colnames(counts)) {
    counts$mean <- as.numeric(counts$mean)
  }
  # ============================================================

  if (p <= 0 || Tnum <= 0 || D <= 0 || numrefs <= 0)
    stop("Parameters p, Tnum, D, and numrefs must all be positive")

  chr <- as.character(counts$chromosome)
  if (all(grepl("^chr", chr))) {
    chr <- dplyr::case_when(
      chr == "chrX" ~ "23", chr == "chrY" ~ "24", TRUE ~ sub("^chr", "", chr)
    )
  } else {
    chr <- dplyr::case_when(chr == "X" ~ "23", chr == "Y" ~ "24", TRUE ~ chr)
  }
  counts$chromosome <- suppressWarnings(as.integer(chr))
  unmapped <- is.na(counts$chromosome)
  if (any(unmapped)) {
    message(sprintf(
      "[WARNING] Dropping %d target(s) with unrecognised chromosome label(s) (e.g. chrM); these cannot be placed on the autosome/X/Y numeric scale used by the HMM.",
      sum(unmapped)
    ))
    counts <- counts[!unmapped, , drop = FALSE]
  }
  # The HMM engines (hmm_engine.R / canoes_legacy_engine.R) step through
  # targets with constructs like `for (i in 2:n)`; with n == 1 that becomes
  # `2:1`, R's classic descending-sequence pitfall, which would silently
  # index row 0 instead of erroring. Requiring >= 2 targets here turns that
  # into a clear, immediate error instead.
  if (nrow(counts) < 2)
    stop("At least 2 targets are required after removing unrecognised chromosomes (found ",
         nrow(counts), ")")

  counts <- dplyr::arrange(counts, chromosome, start)

  meta_cols <- c("target", "chromosome", "start", "end", "gc")
  sample_names <- setdiff(colnames(counts), c(meta_cols, "GENE", "mean"))

  col_medians <- vapply(sample_names, function(nm) median(counts[[nm]], na.rm = TRUE), numeric(1))
  global_med <- median(col_medians, na.rm = TRUE)
  for (nm in sample_names) {
    m <- col_medians[nm]
    if (is.finite(m) && m > 0) counts[[nm]] <- round(counts[[nm]] * global_med / m)
  }

  log_mat <- log2(as.matrix(counts[, sample_names, drop = FALSE]) + 0.5)
  log_mat <- sweep(log_mat, 2, apply(log_mat, 2, median, na.rm = TRUE), "-")
  cor_mat <- stats::cor(log_mat, method = "spearman", use = "pairwise.complete.obs")

  if (is.null(refsample_names) || length(refsample_names) == 0) {
    ref_pool <- setdiff(sample_names, sample_name)
  } else {
    ref_pool <- intersect(refsample_names, sample_names)
  }
  if (length(ref_pool) == 0) stop("No valid reference samples")

  cors <- cor_mat[sample_name, ref_pool]
  if (is.null(names(cors))) names(cors) <- ref_pool

  if (!is.null(min_cor)) {
    pass <- !is.na(cors) & cors >= min_cor
    if (sum(pass) >= 3) cors <- cors[pass] else {
      warning(sprintf(
        "Too few references passed min_cor=%.3f for %s; using top-%d.",
        min_cor, sample_name, numrefs
      ), immediate. = TRUE)
    }
  }

  reference_samples <- names(head(sort(cors, decreasing = TRUE), min(numrefs, length(cors))))
  if (length(reference_samples) < 3)
    stop("Too few valid reference samples for ", sample_name)

  # ============================================================
  # 2. NORMALISE REFERENCE SAMPLES – WITH SAFETY
  # ============================================================
  samp_med <- median(counts[[sample_name]], na.rm = TRUE)
  for (nm in reference_samples) {
    ref_m <- median(counts[[nm]], na.rm = TRUE)
    if (is.finite(ref_m) && ref_m > 0) {
      counts[[nm]] <- as.numeric(counts[[nm]])
      counts[[nm]] <- round(counts[[nm]] * samp_med / ref_m)
    }
  }

  b <- as.numeric(counts[[sample_name]])
  # Force reference columns to be numeric before creating matrix A
  for (nm in reference_samples) {
    counts[[nm]] <- as.numeric(counts[[nm]])
  }
  A <- as.matrix(counts[, reference_samples, drop = FALSE])
  storage.mode(A) <- "numeric"  # ensure matrix is numeric

  # ============================================================
  # 3. NNLS WEIGHTING
  # ============================================================
  set.seed(1L)
  boot_weights <- matrix(0, nrow = 50L, ncol = length(reference_samples))
  for (i in seq_len(50L)) {
    idx <- sample(nrow(A), min(500L, nrow(A)))
    boot_weights[i, ] <- nnls::nnls(A[idx, , drop = FALSE], b[idx])$x
  }
  weights <- colMeans(boot_weights)
  weights[weights < 0] <- 0
  sample_weights <- if (sum(weights) > 0) weights / sum(weights) else
    rep(1 / length(weights), length(weights))

  # ============================================================
  # 4. COMPUTE MEAN AND FILTER
  # ============================================================
  counts$mean <- apply(
    counts[, reference_samples, drop = FALSE], 1,
    function(x) matrixStats::weightedMedian(x, w = sample_weights, na.rm = TRUE)
  )
  counts$mean <- as.numeric(counts$mean)

  keep <- counts$mean >= 10 & is.finite(counts$mean) & counts[[sample_name]] >= 5
  counts <- counts[keep, , drop = FALSE]
  # See the matching guard above -- the HMM engines require at least 2
  # targets to step through safely.
  if (nrow(counts) < 2)
    stop("At least 2 valid targets are required after coverage filtering for ",
         sample_name, " (found ", nrow(counts), ")")
  rownames(counts) <- NULL

  distances <- if (engine == "legacy_canoes") legacy_get_distances(counts) else get_distances(counts)
  ref_data <- counts[, reference_samples, drop = FALSE]
  robust_var <- apply(ref_data, 1, function(x) {
    v <- stats::mad(x, constant = 1.4826, na.rm = TRUE)^2
    if (!is.finite(v) || v <= 0) v <- stats::var(x, na.rm = TRUE)
    v
  })
  robust_var <- pmax(robust_var, counts$mean + 1)

  test_counts <- pmax(round(as.numeric(counts[[sample_name]])), 0L)
  ref_means   <- pmax(round(as.numeric(counts$mean)), 1L)

  em_probs <- if (engine == "legacy_canoes") {
    legacy_emission_probs(test_counts, ref_means, robust_var, counts$target)
  } else {
    emission_probs(test_counts, ref_means, robust_var, counts$target)
  }

  if (get_dfs) return(list(emission_probs = em_probs, distances = distances))

  if (engine == "legacy_canoes") {
    viterbi_df <- legacy_viterbi(em_probs, distances, p, Tnum, D)
  } else {
    viterbi_df <- decode_hmm_states(
      em_probs, distances, hmm_par, transition_model = transition_model
    )
  }

  cnvs <- print_cnvs(sample_name, viterbi_df, counts, reference_samples)

  if (nrow(cnvs) > 0 && nrow(cnvs) <= max_genotype) {
    if (engine == "legacy_canoes") {
      fwd_m <- legacy_forward_matrix(em_probs, distances, p, Tnum, D)
      bwd_m <- legacy_backward_matrix(em_probs, distances, p, Tnum, D)
    } else {
      fwd_m <- hmm_forward(em_probs, distances, hmm_par)
      bwd_m <- hmm_backward(em_probs, distances, hmm_par)
    }
    qualities <- genotype_cnvs(
      cnvs, sample_name, counts, p, Tnum, D, numrefs,
      emission_probs = em_probs, distances = distances,
      forward_matrix = fwd_m, backward_matrix = bwd_m, min_cor = min_cor,
      engine = engine
    )
    for (i in seq_len(nrow(cnvs))) {
      cnvs$Q_SOME[i] <- if (cnvs$CNV[i] == "DEL") qualities[i, "SQDel"] else qualities[i, "SQDup"]
    }
  }

  if (nrow(cnvs) > 0) {
    data_for_cn <- data.frame(
      target = counts$target,
      countsmean = counts$mean,
      varestimate = robust_var,
      sample = test_counts,
      chromosome = counts$chromosome,
      start = counts$start,
      end = counts$end
    )
    cnvs <- calc_copy_number(data_for_cn, cnvs, homdel_mean)
  }

  if (full_output) {
    return(list(
      cnvs = cnvs,
      reference_samples = reference_samples,
      targets = counts$target,
      mean = counts$mean,
      var_estimate = robust_var,
      test_counts = test_counts,
      ref_matrix = as.matrix(ref_data),
      sample_weights = sample_weights,
      emission_probs = em_probs,
      distances = distances,
      engine = engine
    ))
  }
  cnvs
}


#' Genotype Quality Scores for Called CNVs
#'
#' Computes Phred-scaled "no-call" (NQ) and "some-call" (SQ) quality scores
#' for each CNV using the HMM forward/backward matrices. Delegates all
#' likelihood math to \code{hmm_constrained_likelihood()} in
#' \code{hmm_engine.R} (or, under \code{engine = "legacy_canoes"}, to
#' \code{legacy_modified_likelihood()} in \code{canoes_legacy_engine.R} —
#' both files are unmodified by this function either way).
#'
#' @inheritParams call_cnvs
#' @param xcnvs Data frame of called CNVs (the \code{xcnv}-style table).
#' @param emission_probs Optional precomputed emission matrix (avoids recompute).
#' @param distances Optional precomputed distances data frame.
#' @param forward_matrix Optional precomputed HMM forward matrix.
#' @param backward_matrix Optional precomputed HMM backward matrix.
#'
#' @return Data frame: \code{INTERVAL, NQDel, SQDel, NQDup, SQDup}.
#' @export
genotype_cnvs <- function(
    xcnvs, sample_name, counts, p, Tnum, D, numrefs,
    emission_probs = NULL, distances = NULL,
    forward_matrix = NULL, backward_matrix = NULL,
    min_cor = NULL,
    engine = c("new", "legacy_canoes")
) {
  engine <- match.arg(engine)
  if (!sample_name %in% colnames(counts))
    stop("No column for sample ", sample_name)

  if (is.null(emission_probs) || is.null(distances)) {
    l <- call_cnvs(sample_name, counts, p, Tnum, D, numrefs,
                   get_dfs = TRUE, homdel_mean = 0.2, min_cor = min_cor,
                   engine = engine)
    emission_probs <- l$emission_probs
    distances <- l$distances
  }

  if (engine == "legacy_canoes") {
    fwd_m <- if (!is.null(forward_matrix)) forward_matrix else
      legacy_forward_matrix(emission_probs, distances, p, Tnum, D)
    bwd_m <- if (!is.null(backward_matrix)) backward_matrix else
      legacy_backward_matrix(emission_probs, distances, p, Tnum, D)
  } else {
    hmm_par <- hmm_params(p, Tnum, D)
    fwd_m <- if (!is.null(forward_matrix)) forward_matrix else
      hmm_forward(emission_probs, distances, hmm_par)
    bwd_m <- if (!is.null(backward_matrix)) backward_matrix else
      hmm_backward(emission_probs, distances, hmm_par)
  }

  num_cnvs <- nrow(xcnvs)
  qualities <- data.frame(
    INTERVAL = as.character(xcnvs$INTERVAL),
    NQDel = integer(num_cnvs),
    SQDel = integer(num_cnvs),
    NQDup = integer(num_cnvs),
    SQDup = integer(num_cnvs),
    stringsAsFactors = FALSE
  )

  # Phred scoring: engine = "new" keeps the numerical floor added on top of
  # the original formula (see call_cnvs.R history); engine = "legacy_canoes"
  # uses the literal, unguarded original formula (legacy_phred()) on purpose.
  phred <- if (engine == "legacy_canoes") {
    legacy_phred
  } else {
    function(prob) {
      if (is.na(prob) || !is.finite(prob)) return(NA_integer_)
      round(min(99, -10 * log10(max(1 - prob, 1e-10))))
    }
  }

  # Constrained-likelihood dispatch: same c(modified, unmodified) return
  # order from either engine, so the call sites below don't need to change.
  get_lik <- function(forbidden_new, forbidden_legacy, start_target, end_target) {
    if (engine == "legacy_canoes") {
      legacy_modified_likelihood(
        fwd_m, bwd_m, emission_probs, distances,
        start_target, end_target, forbidden_legacy, p, Tnum, D
      )
    } else {
      hmm_constrained_likelihood(
        fwd_m, bwd_m, emission_probs, distances,
        start_target, end_target, forbidden_new, hmm_par
      )
    }
  }

  for (i in seq_len(num_cnvs)) {
    targets <- parse_canope_targets(xcnvs$TARGETS[i])
    if (length(targets) < 2) {
      warning(sprintf("Skipping malformed TARGETS: %s", xcnvs$TARGETS[i]))
      next
    }
    start_target <- min(targets)
    end_target <- max(targets)

    lik_all <- get_lik(
      c(state_deletion, state_duplication),
      c(legacy_state_deletion, legacy_state_duplication),
      start_target, end_target
    )
    prob_all_normal <- exp(lik_all[1] - lik_all[2])

    lik_no_del <- get_lik(state_deletion, legacy_state_deletion, start_target, end_target)
    prob_no_deletion <- exp(lik_no_del[1] - lik_no_del[2])

    lik_no_dup <- get_lik(state_duplication, legacy_state_duplication, start_target, end_target)
    prob_no_duplication <- exp(lik_no_dup[1] - lik_no_dup[2])

    qualities$NQDel[i] <- phred(prob_no_deletion)
    qualities$SQDel[i] <- phred(prob_no_duplication - prob_all_normal)
    qualities$NQDup[i] <- phred(prob_no_duplication)
    qualities$SQDup[i] <- phred(prob_no_deletion - prob_all_normal)
  }
  qualities
}


#' Summarise a Run of Same-State Targets into One CNV Record
#' @noRd
summarize_cnvs <- function(cnv_targets, counts, sample_name, state) {
  cnv_type <- if (state == state_duplication) "DUP" else "DEL"
  row_idx <- targets_to_rows(cnv_targets$target, counts$target)
  if (any(is.na(row_idx))) stop("CNV targets not found in counts table")
  rows <- counts[row_idx, , drop = FALSE]

  start_base <- min(rows$start)
  end_base <- max(rows$end)
  start_target <- min(rows$target)
  end_target <- max(rows$target)
  chromosome <- format_chr_label(rows$chromosome[1])

  gene_str <- NA_character_
  if ("GENE" %in% colnames(counts)) {
    genes <- unique(rows$GENE)
    genes <- genes[!is.na(genes) & genes != "."]
    if (length(genes) > 0) gene_str <- paste(genes, collapse = ";")
  }

  data.frame(
    sample_name = sample_name,
    cnv_type = cnv_type,
    cnv_interval = paste0(chromosome, ":", start_base, "-", end_base),
    cnv_kbs = (end_base - start_base) / 1000,
    cnv_chromosome = chromosome,
    cnv_midbp = round((end_base + start_base) / 2),
    cnv_targets = paste0(start_target, "..", end_target),
    num_targets = end_target - start_target + 1,
    cnv_gene = gene_str,
    stringsAsFactors = FALSE
  )
}


#' Convert Viterbi State Path into CNV Call Records
#'
#' @noRd
print_cnvs <- function(test_sample_name, viterbi_df, nonzero_counts, reference_samples) {
  consecutive_groups <- function(idx, chrom) {
    g <- integer(length(idx))
    g[1] <- 1L
    if (length(idx) > 1) {
      for (i in 2:length(idx)) {
        same_chrom  <- identical(chrom[i], chrom[i - 1])
        contiguous  <- idx[i] == idx[i - 1] + 1L
        g[i] <- g[i - 1] + as.integer(!(same_chrom && contiguous))
      }
    }
    g
  }

  results <- list()
  for (state in c(state_deletion, state_duplication)) {
    row_idx <- which(viterbi_df$viterbi_state == state)
    if (length(row_idx) > 0) {
      chrom_at_idx <- nonzero_counts$chromosome[row_idx]
      df <- data.frame(
        target = viterbi_df$target[row_idx],
        group = consecutive_groups(row_idx, chrom_at_idx)
      )
      results <- c(results, lapply(split(df, df$group), function(g)
        summarize_cnvs(g, nonzero_counts, test_sample_name, state)))
    }
  }

  empty_df <- data.frame(
    SAMPLE = character(0), CNV = character(0), INTERVAL = character(0),
    KB = numeric(0), CHR = character(0), MID_BP = numeric(0),
    TARGETS = character(0), NUM_TARG = integer(0), GENE = character(0),
    MLCN = integer(0), Q_SOME = integer(0),
    NUM_REFS = integer(0), REF_SAMPLES = character(0),
    stringsAsFactors = FALSE
  )
  if (length(results) == 0) return(empty_df)

  cnvs_df <- do.call(rbind, results)
  rownames(cnvs_df) <- NULL

  xcnv <- data.frame(
    SAMPLE = cnvs_df$sample_name,
    CNV = cnvs_df$cnv_type,
    INTERVAL = cnvs_df$cnv_interval,
    KB = cnvs_df$cnv_kbs,
    CHR = cnvs_df$cnv_chromosome,
    MID_BP = cnvs_df$cnv_midbp,
    TARGETS = cnvs_df$cnv_targets,
    NUM_TARG = cnvs_df$num_targets,
    GENE = if ("cnv_gene" %in% colnames(cnvs_df)) cnvs_df$cnv_gene else NA_character_,
    MLCN = 0L,
    Q_SOME = NA_integer_,
    NUM_REFS = length(reference_samples),
    REF_SAMPLES = paste(reference_samples, collapse = ";"),
    stringsAsFactors = FALSE
  )
  cat(nrow(xcnv), "CNVs called in sample", test_sample_name, "\n")
  xcnv
}


#' Estimate Most-Likely Copy Number for Each CNV (negative-binomial states C0-C6)
#' @export
calc_copy_number <- function(data, cnvs, homdel_mean) {
  for (i in seq_len(nrow(cnvs))) {
    cnv_row <- cnvs[i, ]
    tb <- parse_canope_targets(cnv_row$TARGETS)
    if (length(tb) < 2) next

    cd <- data[data$target >= tb[1] & data$target <= tb[2], , drop = FALSE]
    if (nrow(cd) == 0) next

    disp <- cd$varestimate - cd$countsmean
    size <- ifelse(disp > 0, cd$countsmean^2 / disp, Inf)
    state_means <- t(apply(data.frame(x = cd$countsmean), 1,
                           function(x) c(x / 2, x, x * 3 / 2, x * 2, x * 5 / 2, x * 3)))

    em <- matrix(NA_real_, nrow(cd), 7L)
    colnames(em) <- paste0("C", 0:6)
    em[, "C0"] <- stats::dpois(cd$sample, homdel_mean, log = TRUE)
    for (s in 1:6) {
      em[, s + 1] <- stats::dnbinom(
        cd$sample, mu = state_means[, s], size = size * (s / 2), log = TRUE
      )
    }

    ml_state <- which.max(colSums(em, na.rm = TRUE)) - 1L
    if (ml_state == 2L) ml_state <- if (cnv_row$CNV == "DEL") 1L else 3L
    cnvs$MLCN[i] <- ml_state
  }
  cnvs
}