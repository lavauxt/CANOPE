#' Stop if a value is NULL or empty
#' @param val Object to validate.
#' @param msg Error message.
#' @export
stop_if_missing <- function(val, msg) {
  if (is.null(val) || length(val) == 0) stop(msg)
}

#' Stop if a file does not exist
#' @param path File path.
#' @param msg Error message.
#' @export
stop_if_not_file <- function(path, msg) {
  if (is.null(path) || !file.exists(path)) stop(msg)
}

#' Load an RData file into an isolated environment
#' @param path Path to the .RData file.
#' @param required Optional character vector of required object names.
#' @return Named list of objects loaded from the file.
#' @export
canope_load_rdata <- function(path, required = NULL) {
  stop_if_not_file(path, paste0("[ERROR] RData not found: ", path))
  env <- new.env()
  load(path, envir = env)
  if (!is.null(required)) {
    missing_vars <- setdiff(required, ls(env))
    if (length(missing_vars) > 0)
      stop("[ERROR] Missing variables in RData: ", paste(missing_vars, collapse = ", "),
           "\nAvailable: ", paste(ls(env), collapse = ", "))
  }
  as.list(env)
}

#' Null-coalescing operator
#' @param a Primary value.
#' @param b Fallback value.
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Sanitize a string for use as a filename
#' @param name Character string.
#' @export
sanitize_filename <- function(name) gsub("[^[:alnum:]._-]", "_", as.character(name))

#' Filter a data frame by a chromosome column
#'
#' Handles both \code{CHR} and \code{chromosome} column names and normalises
#' chr-prefix mismatches automatically.
#'
#' @param df Data frame with a \code{CHR} or \code{chromosome} column.
#' @param include Optional chromosome names to keep.
#' @param exclude Optional chromosome names to remove.
#' @return Filtered data frame.
#' @export
canope_filter_chromosomes <- function(df, include = NULL, exclude = NULL) {
  if (is.null(df) || nrow(df) == 0) return(df)
  chrom_col <- intersect(c("CHR", "chromosome", "Chromosome"), colnames(df))[1]
  if (is.na(chrom_col)) stop("Data frame must have a CHR / chromosome column.")
  norm <- function(x) unique(c(x, sub("^chr", "", x), paste0("chr", sub("^chr", "", x))))
  if (!is.null(include)) df <- df[df[[chrom_col]] %in% norm(include), , drop = FALSE]
  if (!is.null(exclude)) df <- df[!df[[chrom_col]] %in% norm(exclude), , drop = FALSE]
  df
}

#' Assign sequential exon numbers within each gene
#'
#' Numbers exons 1..n per gene in genomic order (chromosome, then start, then
#' end), regardless of the row order of \code{bed_df}. Duplicate intervals
#' (same chrom + start + end + gene) are \emph{not} counted as separate
#' exons -- they share the \code{exon_number} of their first occurrence --
#' but they are NOT dropped from the output. The returned data frame always
#' has exactly \code{nrow(bed_df)} rows, in the same order as the input; the
#' genomic sort used for numbering is internal only, so callers can safely
#' bind \code{exon_number} back onto \code{bed_df} (or any table with the
#' same row identity, e.g. a matching counts table) purely by position.
#'
#' @param bed_df Data frame with chromosome/Chr, start/Start, end/End, gene/Gene/GENE columns.
#' @return The input data frame (original row order and row count preserved)
#'   with an added \code{exon_number} integer column.
#' @importFrom data.table := as.data.table setnames setorder .N .I
#' @export
assign_exon_numbers_per_gene <- function(bed_df) {
  chrom_col <- intersect(c("chromosome", "Chr", "CHROM"), names(bed_df))[1]
  start_col <- intersect(c("start", "Start", "START"), names(bed_df))[1]
  end_col   <- intersect(c("end", "End", "END"), names(bed_df))[1]
  gene_col  <- intersect(c("GENE", "gene", "Gene"), names(bed_df))[1]
  stopifnot(!is.na(chrom_col), !is.na(start_col), !is.na(end_col), !is.na(gene_col))

  n_in <- nrow(bed_df)

  dt <- data.table::as.data.table(bed_df)
  data.table::setnames(dt, c(chrom_col, start_col, end_col, gene_col),
                       c("._chrom", "._start", "._end", "._gene"))

  # Track original row position so input order can be restored after sorting.
  dt[, ._orig_row := .I]
  dt[, ._key := paste(`._chrom`, `._start`, `._end`, `._gene`, sep = "\r")]

  dup_rows <- duplicated(dt, by = "._key")
  if (any(dup_rows)) {
    warning(sprintf(
      "Found %d duplicate BED row(s) (same chrom, start, end, gene); duplicates will share the exon_number of their first occurrence.",
      sum(dup_rows)), immediate. = TRUE)
  }
  # Number only the unique (chrom, start, end, gene) combinations -- a
  # duplicate row must NOT get its own exon_number bumped from the count,
  # or a gene with a repeated interval would appear to have more exons than
  # it does. Crucially, this must not shrink the *output* row count: every
  # caller relies on binding exon_number back onto bed_df (or any table
  # sharing its row identity) purely by position, so duplicates are matched
  # back onto every original row below rather than dropped.
  dt_unique <- dt[!dup_rows]

  chrom_levels <- c(paste0("chr", c(1:22, "X", "Y", "M")), c(as.character(1:22), "X", "Y", "M"))
  dt_unique[, .chrom_fac := factor(`._chrom`, levels = unique(c(chrom_levels, unique(`._chrom`))))]
  data.table::setorder(dt_unique, .chrom_fac, `._start`, `._end`)
  dt_unique[, .chrom_fac := NULL]

  dt_unique[, exon_number := seq_len(.N), by = "._gene"]

  # Map exon_number back onto every original row (including any duplicates
  # identified above) by key, then restore the original row order. This
  # guarantees nrow(output) == nrow(bed_df) always, regardless of duplicates.
  dt[, exon_number := dt_unique$exon_number[match(`._key`, dt_unique$`._key`)]]

  data.table::setorder(dt, ._orig_row)
  dt[, c("._orig_row", "._key") := NULL]

  data.table::setnames(dt, c("._chrom", "._start", "._end", "._gene"),
                       c(chrom_col, start_col, end_col, gene_col))
  out <- as.data.frame(dt)
  stopifnot(nrow(out) == n_in)
  out
}

#' Compute within-gene exon index from a BED-like data frame
#'
#' Returns the exon_number column if already present, otherwise calls
#' \code{\link{assign_exon_numbers_per_gene}}.
#'
#' @param bed_df Data frame (as above).
#' @return Integer vector, one per row.
#' @export
compute_exon_index <- function(bed_df) {
  if ("exon_number" %in% names(bed_df) && !all(is.na(bed_df$exon_number)))
    return(bed_df$exon_number)
  assign_exon_numbers_per_gene(bed_df)$exon_number
}

#' Reclassify off-target / filler BED intervals before exon numbering
#'
#' Ported from ECHO. Some target panels include intervals that were never a
#' real gene exon at all: normalization/backbone probes placed off-target
#' for coverage calibration, commonly named things like \code{"HorsROI"}
#' ("hors ROI" is French for "outside the region of interest"),
#' \code{"OffTarget"}, \code{"Backbone"}, and so on. Left alone, the BED-name
#' parser extracts whatever the first token of that name happens to be
#' (e.g. \code{"HorsROI"}) as if it were a gene symbol, and
#' \code{\link{assign_exon_numbers_per_gene}} then numbers it 1..n exactly
#' like a real gene with its own exons -- so a plot window that happens to
#' straddle one of these intervals shows it interleaved with the real
#' gene's exons under its own (fake) "gene" tile. This function catches
#' those rows, by name, \strong{before} any numbering happens, and lets the
#' caller choose what should happen to them.
#'
#' @param bed_df data.frame with chromosome/Chr, start/Start, end/End,
#'   GENE/gene/Gene columns (1-based coordinates; genomic order not
#'   required -- \code{handling = "merge"} sorts internally).
#' @param pattern Character. A regular expression (matched against the
#'   gene column, case-sensitively, consistent with the exon-name parser)
#'   identifying off-target/filler rows -- e.g. \code{"^HorsROI"}, or
#'   \code{"^(HorsROI|OffTarget|Backbone)$"} for a panel using several such
#'   labels. \code{NULL} or \code{""} disables this feature entirely
#'   (returns \code{bed_df} unchanged). Default \code{"^HorsROI"}.
#' @param handling One of:
#'   \itemize{
#'     \item \code{"na"} (default) -- keep the interval (it still gets
#'       coverage extracted and still contributes background signal) but
#'       set the gene column to \code{NA} on those rows, so
#'       \code{assign_exon_numbers_per_gene()} and everything downstream
#'       (plots, VCF GENE column, confidence scoring) leaves them
#'       un-numbered and out of any gene-based grouping, instead of
#'       numbering the filler label as if it were a gene.
#'     \item \code{"remove"} -- drop those rows entirely.
#'     \item \code{"merge"} -- reassign the gene column to the nearest
#'       neighbouring \emph{real} gene on the same chromosome (ties go to
#'       the preceding gene), so the interval is treated as one of that
#'       gene's own exons and gets numbered like any other exon. A row
#'       with no real gene anywhere on its chromosome is left unchanged.
#'   }
#' @param verbose Logical. Print a one-line summary. Default \code{TRUE}.
#' @return The same data.frame: column set and row order preserved for
#'   \code{"na"}/\code{"merge"}; row count reduced for \code{"remove"}.
#' @export
handle_off_target_regions <- function(bed_df, pattern = "^HorsROI",
                                      handling = c("na", "remove", "merge"),
                                      verbose = TRUE) {
  handling <- match.arg(handling)
  if (is.null(pattern) || !nzchar(pattern)) return(bed_df)
  if (is.null(bed_df) || nrow(bed_df) == 0) return(bed_df)

  chrom_col <- intersect(c("chromosome", "Chr", "CHROM"), names(bed_df))[1]
  start_col <- intersect(c("start", "Start", "START"), names(bed_df))[1]
  end_col   <- intersect(c("end", "End", "END"), names(bed_df))[1]
  gene_col  <- intersect(c("GENE", "gene", "Gene"), names(bed_df))[1]
  if (is.na(gene_col)) return(bed_df)

  gene_vals <- as.character(bed_df[[gene_col]])
  is_off <- !is.na(gene_vals) & grepl(pattern, gene_vals, perl = TRUE)
  n_off <- sum(is_off)
  if (n_off == 0) return(bed_df)

  if (handling == "remove") {
    out <- bed_df[!is_off, , drop = FALSE]
    if (verbose) {
      message(sprintf("[INFO] handle_off_target_regions: removed %d off-target interval(s) matching /%s/.",
                      n_off, pattern))
    }
    return(out)
  }

  if (handling == "na") {
    bed_df[[gene_col]][is_off] <- NA_character_
    if (verbose) {
      message(sprintf(
        "[INFO] handle_off_target_regions: set gene = NA for %d off-target interval(s) matching /%s/ (kept, excluded from exon numbering).",
        n_off, pattern))
    }
    return(bed_df)
  }

  # handling == "merge": walk the off-target rows in genomic order and
  # attach each one to whichever real gene -- the previous one or the
  # next one, on the same chromosome -- sits closer. A simple forward/
  # backward carry-forward pass (O(n), two linear scans) rather than a
  # per-row search, since a panel BED can run into the tens of thousands
  # of rows.
  ord     <- order(bed_df[[chrom_col]], bed_df[[start_col]], bed_df[[end_col]])
  n       <- length(ord)
  chrom_s <- as.character(bed_df[[chrom_col]])[ord]
  start_s <- bed_df[[start_col]][ord]
  end_s   <- bed_df[[end_col]][ord]
  gene_s  <- gene_vals[ord]
  off_s   <- is_off[ord]

  prev_gene <- character(n); prev_end <- numeric(n)
  g <- NA_character_; e <- NA_real_; last_chrom <- NA_character_
  for (i in seq_len(n)) {
    if (!identical(chrom_s[i], last_chrom)) { g <- NA_character_; e <- NA_real_; last_chrom <- chrom_s[i] }
    prev_gene[i] <- g; prev_end[i] <- e
    if (!off_s[i]) { g <- gene_s[i]; e <- end_s[i] }
  }

  next_gene <- character(n); next_start <- numeric(n)
  g <- NA_character_; s <- NA_real_; last_chrom <- NA_character_
  for (i in rev(seq_len(n))) {
    if (!identical(chrom_s[i], last_chrom)) { g <- NA_character_; s <- NA_real_; last_chrom <- chrom_s[i] }
    next_gene[i] <- g; next_start[i] <- s
    if (!off_s[i]) { g <- gene_s[i]; s <- start_s[i] }
  }

  new_gene_s <- gene_s
  n_merged   <- 0L
  for (i in which(off_s)) {
    has_prev <- !is.na(prev_gene[i])
    has_next <- !is.na(next_gene[i])
    if (has_prev && has_next) {
      d_prev <- start_s[i] - prev_end[i]
      d_next <- next_start[i] - end_s[i]
      new_gene_s[i] <- if (d_prev <= d_next) prev_gene[i] else next_gene[i]
      n_merged <- n_merged + 1L
    } else if (has_prev) {
      new_gene_s[i] <- prev_gene[i]; n_merged <- n_merged + 1L
    } else if (has_next) {
      new_gene_s[i] <- next_gene[i]; n_merged <- n_merged + 1L
    } # else: no real gene anywhere on this chromosome -- leave as-is
  }

  bed_df[[gene_col]][ord] <- new_gene_s
  if (verbose) {
    message(sprintf(
      "[INFO] handle_off_target_regions: merged %d/%d off-target interval(s) matching /%s/ into their nearest neighbouring gene (%d had no real gene on their chromosome to attach to).",
      n_merged, n_off, pattern, n_off - n_merged))
  }
  bed_df
}

#' Pad the outer edge of each gene's first and last exon
#'
#' Ported from ECHO. Capture-based coverage often drops off right at the
#' true edge of a target interval (probe/bait tiling is rarely perfect
#' exactly at the boundary, and reads whose alignment barely spans the edge
#' get soft-clipped or excluded). For an internal exon this is usually
#' harmless -- its neighbours carry the signal -- but for a gene's *first*
#' (lowest-coordinate) or *last* (highest-coordinate) exon there is no such
#' neighbour on the outward side, so a thin sliver of low/zero coverage
#' right at that edge can pull the whole exon's count down. This function
#' extends only the outward-facing edge of those two terminal exons per
#' gene (both edges, for a single-exon gene) by \code{padding} bases.
#'
#' "First"/"last" follows the same purely-genomic (chrom, start, end)
#' ordering that \code{\link{assign_exon_numbers_per_gene}} uses everywhere
#' else in CANOPE. Internal exons (and the inward-facing edge of a
#' terminal exon) are left untouched.
#'
#' Padding is applied on a best-effort basis ("if possible"): the function
#' never creates an overlap with whatever interval sits next to it on the
#' same chromosome (a neighbouring exon of the same gene or of a different
#' one), and never pushes a coordinate below 1 or past the contig length
#' (when \code{chr_lengths} is supplied). Where the available gap is
#' narrower than \code{padding} -- including the case where two different
#' genes' terminal exons sit right next to each other and both want a
#' share of the same small gap -- that gap is split between the two
#' competing sides rather than let either one overlap the other.
#'
#' @param bed_df data.frame with columns chromosome/Chr, start/Start,
#'   end/End, GENE/gene/Gene (1-based, inclusive coordinates).
#' @param padding Integer >= 0. Bases to add to the outward edge of each
#'   gene's first and last exon. \code{0} (the default) disables padding
#'   and returns \code{bed_df} unchanged.
#' @param chr_lengths Optional named numeric vector (names = chromosome,
#'   values = contig length) used to cap the last exon's End at the contig
#'   boundary. If \code{NULL}, no contig-length clamp is applied (only the
#'   neighbouring-interval clamp).
#' @param verbose Logical. Print a one-line summary. Default \code{TRUE}.
#' @return The same data.frame (original row order and row count
#'   preserved), with Start/End adjusted for terminal-exon rows only.
#' @export
pad_gene_terminal_exons <- function(bed_df, padding = 0, chr_lengths = NULL, verbose = TRUE) {
  if (is.null(padding) || length(padding) != 1 || is.na(padding) || padding <= 0) {
    return(bed_df)
  }
  padding <- as.integer(round(padding))

  chrom_col <- intersect(c("chromosome", "Chr", "CHROM"), names(bed_df))[1]
  start_col <- intersect(c("start", "Start", "START"), names(bed_df))[1]
  end_col   <- intersect(c("end", "End", "END"), names(bed_df))[1]
  gene_col  <- intersect(c("GENE", "gene", "Gene"), names(bed_df))[1]
  stopifnot(!is.na(chrom_col), !is.na(start_col), !is.na(end_col), !is.na(gene_col))

  n_in <- nrow(bed_df)
  if (n_in == 0) return(bed_df)

  # Reuse the pipeline's own per-gene ordering so "first"/"last" here
  # always agrees with exon_number everywhere else in CANOPE.
  numbered <- assign_exon_numbers_per_gene(bed_df)

  dt <- data.table::as.data.table(numbered)
  dt[, .orig_row := .I]
  data.table::setnames(dt, c(chrom_col, start_col, end_col, gene_col),
                       c("chrom", "start", "end", "gene"))

  dt[, is_first := exon_number == 1L]
  dt[, is_last  := exon_number == max(exon_number), by = "gene"]
  no_gene <- is.na(dt$gene) | dt$gene %in% c("", ".", "Unknown")
  dt[no_gene, c("is_first", "is_last") := FALSE]

  # Sort a copy by genomic position (per chromosome) so each terminal
  # exon can see its nearest neighbour on either side -- regardless of
  # which gene that neighbour belongs to -- and never be padded into it.
  chrom_levels <- c(paste0("chr", c(1:22, "X", "Y", "M")),
                    c(as.character(1:22), "X", "Y", "M"))
  dt[, .chrom_fac := factor(chrom, levels = unique(c(chrom_levels, unique(chrom))))]
  data.table::setorder(dt, .chrom_fac, start, end)

  n        <- nrow(dt)
  chrom_id <- as.integer(dt$.chrom_fac)
  start_v  <- dt$start
  end_v    <- dt$end
  is_first_v <- dt$is_first
  is_last_v  <- dt$is_last

  right_extend <- integer(n)  # applies to is_last rows: bp added to end
  left_extend  <- integer(n)  # applies to is_first rows: bp subtracted from start

  # Gap k (k = 1..n-1) sits between sorted row k and row k+1. Both may
  # want a share of it at once -- row k if it's a last exon growing
  # rightward, row k+1 if it's a first exon growing leftward (this is
  # the one place two *different* genes' terminal exons can compete for
  # the same free space). Give each what it asks for if the gap is big
  # enough for both; otherwise split the gap between them so neither
  # padded interval ever crosses into the other's.
  if (n > 1) {
    same_chr_pair <- chrom_id[-n] == chrom_id[-1]
    gap        <- pmax(start_v[-1] - end_v[-n] - 1L, 0L)
    want_left  <- ifelse(is_last_v[-n],  padding, 0L)  # row k wants to grow right
    want_right <- ifelse(is_first_v[-1], padding, 0L)  # row k+1 wants to grow left
    demand     <- want_left + want_right

    grant_left  <- integer(n - 1L)
    grant_right <- integer(n - 1L)
    has_demand  <- same_chr_pair & demand > 0
    fits        <- has_demand & demand <= gap
    tight       <- has_demand & demand > gap

    grant_left[fits]  <- want_left[fits]
    grant_right[fits] <- want_right[fits]
    grant_left[tight]  <- as.integer(floor(gap[tight] * want_left[tight] / demand[tight]))
    grant_right[tight] <- gap[tight] - grant_left[tight]

    right_extend[-n] <- grant_left
    left_extend[-1]  <- grant_right
  }

  # Rows at a chromosome boundary (no same-chromosome neighbour on the
  # relevant side) have no interval to compete with there, so they fall
  # back to the contig start (position 1) / contig length instead.
  has_prev <- c(FALSE, if (n > 1) chrom_id[-1] == chrom_id[-n] else logical(0))
  has_next <- c(if (n > 1) chrom_id[-n] == chrom_id[-1] else logical(0), FALSE)

  no_prev_first <- is_first_v & !has_prev
  if (any(no_prev_first)) {
    left_extend[no_prev_first] <- pmin(padding, pmax(start_v[no_prev_first] - 1L, 0L))
  }

  no_next_last <- is_last_v & !has_next
  if (any(no_next_last)) {
    chr_len_here <- if (!is.null(chr_lengths)) {
      unname(chr_lengths[as.character(dt$chrom[no_next_last])])
    } else {
      rep(NA_real_, sum(no_next_last))
    }
    avail <- ifelse(!is.na(chr_len_here), chr_len_here - end_v[no_next_last], Inf)
    right_extend[no_next_last] <- pmin(padding, pmax(avail, 0))
  }

  n_padded_start <- sum(left_extend > 0L)
  n_padded_end   <- sum(right_extend > 0L)
  n_clamped      <- sum(is_first_v & left_extend  < padding) +
                     sum(is_last_v  & right_extend < padding)

  dt[, start := start_v - left_extend]
  dt[, end   := end_v   + right_extend]

  data.table::setorder(dt, .orig_row)  # restore original (input) row order
  dt[, c(".orig_row", ".chrom_fac", "is_first", "is_last", "exon_number") := NULL]
  data.table::setnames(dt, c("chrom", "start", "end", "gene"),
                       c(chrom_col, start_col, end_col, gene_col))
  out <- as.data.frame(dt)
  stopifnot(nrow(out) == n_in)

  if (verbose) {
    message(sprintf(
      "[INFO] pad_gene_terminal_exons: requested %d bp padding -- extended %d gene start(s) and %d gene end(s); %d side(s) received less than the full request (shared gap with a neighbouring interval, or a contig boundary).",
      padding, n_padded_start, n_padded_end, n_clamped))
  }
  out
}

#' Pad a BED \emph{file}'s gene-terminal exons (file-level wrapper)
#'
#' CANOPE's coverage/GC-content functions (\code{\link{get_coverage_from_bams}},
#' \code{\link{get_coverage_from_bams_megadepth}}, \code{\link{compute_gc_from_fasta}},
#' \code{\link{compute_gc_from_bed}}) all take a BED \emph{path}, re-reading it
#' from disk independently -- unlike ECHO, where a single in-memory
#' \code{bed_file} is padded once and reused for everything. This wrapper
#' reads a BED file, applies \code{\link{pad_gene_terminal_exons}}, and
#' writes the result back out to \code{output_bed} in the same column
#' layout as the input, so \code{run_canope()} can point every downstream
#' consumer at one consistently-padded file.
#'
#' Internally converts the (0-based, half-open) BED coordinates to 1-based
#' inclusive before padding, and back to 0-based on write-out, so the
#' padding math is identical to ECHO's (which operates on 1-based BED
#' output from \code{process_bed_file()}).
#'
#' @param input_bed Path to a BED file (>= 3 columns; a 4th "name"/gene
#'   column is required for per-gene padding to be meaningful -- without
#'   one, every row is treated as its own single-exon "gene").
#' @param output_bed Path to write the padded BED to.
#' @param padding Integer >= 0 bp of padding (see \code{\link{pad_gene_terminal_exons}}).
#'   \code{0} just copies \code{input_bed} to \code{output_bed} unchanged.
#' @param chr_lengths Optional named numeric vector of contig lengths (see
#'   \code{\link{pad_gene_terminal_exons}}).
#' @param verbose Logical. Default \code{TRUE}.
#' @return Invisibly returns \code{output_bed}.
#' @export
pad_bed_file <- function(input_bed, output_bed, padding = 0, chr_lengths = NULL, verbose = TRUE) {
  if (is.null(padding) || is.na(padding) || padding <= 0) {
    if (!identical(normalizePath(input_bed, mustWork = FALSE),
                   normalizePath(output_bed, mustWork = FALSE))) {
      file.copy(input_bed, output_bed, overwrite = TRUE)
    }
    return(invisible(output_bed))
  }

  raw <- utils::read.table(input_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (ncol(raw) < 3) stop("[ERROR] pad_bed_file: input BED must have at least 3 columns: ", input_bed)
  colnames(raw)[1:3] <- c("chromosome", "start", "end")
  if (ncol(raw) >= 4) colnames(raw)[4] <- "GENE" else raw$GENE <- paste0("Target_", seq_len(nrow(raw)))

  # 0-based half-open -> 1-based inclusive, so this reuses the exact same
  # (already-verified) gap math as ECHO's 1-based pipeline.
  raw$start <- raw$start + 1L

  padded <- pad_gene_terminal_exons(raw, padding = padding, chr_lengths = chr_lengths, verbose = verbose)

  # 1-based inclusive -> back to 0-based half-open for BED output.
  padded$start <- padded$start - 1L

  # Keep every original column (including any 5th+ columns beyond
  # chrom/start/end/gene, e.g. score/strand), but drop the synthetic GENE
  # column again if the input never had a 4th column to begin with.
  keep   <- if (ncol(raw) >= 4) colnames(raw) else setdiff(colnames(raw), "GENE")
  out_df <- padded[, keep, drop = FALSE]

  dir.create(dirname(output_bed), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(out_df, file = output_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  if (verbose) message("[INFO] pad_bed_file: padded BED written: ", output_bed)
  invisible(output_bed)
}

#' Compute gap-inserted x-axis positions for CNV window plots
#'
#' Ported from ECHO. All of CANOPE's per-call plots (the four PDF panels in
#' \code{generate_plots.R} and the three interactive panels in
#' \code{CANOPE_report.Rmd}) lay a window of exons out along a single
#' x-axis. Plotted at plain 1..n integer positions, a gene boundary inside
#' that window looks identical to an ordinary intron between two exons of
#' the *same* gene -- there's nothing to tell a reader "these two points
#' belong to different genes" other than the tile-track colour (PDF only;
#' the HTML report has no tile track at all). This function computes an
#' alternative x-position (\code{px}) for each exon in the window that
#' inserts \code{gap} extra, unlabelled axis units wherever the gene
#' column changes between consecutive exons -- i.e. between a gene's last
#' exon and the next gene's first exon -- while keeping ordinary
#' within-gene spacing at a plain 1 unit.
#'
#' It also returns \code{gene_group}, a per-position integer that
#' increments at every such boundary. Passing this as the \code{group}
#' aesthetic on a \code{geom_line()}/\code{geom_ribbon()} keeps that visual
#' gap genuinely blank -- otherwise ggplot draws a single connected
#' line/ribbon straight across it.
#'
#' @param bed_df data.frame with a GENE/gene/Gene column. \code{exon_range}
#'   values are row indices into this data.frame.
#' @param exon_range Integer vector of \code{bed_df} row indices, in the
#'   order they'll be plotted along the x-axis (ascending genomic order).
#' @param gap Numeric >= 0. Extra x-axis units inserted at each gene
#'   boundary. \code{0} falls back to plain 1..n spacing (no visual gap,
#'   but \code{gene_group} is still computed correctly). Default \code{1}.
#' @return data.frame with one row per element of \code{exon_range}:
#'   \code{idx}, \code{px}, \code{gene_break}, \code{gene_group}.
#' @export
compute_gene_gap_positions <- function(bed_df, exon_range, gap = 1) {
  if (length(exon_range) == 0) {
    return(data.frame(idx = integer(0), px = numeric(0),
                      gene_break = logical(0), gene_group = integer(0)))
  }
  gap <- if (is.null(gap) || is.na(gap) || gap < 0) 1 else gap

  gene_col <- intersect(c("GENE", "gene", "Gene"), names(bed_df))[1]
  genes    <- as.character(bed_df[[gene_col]][exon_range])
  genes[is.na(genes)] <- ""

  gene_break <- c(FALSE, genes[-1] != genes[-length(genes)])
  step       <- ifelse(gene_break, 1 + gap, 1)
  step[1]    <- 0
  px         <- cumsum(step) + 1

  data.frame(idx = exon_range, px = px, gene_break = gene_break,
            gene_group = cumsum(gene_break) + 1L, stringsAsFactors = FALSE)
}

#' Normalise chromosome names to match a reference naming style (vectorised)
#'
#' @param chr_vec Character vector of chromosome names to normalise.
#' @param ref_chromosomes Reference vector that defines which style to adopt.
#' @return Character vector with chr-prefix aligned to the reference style.
#' @export
normalize_chromosome_vec <- function(chr_vec, ref_chromosomes) {
  has_chr <- any(grepl("^chr", ref_chromosomes))
  if (has_chr) ifelse(grepl("^chr", chr_vec), chr_vec, paste0("chr", chr_vec))
  else          ifelse(grepl("^chr", chr_vec), sub("^chr", "", chr_vec), chr_vec)
}

#' Check whether a BAM file has an index on disk
#' @param bam Path to BAM file.
#' @return Logical.
#' @noRd
bam_has_index <- function(bam) {
  file.exists(paste0(bam, ".bai")) || file.exists(paste0(bam, ".bam.bai"))
}

#' Parse a CANOPE INTERVAL string into components
#'
#' @param interval Character string of the form \code{"chr1:100-200"}.
#' @return Named list with \code{chrom}, \code{start}, \code{end}.
#' @export
parse_canope_interval <- function(interval) {
  interval <- as.character(interval)
  colon    <- regexpr(":", interval, fixed = TRUE)
  chrom    <- substr(interval, 1L, colon - 1L)
  rest     <- substr(interval, colon + 1L, nchar(interval))
  dash     <- regexpr("-", rest, fixed = TRUE)
  list(
    chrom = chrom,
    start = as.integer(substr(rest, 1L, dash - 1L)),
    end   = as.integer(substr(rest, dash + 1L, nchar(rest)))
  )
}

#' Parse a CANOPE TARGETS string into target ID range
#'
#' @param targets Character string of the form \code{"42..47"}.
#' @return Integer vector \code{c(start_target, end_target)}.
#' @export
parse_canope_targets <- function(targets) {
  parts <- as.integer(strsplit(as.character(targets), "..", fixed = TRUE)[[1]])
  parts[!is.na(parts)]
}

# ─── PCA of coverage profiles ──────────────────────────────────────────────

#' Plot PCA of Sample Coverage Profiles (CANOPE version)
#'
#' Mirrors ECHO's \code{plot_coverage_pca} but uses CANOPE's column layout
#' where sample columns begin after the metadata columns.
#'
#' @param counts      Data frame (rows = targets, columns = metadata + samples).
#' @param sample_names Character vector of sample column names.
#' @param output_pdf  Optional path for PDF output (NULL for on-screen).
#' @param color_by    Optional factor/character vector of group labels per sample.
#' @param scale       Logical; scale variables before PCA. Default TRUE.
#' @return Invisibly returns list(pca, var_exp).
#' @export
plot_coverage_pca <- function(counts, sample_names,
                              output_pdf = NULL, color_by = NULL, scale = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for PCA plot.")

  count_mat <- as.matrix(counts[, sample_names, drop = FALSE])
  log_mat   <- log2(count_mat + 1)
  row_var   <- apply(log_mat, 1, var, na.rm = TRUE)
  keep      <- row_var > 0 & !is.na(row_var)
  if (sum(keep) < 2) {
    warning("Insufficient variation for PCA (< 2 informative rows).")
    return(invisible(NULL))
  }
  log_mat <- log_mat[keep, , drop = FALSE]

  pca     <- stats::prcomp(t(log_mat), scale. = scale, center = TRUE)
  var_exp <- summary(pca)$importance[2, ] * 100
  pca_df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         Sample = sample_names, stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = Sample))

  if (!is.null(color_by) && length(color_by) == nrow(pca_df)) {
    pca_df$Group <- as.factor(color_by)
    p <- p + ggplot2::aes(colour = Group)
    if (length(unique(color_by)) >= 2)
      p <- p + ggplot2::stat_ellipse(ggplot2::aes(colour = Group),
                                      type = "norm", linetype = "dashed")
  }

  p <- p +
    ggplot2::geom_point(size = 3,
                        colour = if (is.null(color_by)) "steelblue" else NULL) +
    ggplot2::labs(
      x     = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y     = paste0("PC2 (", round(var_exp[2], 1), "%)"),
      title = "PCA – CANOPE Coverage Profiles"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.title  = ggplot2::element_text(face = "bold")
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(size = 3, max.overlaps = 15)
  } else {
    p <- p + ggplot2::geom_text(size = 3, hjust = -0.2, check_overlap = TRUE)
  }

  if (!is.null(output_pdf)) {
    ggplot2::ggsave(output_pdf, p, width = 8, height = 6)
    message("[INFO] PCA plot saved: ", output_pdf)
  } else {
    print(p)
  }
  invisible(list(pca = pca, var_exp = var_exp))
}


#' Flag Background-Exon Calibration Issues for a CNV Call Window
#'
#' Tests whether the test sample's ratio at *non-called* ("background")
#' exons in a plotted window falls outside the modelled 95% predictive
#' interval more often than the ~5% a well-calibrated interval implies.
#' A high background "outside" rate for a specific call can indicate real
#' signal extending beyond the called boundary, an atypical reference
#' match, or a technical/batch difference for that sample.
#'
#' This is a diagnostic flag only — it doesn't change the interval, the
#' call, or the confidence score. It's meant to prompt a manual look at
#' specific calls.
#'
#' @param ratio Numeric vector of log2(observed/expected) for every exon in
#'   the plotted window (background and affected together).
#' @param lo,hi Numeric vectors (same length as \code{ratio}) giving the
#'   95% predictive interval bounds at each exon.
#' @param is_affected Logical vector (same length); \code{TRUE} for exons
#'   already called as part of this CNV — excluded from the check, since
#'   those are expected to sit outside the interval.
#' @param min_n Minimum number of background exons required before
#'   flagging (default 5) — below this the percentage is too noisy on its
#'   own to test meaningfully.
#' @return A list with \code{n_background}, \code{n_outside},
#'   \code{pct_outside}, and \code{flag} — \code{flag} is \code{TRUE} when a
#'   one-sided binomial test of \code{n_outside} against a 5% null rate is
#'   significant at p < 0.05.
#' @export
check_background_calibration <- function(ratio, lo, hi, is_affected, min_n = 5) {
  bg <- !is_affected
  n_bg <- sum(bg, na.rm = TRUE)
  if (n_bg == 0) {
    return(list(n_background = 0L, n_outside = 0L, pct_outside = NA_real_, flag = FALSE))
  }
  outside <- (ratio[bg] < lo[bg]) | (ratio[bg] > hi[bg])
  outside[is.na(outside)] <- FALSE
  n_outside <- sum(outside)
  pct_outside <- 100 * n_outside / n_bg
  flag <- n_bg >= min_n &&
    stats::pbinom(n_outside - 1L, size = n_bg, prob = 0.05, lower.tail = FALSE) < 0.05
  list(n_background = n_bg, n_outside = n_outside, pct_outside = pct_outside, flag = flag)
}
