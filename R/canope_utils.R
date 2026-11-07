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

  dt[, ._orig_row := .I]
  dt[, ._key := paste(`._chrom`, `._start`, `._end`, `._gene`, sep = "\r")]

  dup_rows <- duplicated(dt, by = "._key")
  if (any(dup_rows)) {
    warning(sprintf(
      "Found %d duplicate BED row(s) (same chrom, start, end, gene); duplicates will share the exon_number of their first occurrence.",
      sum(dup_rows)), immediate. = TRUE)
  }

  dt_unique <- dt[!dup_rows]

  chrom_levels <- c(paste0("chr", c(1:22, "X", "Y", "M")), c(as.character(1:22), "X", "Y", "M"))
  dt_unique[, .chrom_fac := factor(`._chrom`, levels = unique(c(chrom_levels, unique(`._chrom`))))]
  data.table::setorder(dt_unique, .chrom_fac, `._start`, `._end`)
  dt_unique[, .chrom_fac := NULL]

  dt_unique[, exon_number := seq_len(.N), by = "._gene"]
  # Rows with no real gene (NA / "" / "." -- e.g. off-target filler intervals
  # handled with off_target_handling = "na") must stay un-numbered, not be
  # numbered 1..k as if "NA" were a gene.
  no_gene_u <- is.na(dt_unique$`._gene`) | as.character(dt_unique$`._gene`) %in% c("", ".", "NA")
  dt_unique[no_gene_u, exon_number := NA_integer_]

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
#' Capture-based coverage often drops off right at the
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

  numbered <- assign_exon_numbers_per_gene(bed_df)

  dt <- data.table::as.data.table(numbered)
  dt[, .orig_row := .I]
  data.table::setnames(dt, c(chrom_col, start_col, end_col, gene_col),
                       c("chrom", "start", "end", "gene"))

  dt[, is_first := exon_number == 1L]
  dt[, is_last  := exon_number == max(exon_number), by = "gene"]
  no_gene <- is.na(dt$gene) | dt$gene %in% c("", ".", "Unknown")
  dt[no_gene, c("is_first", "is_last") := FALSE]

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

  right_extend <- integer(n)  
  left_extend  <- integer(n)  
  if (n > 1) {
    same_chr_pair <- chrom_id[-n] == chrom_id[-1]
    gap        <- pmax(start_v[-1] - end_v[-n] - 1L, 0L)
    want_left  <- ifelse(is_last_v[-n],  padding, 0L)  
    want_right <- ifelse(is_first_v[-1], padding, 0L) 
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

  data.table::setorder(dt, .orig_row)  
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

  raw$start <- raw$start + 1L
  padded <- pad_gene_terminal_exons(raw, padding = padding, chr_lengths = chr_lengths, verbose = verbose)
  padded$start <- padded$start - 1L
  keep   <- if (ncol(raw) >= 4) colnames(raw) else setdiff(colnames(raw), "GENE")
  out_df <- padded[, keep, drop = FALSE]

  dir.create(dirname(output_bed), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(out_df, file = output_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  if (verbose) message("[INFO] pad_bed_file: padded BED written: ", output_bed)
  invisible(output_bed)
}

#' Compute gap-inserted x-axis positions for CNV window plots
#'
#' It also returns \code{gene_group}, a per-position integer that
#' increments at every such boundary. Passing this as the \code{group}
#' aesthetic on a \code{geom_ribbon()} keeps that visual gap genuinely
#' blank -- otherwise ggplot draws a single connected ribbon straight
#' across it. (Points are never joined by a line in these plots, so this
#' only matters for the ribbon.)
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
  file.exists(paste0(bam, ".bai")) || file.exists(sub("\\.bam$", ".bai", bam))
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
  use_groups <- !is.null(color_by) && length(color_by) == nrow(pca_df)
  # The grouping column has to exist *before* ggplot() captures the data frame;
  # adding it afterwards leaves the plot looking for a column that isn't there.
  if (use_groups) pca_df$Group <- as.factor(ifelse(is.na(color_by), "NA", as.character(color_by)))

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2))

  if (use_groups) {
    p <- p + ggplot2::aes(colour = Group)
    # stat_ellipse needs >= 3 points per group to be defined
    if (any(table(pca_df$Group) >= 3L))
      p <- p + suppressWarnings(ggplot2::stat_ellipse(type = "norm", linetype = "dashed"))
    p <- p + ggplot2::geom_point(size = 3)
  } else {
    p <- p + ggplot2::geom_point(size = 3, colour = "steelblue")
  }

  p <- p +
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
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = Sample), size = 3, max.overlaps = 15)
  } else {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = Sample), size = 3, hjust = -0.2, check_overlap = TRUE)
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


#' Classify chromosome labels as autosome / X / Y / MT / other
#'
#' @param chr Character vector of chromosome labels (with or without a
#'   \code{chr} prefix; \code{23}/\code{24} are read as X/Y).
#' @return Character vector of the same length: \code{"autosome"}, \code{"X"},
#'   \code{"Y"}, \code{"MT"} or \code{"other"}.
#' @export
canope_chrom_class <- function(chr) {
  s <- toupper(sub("^chr", "", as.character(chr), ignore.case = TRUE))
  out <- rep("other", length(s))
  out[s %in% as.character(1:22)] <- "autosome"
  out[s %in% c("X", "23")] <- "X"
  out[s %in% c("Y", "24")] <- "Y"
  out[s %in% c("M", "MT")] <- "MT"
  out[is.na(s)] <- NA_character_
  out
}


#' Which chromosome classes were (not) analysed in a CANOPE run
#'
#' One run analyses a single \code{modechrom}: \code{"A"} (autosomes only),
#' \code{"XX"} (chrX only) or \code{"XY"} (chrX + chrY). A run that covers only
#' autosomes silently leaves the gonosomes un-analysed (and vice versa), which
#' is easy to miss when reading the calls — this makes it explicit.
#'
#' @param input_chromosomes Chromosome labels present in the BED / count table
#'   \emph{before} any \code{modechrom} filtering.
#' @param analysed_chromosomes Chromosome labels of the targets that actually
#'   reached the CNV caller.
#' @param modechrom,removeY The run settings.
#' @return A list with \code{status} (one row per class: autosome, X, Y with
#'   \code{in_input}, \code{analysed}, \code{reason}), \code{messages}
#'   (character vector of human-readable warnings) and
#'   \code{autosomes_and_gonosomes} (\code{TRUE} only if autosomes and
#'   chrX were both analysed).
#' @export
canope_chromosome_status <- function(input_chromosomes, analysed_chromosomes,
                                     modechrom = "A", removeY = TRUE) {
  cls_in  <- canope_chrom_class(input_chromosomes)
  cls_ana <- canope_chrom_class(analysed_chromosomes)
  classes <- c("autosome", "X", "Y")
  label   <- c(autosome = "Autosomes", X = "chrX", Y = "chrY")
  in_in   <- vapply(classes, function(k) any(cls_in  == k, na.rm = TRUE), logical(1))
  in_ana  <- vapply(classes, function(k) any(cls_ana == k, na.rm = TRUE), logical(1))

  reason <- setNames(rep("", 3L), classes)
  for (k in classes) {
    if (!in_in[[k]]) {
      reason[[k]] <- "not present in the BED / count table"
    } else if (!in_ana[[k]]) {
      reason[[k]] <- if (k == "autosome" && modechrom != "A") {
        sprintf("excluded by modechrom = '%s' (sex-chromosome mode)", modechrom)
      } else if (k == "X" && modechrom == "A") {
        "excluded by modechrom = 'A' (autosomes only)"
      } else if (k == "Y" && modechrom == "A") {
        "excluded by modechrom = 'A' (autosomes only)"
      } else if (k == "Y" && modechrom == "XX") {
        "excluded by modechrom = 'XX' (chrX only)"
      } else if (k == "Y" && isTRUE(removeY)) {
        "excluded by removeY = TRUE"
      } else {
        "all targets were removed by QC / coverage filters"
      }
    }
  }

  msgs <- character(0)
  for (k in classes) {
    if (!in_in[[k]]) {
      msgs <- c(msgs, sprintf("%s: no targets in the BED / count table, so %s cannot be called.",
                              label[[k]], if (k == "autosome") "autosomal CNVs" else paste0(label[[k]], " CNVs")))
    } else if (!in_ana[[k]]) {
      msgs <- c(msgs, sprintf("%s targets are present in the BED but were NOT analysed: %s.",
                              if (k == "autosome") "Autosomal" else label[[k]], reason[[k]]))
    }
  }
  if (in_in[["autosome"]] && !in_ana[["autosome"]] || (in_in[["X"]] && !in_ana[["X"]])) {
    msgs <- c(msgs, sprintf(
      "Autosomes and gonosomes were not both analysed in this run (modechrom = '%s'). Run again with %s to cover the rest; sex chromosomes need sex-matched references.",
      modechrom, if (modechrom == "A") "modechrom = 'XX' (females) and/or 'XY' (males)" else "modechrom = 'A'"))
  }

  list(
    status = data.frame(class = classes, label = unname(label[classes]),
                        in_input = unname(in_in), analysed = unname(in_ana),
                        reason = unname(reason), stringsAsFactors = FALSE),
    messages = msgs,
    autosomes_and_gonosomes = isTRUE(in_ana[["autosome"]] && in_ana[["X"]])
  )
}


#' Infer each sample's sex from sex-chromosome coverage
#'
#' Heuristic on the depth of chrX (and chrY when present) relative to the
#' autosomal median: a female has ~1x chrX and ~0x chrY, a male ~0.5x chrX and
#' a clearly non-zero chrY. Only meant to flag an annotation / \code{modechrom}
#' mismatch, not to replace a proper sex check.
#'
#' @param chrom_depth Numeric matrix, rows = chromosomes (named), columns =
#'   samples: median target depth on the chromosome divided by the sample's
#'   autosomal median.
#' @param x_threshold chrX ratio below which a sample looks male. Default 0.75.
#' @param y_threshold chrY ratio at/above which a sample looks male. Default 0.1.
#' @return Data frame \code{sample, x_ratio, y_ratio, inferred_sex}
#'   (\code{"F"}, \code{"M"} or \code{"unclear"}), or \code{NULL} when the
#'   table has no chrX.
#' @export
canope_infer_sex <- function(chrom_depth, x_threshold = 0.75, y_threshold = 0.1) {
  if (is.null(chrom_depth) || is.null(dim(chrom_depth))) return(NULL)
  cls <- canope_chrom_class(rownames(chrom_depth))
  if (!any(cls == "X", na.rm = TRUE)) return(NULL)
  xr <- as.numeric(chrom_depth[which(cls == "X")[1L], ])
  yr <- if (any(cls == "Y", na.rm = TRUE)) as.numeric(chrom_depth[which(cls == "Y")[1L], ]) else rep(NA_real_, length(xr))
  male   <- !is.na(xr) & xr < x_threshold  & (is.na(yr) | yr >= y_threshold)
  female <- !is.na(xr) & xr >= x_threshold & (is.na(yr) | yr <  y_threshold)
  data.frame(sample = colnames(chrom_depth), x_ratio = xr, y_ratio = yr,
             inferred_sex = ifelse(male, "M", ifelse(female, "F", "unclear")),
             stringsAsFactors = FALSE)
}


#' Normalise a sex / gender annotation to "F" / "M"
#' @param x Character vector (F, M, female, male, XX, XY ...).
#' @return \code{"F"}, \code{"M"} or \code{NA}.
#' @noRd
.canope_norm_sex <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(x %in% c("F", "FEMALE", "XX", "W"), "F",
         ifelse(x %in% c("M", "MALE", "XY"), "M", NA_character_))
}


#' Guess the genome build label from a BSgenome package / FASTA name
#' @param ... Character strings to inspect (BSgenome package, FASTA path,
#'   \code{genome_version}).
#' @return \code{"hg38"}, \code{"hg19"}, \code{"T2T-CHM13"} or \code{"unknown"}.
#' @export
canope_infer_genome_build <- function(...) {
  x <- tolower(paste(unlist(list(...)), collapse = " "))
  if (grepl("hg38|grch38|b38", x)) "hg38"
  else if (grepl("hg19|grch37|b37", x)) "hg19"
  else if (grepl("t2t|chm13", x)) "T2T-CHM13"
  else "unknown"
}

#' Export CANOPE CNV calls to VCF
#'
#' Converts a CANOPE CNV-call data frame (as returned by \code{\link{run_canope}}
#' and optionally annotated by \code{\link{score_canope_confidence}}) into a
#' standard VCF file.  When \code{sample_name} is \code{NULL}, a single
#' multi-sample VCF is written; otherwise a single-sample VCF is produced.
#'
#' @param cnv_calls   Data frame of CANOPE CNV calls (the \code{xcnv}-style
#'   table produced by \code{print_cnvs()} / \code{run_canope()}).
#' @param output_vcf  Output VCF file path.
#' @param sample_name Optional character string. If supplied, write a
#'   single-sample VCF for that sample only.
#' @param source      Source description for the VCF header. Default \code{"CANOPE"}.
#' @param reference   Reference genome assembly label (e.g. \code{"hg19"}).
#'
#' @return Invisibly returns the path written.
#' @export
export_canope_to_vcf <- function(cnv_calls, output_vcf, sample_name = NULL,
                                 source = "CANOPE", reference = "hg19") {
  if (is.null(cnv_calls) || nrow(cnv_calls) == 0) {
    message("[INFO] No CNVs to export.")
    return(invisible(NULL))
  }

  if (!is.null(sample_name)) {
    cnv_calls <- cnv_calls[cnv_calls$SAMPLE == sample_name, , drop = FALSE]
    if (nrow(cnv_calls) == 0) {
      warning("No CNVs found for sample ", sample_name)
      return(invisible(NULL))
    }
  }

  if (!all(c("Start", "End") %in% colnames(cnv_calls))) {
    parsed <- lapply(cnv_calls$INTERVAL, parse_canope_interval)
    cnv_calls$Start <- vapply(parsed, function(p) p$start, integer(1))
    cnv_calls$End   <- vapply(parsed, function(p) p$end,   integer(1))
  }

  safe_char <- function(x, default = ".") {
    if (length(x) == 0 || is.na(x) || x == "NA") default else as.character(x)
  }
  samples <- sort(unique(as.character(cnv_calls$SAMPLE)))

  header <- c(
    "##fileformat=VCFv4.2",
    paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
    paste0("##source=", source),
    paste0("##reference=", reference),
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant (DUP/DEL)">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=CONFIDENCE,Number=1,Type=String,Description="Best confidence level (HIGH/MEDIUM/LOW) among the carrier samples; see FORMAT/CF per sample">',
    '##INFO=<ID=NUMTARG,Number=1,Type=Integer,Description="Number of targets (exons) spanned">',
    '##INFO=<ID=REFS,Number=.,Type=String,Description="Reference samples used (comma-separated); only written when a single sample carries the call, otherwise they differ per sample">',
    '##INFO=<ID=NCOMP,Number=1,Type=Integer,Description="Largest number of reference samples among the carrier samples">',
    '##INFO=<ID=QSOME,Number=1,Type=Integer,Description="Highest phred-scaled probability some CNV exists in this interval (SQ) among the carrier samples; see FORMAT/QS per sample">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Most-likely copy number (MLCN)">',
    '##FORMAT=<ID=QS,Number=1,Type=Integer,Description="Per-sample Q_SOME phred score">',
    '##FORMAT=<ID=CF,Number=1,Type=String,Description="Per-sample confidence level (HIGH/MEDIUM/LOW)">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=DEL,Description="Deletion">',
    paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
           paste(samples, collapse = "\t"))
  )

  key_df <- unique(data.frame(
    Chromosome = as.character(cnv_calls$CHR),
    Start      = as.integer(cnv_calls$Start),
    End        = as.integer(cnv_calls$End),
    Type       = as.character(cnv_calls$CNV),
    stringsAsFactors = FALSE
  ))

  chrom_order <- function(chr) {
    stripped <- sub("^chr", "", chr)
    num      <- suppressWarnings(as.integer(stripped))
    ifelse(is.na(num), Inf, num)
  }
  key_df <- key_df[order(chrom_order(key_df$Chromosome),
                         key_df$Chromosome,
                         key_df$Start, key_df$End, key_df$Type), , drop = FALSE]

  make_sample_field <- function(row_df) {
    if (nrow(row_df) == 0) return("0/0:2:.:.")
    mlcn <- suppressWarnings(as.integer(row_df$MLCN[1]))
    qs <- suppressWarnings(as.integer(row_df$Q_SOME[1]))
    qs_str <- if (is.na(qs)) "." else as.character(qs)
    cf <- if ("Confidence" %in% colnames(row_df)) safe_char(row_df$Confidence[1]) else "."
    # An unknown copy number must not be written as a normal (CN=2, 0/0) call.
    if (is.na(mlcn)) return(paste0("./.:.:", qs_str, ":", cf))
    gt <- if (mlcn == 0L) "1/1" else if (mlcn < 2L) "0/1" else if (mlcn > 2L) "0/1" else "0/0"
    paste0(gt, ":", mlcn, ":", qs_str, ":", cf)
  }

  records <- character(nrow(key_df))
  for (i in seq_len(nrow(key_df))) {
    key_row  <- key_df[i, ]
    chrom    <- as.character(key_row[["Chromosome"]])
    start    <- as.integer(key_row[["Start"]])
    end      <- as.integer(key_row[["End"]])
    cnv_type <- toupper(as.character(key_row[["Type"]]))
    is_dup   <- cnv_type == "DUP"
    alt      <- if (is_dup) "<DUP>" else "<DEL>"
    sv_type  <- if (is_dup) "DUP" else "DEL"
    sv_len   <- if (is_dup) end - start + 1L else -(end - start + 1L)

    row_hits <- cnv_calls[
      as.character(cnv_calls$CHR) == chrom &
      as.integer(cnv_calls$Start) == start &
      as.integer(cnv_calls$End)   == end   &
      toupper(as.character(cnv_calls$CNV)) == cnv_type, , drop = FALSE]
    representative <- row_hits[1, , drop = FALSE]

    id <- paste0(gsub("^chr", "", chrom), "_", start, "_", end, "_", sv_type)

    info_items <- c(
      paste0("SVTYPE=", sv_type),
      paste0("SVLEN=", sv_len),
      paste0("END=", end)
    )

    # Cohort-level INFO fields summarise *all* carriers of this exact call
    # (previously they were read from the first carrier only).
    if ("Confidence" %in% colnames(row_hits)) {
      rank <- c(HIGH = 3L, MEDIUM = 2L, LOW = 1L)
      cf_all <- toupper(as.character(row_hits$Confidence))
      cf_all <- cf_all[cf_all %in% names(rank)]
      if (length(cf_all) > 0)
        info_items <- c(info_items, paste0("CONFIDENCE=", names(rank)[max(rank[cf_all])]))
    }

    if ("NUM_TARG" %in% colnames(representative)) {
      nt <- safe_char(representative$NUM_TARG)
      if (nt != ".") info_items <- c(info_items, paste0("NUMTARG=", nt))
    }
    if ("REF_SAMPLES" %in% colnames(representative) && nrow(row_hits) == 1L) {
      refs_val <- safe_char(representative$REF_SAMPLES, "")
      if (nzchar(refs_val)) {
        refs_vcf <- gsub(";", ",", refs_val, fixed = TRUE)
        info_items <- c(info_items, paste0("REFS=", refs_vcf))
      }
    }
    if ("NUM_REFS" %in% colnames(row_hits)) {
      nr <- suppressWarnings(max(as.integer(row_hits$NUM_REFS), na.rm = TRUE))
      if (is.finite(nr)) info_items <- c(info_items, paste0("NCOMP=", nr))
    }
    if ("Q_SOME" %in% colnames(row_hits)) {
      qm <- suppressWarnings(max(as.integer(row_hits$Q_SOME), na.rm = TRUE))
      if (is.finite(qm)) info_items <- c(info_items, paste0("QSOME=", qm))
    }

    format_field  <- "GT:CN:QS:CF"
    sample_fields <- vapply(samples, function(s) {
      sample_row <- row_hits[row_hits$SAMPLE == s, , drop = FALSE]
      make_sample_field(sample_row)
    }, character(1))

    records[i] <- paste(chrom, start, id, "N", alt, ".", "PASS",
                        paste(info_items, collapse = ";"),
                        format_field,
                        paste(sample_fields, collapse = "\t"),
                        sep = "\t")
  }

  dir.create(dirname(output_vcf), showWarnings = FALSE, recursive = TRUE)
  con <- file(output_vcf, "w")
  on.exit(close(con), add = TRUE)
  writeLines(header,  con)
  writeLines(records, con)
  message("[INFO] Wrote ", length(records), " CNV records to ", output_vcf)
  invisible(output_vcf)
}
