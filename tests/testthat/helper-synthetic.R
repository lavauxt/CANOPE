# Small synthetic exome: 3 chromosomes (chr2, chr10, chrX), 14 samples, one
# heterozygous deletion in S03 (chr2, gene G2_2) and a duplication in S07
# (chr10, gene G10_1).  Counts are negative-binomial with GC / capture bias.
# `order = "lex"` writes the table sorted as "chr10, chr2, chrX" (what
# `sort()` gives), `"karyo"` as "chr2, chr10, chrX" -- the *same* data.
make_synthetic_cohort <- function(dir, order = c("karyo", "lex"), ns = 14L, seed = 1L,
                                  spike = TRUE) {
  order <- match.arg(order)
  set.seed(seed)
  dir.create(file.path(dir, "bams"), recursive = TRUE, showWarnings = FALSE)
  chrs <- c("chr2", "chr10", "chrX")
  rows <- list(); pos0 <- 100000
  for (ch in chrs) {
    ng <- if (ch == "chrX") 3L else 5L
    pos <- pos0
    for (g in seq_len(ng)) {
      st <- pos + cumsum(sample(2000:6000, 8L)); en <- st + sample(120:260, 8L, TRUE)
      rows[[length(rows) + 1L]] <- data.frame(chromosome = ch, start = st, end = en,
                                              GENE = paste0("G", sub("chr", "", ch), "_", g))
      pos <- max(en) + 150000
    }
  }
  bed <- do.call(rbind, rows)
  n <- nrow(bed)
  snames <- sprintf("S%02d", seq_len(ns))
  sex <- rep(c("F", "M"), length.out = ns)
  gc <- pmin(pmax(rnorm(n, 0.5, 0.08), 0.2), 0.8)
  eff <- exp(rnorm(n, 0, 0.3)) * exp(-((gc - 0.5)^2) * 6)
  depth <- exp(rnorm(ns, 0, 0.1)) * 600
  M <- matrix(0, n, ns, dimnames = list(NULL, snames))
  for (j in seq_len(ns)) {
    mu <- eff * depth[j]
    isX <- bed$chromosome == "chrX"
    if (sex[j] == "M") mu[isX] <- mu[isX] * 0.5
    M[, j] <- rnbinom(n, mu = mu, size = 60)
  }
  if (spike) {
    sp <- function(s, gene, f) {
      i <- which(bed$GENE == gene)
      M[i, s] <<- rnbinom(length(i), mu = eff[i] * depth[match(s, snames)] * f, size = 60)
    }
    sp("S03", "G2_2", 0.5); sp("S07", "G10_1", 2)   # 1 and 4 copies: strong enough for a 14-sample cohort
  }
  perm <- if (order == "lex") order(bed$chromosome, bed$start) else seq_len(n)
  bed <- bed[perm, ]; M <- M[perm, ]; gc <- gc[perm]; rownames(bed) <- NULL
  utils::write.table(cbind(bed, as.data.frame(M)), file.path(dir, "counts.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(cbind(bed[, 1:4], GC_CONTENT = gc), file.path(dir, "gc.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(bed$chromosome, bed$start - 1L, bed$end, bed$GENE),
                     file.path(dir, "targets.bed"), sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  utils::write.table(data.frame(sample_name = snames, gender = sex), file.path(dir, "sample_table.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  for (s in snames) file.create(file.path(dir, "bams", paste0(s, ".bam")))
  invisible(dir)
}

run_synthetic <- function(dir, out, modechrom = "A", report = FALSE, ...) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  args <- list(
    reads_file = file.path(dir, "counts.tsv"), gc_file = file.path(dir, "gc.tsv"),
    samples = file.path(dir, "bams"), bed_file = file.path(dir, "targets.bed"),
    modechrom = modechrom, output_file = file.path(out, "CANOPE_T_CNVCall.csv"),
    rdata_output = file.path(out, "CANOPE_T_workspace.RData"), output_prefix = "CANOPE_T",
    p_value = 1e-8, Tnum = 6, generate_plots_flag = FALSE, pca_plot = FALSE,
    report = report, sample_table = file.path(dir, "sample_table.tsv"))
  args[names(list(...))] <- list(...)
  suppressMessages(suppressWarnings(do.call(run_canope, args)))
}
