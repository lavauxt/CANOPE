# Tests for canope_bed_process.R.
#
# process_bed_file()'s STANDARD (panel-based) and REGEN (TxDb-based) modes
# need real annotation resources (a reference panel BED / TxDb.* +
# org.Hs.eg.db) to exercise properly and are integration-level; not unit
# tested here. "NO" mode is fully self-contained (just parses the BED name
# column and re-sorts), so it's covered as a deterministic smoke test below.
#
# The REGEN mode's per-gene NM_-transcript-preference fix (see the review
# notes) is exercised in isolation further down, since the surrounding
# TxDb/org.Hs.eg.db annotation machinery can't run without those packages.

test_that("process_bed_file NO mode extracts gene/exon from the name column and sorts naturally", {
  input_bed <- tempfile(fileext = ".bed")
  output_bed <- tempfile(fileext = ".bed")
  on.exit(unlink(c(input_bed, output_bed)))

  writeLines(c(
    "chr1\t1000\t1100\tBRCA1_ex2",
    "chr1\t100\t200\tBRCA1_ex1",
    "chr2\t50\t60\tTP53_ex1"
  ), input_bed)

  process_bed_file(
    input_bed, output_bed,
    bed_process = "NO",
    exon_sep = "_",
    auto_exon_number = FALSE
  )

  out <- read.table(output_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  expect_equal(out$V1, c("chr1", "chr1", "chr2"))
  expect_equal(out$V2, c(100L, 1000L, 50L))
  expect_equal(out$V3, c(200L, 1100L, 60L))
  expect_equal(out$V4, c("BRCA1", "BRCA1", "TP53"))
  expect_equal(out$V5, c(1L, 2L, 1L))
})

test_that("the REGEN per-gene NM_-preference logic keeps genes whose only transcripts are non-coding", {
  # Isolated reproduction of the fix: previously
  # `if (any(grepl("^NM_", ref_df$Transcript))) ref_df <- ref_df[grepl(...), ]`
  # applied globally, so a single NM_ transcript anywhere would wipe out
  # every gene that only has NR_/XM_ transcripts. The fix scopes the
  # preference to each gene individually.
  ref_df <- data.frame(
    Gene       = c("GENE_CODING", "GENE_CODING", "GENE_NONCODING", "GENE_NONCODING"),
    Transcript = c("NM_000001", "XM_000002", "NR_000003", "NR_000004"),
    stringsAsFactors = FALSE
  )

  is_nm <- grepl("^NM_", ref_df$Transcript)
  gene_has_nm <- stats::ave(is_nm, ref_df$Gene, FUN = any)
  fixed <- ref_df[is_nm | !gene_has_nm, ]

  # GENE_CODING: only its NM_ transcript should survive.
  expect_equal(fixed$Transcript[fixed$Gene == "GENE_CODING"], "NM_000001")
  # GENE_NONCODING has no NM_ transcript at all -- both of its (non-coding)
  # transcripts must be kept, not silently dropped.
  expect_setequal(fixed$Transcript[fixed$Gene == "GENE_NONCODING"],
                   c("NR_000003", "NR_000004"))

  # Demonstrate what the OLD (buggy) global filter would have done, for
  # contrast: it would have wiped out GENE_NONCODING entirely.
  old_buggy <- ref_df[grepl("^NM_", ref_df$Transcript), ]
  expect_false("GENE_NONCODING" %in% old_buggy$Gene)
})
