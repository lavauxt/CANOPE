# Tests for canope_vcf_export.R

make_cnv_calls <- function() {
  data.frame(
    SAMPLE      = c("S1", "S2"),
    CNV         = c("DEL", "DUP"),
    INTERVAL    = c("chr1:1000-2000", "chr2:5000-6000"),
    KB          = c(1.0, 1.0),
    CHR         = c("chr1", "chr2"),
    MID_BP      = c(1500, 5500),
    TARGETS     = c("1..5", "10..12"),
    NUM_TARG    = c(5, 3),
    GENE        = c("BRCA1", "TP53"),
    MLCN        = c(1, 3),
    Q_SOME      = c(45, 30),
    NUM_REFS    = c(10, 8),
    REF_SAMPLES = c("S2;S3", "S1;S3"),
    Confidence  = c("HIGH", "MEDIUM"),
    stringsAsFactors = FALSE
  )
}

test_that("export_canope_to_vcf writes a well-formed multi-sample VCF", {
  calls <- make_cnv_calls()
  out_vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(out_vcf))

  expect_message(export_canope_to_vcf(calls, out_vcf), "Wrote 2 CNV records")
  lines <- readLines(out_vcf)

  expect_true(any(lines == "##fileformat=VCFv4.2"))
  expect_true(any(lines == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"))

  rec1 <- "chr1\t1000\t1_1000_2000_DEL\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-1001;END=2000;CONFIDENCE=HIGH;NUMTARG=5;REFS=S2,S3;NCOMP=10;QSOME=45\tGT:CN:QS\t0/1:1:45\t0/0:2:."
  rec2 <- "chr2\t5000\t2_5000_6000_DUP\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=1001;END=6000;CONFIDENCE=MEDIUM;NUMTARG=3;REFS=S1,S3;NCOMP=8;QSOME=30\tGT:CN:QS\t0/0:2:.\t0/1:3:30"

  expect_true(rec1 %in% lines)
  expect_true(rec2 %in% lines)

  # records must be ordered chr1 before chr2 (natural chromosome order)
  expect_lt(which(lines == rec1), which(lines == rec2))
})

test_that("REF_SAMPLES is comma-separated in the VCF INFO field, not semicolon-separated (regression)", {
  # Regression test for a VCF-spec violation: REFS is a Number=. INFO field,
  # which must use ',' between multiple values. The pipeline stores
  # REF_SAMPLES ';'-joined (fine in the CSV/report), but writing that
  # separator straight into the INFO field collided with ';', the INFO
  # key=value pair delimiter, corrupting the record for any VCF parser.
  calls <- make_cnv_calls()
  out_vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(out_vcf))
  suppressMessages(export_canope_to_vcf(calls, out_vcf))
  lines <- readLines(out_vcf)

  info_fields <- vapply(lines[!startsWith(lines, "#")], function(l) strsplit(l, "\t")[[1]][8], character(1))
  refs_tokens <- regmatches(info_fields, regexpr("REFS=[^;]*", info_fields))

  expect_equal(refs_tokens, c("REFS=S2,S3", "REFS=S1,S3"))
})

test_that("export_canope_to_vcf filters to a single sample when sample_name is given", {
  calls <- make_cnv_calls()
  out_vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(out_vcf))

  suppressMessages(export_canope_to_vcf(calls, out_vcf, sample_name = "S1"))
  lines <- readLines(out_vcf)

  expect_true(any(lines == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"))
  expect_equal(sum(!startsWith(lines, "#")), 1L)
})

test_that("export_canope_to_vcf warns and writes nothing for a sample with no calls", {
  calls <- make_cnv_calls()
  out_vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(out_vcf))

  expect_warning(res <- export_canope_to_vcf(calls, out_vcf, sample_name = "NOBODY"), "No CNVs found")
  expect_null(res)
  expect_false(file.exists(out_vcf))
})

test_that("export_canope_to_vcf handles empty input without writing a file", {
  out_vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(out_vcf))
  expect_message(res <- export_canope_to_vcf(data.frame(), out_vcf), "No CNVs to export")
  expect_null(res)
  expect_false(file.exists(out_vcf))
})
