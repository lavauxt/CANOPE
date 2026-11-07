test_that("canope_chrom_class classifies labels with or without chr prefix", {
  expect_equal(canope_chrom_class(c("chr1", "22", "chrX", "Y", "chrM", "chrUn_gl000", "23", "24")),
               c("autosome", "autosome", "X", "Y", "MT", "other", "X", "Y"))
  expect_true(is.na(canope_chrom_class(NA_character_)))
})

test_that("autosome-only run warns that gonosomes were not analysed", {
  inp <- c("chr1", "chr2", "chrX", "chrY")
  st <- canope_chromosome_status(inp, c("chr1", "chr2"), modechrom = "A")
  expect_false(st$autosomes_and_gonosomes)
  expect_equal(st$status$analysed, c(TRUE, FALSE, FALSE))
  expect_true(any(grepl("not both analysed", st$messages)))
  expect_true(any(grepl("chrX targets are present in the BED but were NOT analysed", st$messages)))
})

test_that("XX / XY modes warn that the autosomes were not analysed", {
  inp <- c("chr1", "chrX", "chrY")
  xx <- canope_chromosome_status(inp, "chrX", modechrom = "XX")
  expect_false(xx$autosomes_and_gonosomes)
  expect_true(any(grepl("Autosomal targets are present", xx$messages)))
  expect_true(any(grepl("not both analysed", xx$messages)))
  expect_match(xx$status$reason[xx$status$class == "Y"], "chrX only")

  xy <- canope_chromosome_status(inp, "chrX", modechrom = "XY", removeY = TRUE)
  expect_match(xy$status$reason[xy$status$class == "Y"], "removeY")
})

test_that("a chromosome missing from the BED is reported as absent, not as 'not analysed'", {
  st <- canope_chromosome_status(c("chr1", "chr2"), c("chr1", "chr2"), modechrom = "A")
  expect_match(st$status$reason[st$status$class == "X"], "not present")
  expect_false(any(grepl("NOT analysed", st$messages)))
})

test_that("canope_infer_sex flags the expected sex from chrX / chrY depth", {
  m <- rbind(chrX = c(F1 = 1.0, M1 = 0.5, odd = 0.9), chrY = c(0.01, 0.4, 0.5))
  expect_equal(canope_infer_sex(m)$inferred_sex, c("F", "M", "unclear"))
  expect_null(canope_infer_sex(m[2, , drop = FALSE]))
})
