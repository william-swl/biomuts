test_that("BiologySeqSet", {
  expect_error(BiologySeqSet(c("ATG--CCC", "ATG--TTT")))
  bss <- BiologySeqSet(c("ATG---CCC", "ATG---TTT"))
  expect_identical(
    AA(bss),
    Biostrings::AAStringSet(c("1" = "M-P", "2" = "M-F"))
  )

  expect_true(is(bss[[1]], "BiologySeq"))
})

test_that("DNA<-, BiologySeqSet", {
  bss <- BiologySeqSet(c("ATG---CCC", "ATG---TTT"))

  DNA(bss) <- c("ATGCCCCCC", "ATGCCCTTT")
  expect_identical(
    AA(bss),
    Biostrings::AAStringSet(c("1" = "MPP", "2" = "MPF"))
  )

  DNA(bss) <- c("ATG---CCC", "ATG---TTT", "ATG---TTT")
  expect_identical(
    AA(bss),
    Biostrings::AAStringSet(c("1" = "M-P", "2" = "M-F", "3" = "M-F"))
  )
})


test_that("DNAsite_by_AA, BiologySeqSet", {
  bss <- BiologySeqSet(c("ATG---CCC", "ATG---TTT"))
  expect_identical(
    DNAsite_by_AA(bss, 2, 3)$seq_nt,
    c("---CCC", "---TTT")
  )
})


test_that("AAsite_by_DNA, BiologySeqSet", {
  bss <- BiologySeqSet(c("ATG---CCC", "ATG---TTT"))
  expect_identical(AAsite_by_DNA(bss, 1, 9)$seq_aa, c("M-P", "M-F"))
})
