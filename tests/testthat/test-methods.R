test_that("DNA<-, BiologySeq", {
  bs <- BiologySeq("ATGAAA---")

  DNA(bs) <- "ATGAAA"
  expect_identical(AA(bs), Biostrings::AAString("MK"))

  DNA(bs)[1:3] <- "A"
  expect_identical(AA(bs), Biostrings::AAString("KK"))
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


test_that("DNAsite_by_AA, BiologySeq", {
  bs <- BiologySeq("ATGAAA---")
  expect_identical(DNAsite_by_AA(bs, 2, 3)$seq_nt, "AAA---")
})

test_that("DNAsite_by_AA, BiologySeqSet", {
  bss <- BiologySeqSet(c("ATG---CCC", "ATG---TTT"))
  expect_identical(
    DNAsite_by_AA(bss, 2, 3)$seq_nt,
    c("---CCC", "---TTT")
  )
})


test_that("AAsite_by_DNA, BiologySeq", {
  bs <- BiologySeq("ATGAAA---")
  expect_identical(AAsite_by_DNA(bs, 2, 4)$seq_aa, "MK")
})

test_that("AAsite_by_DNA, BiologySeqSet", {
  bss <- BiologySeqSet(c("ATG---CCC", "ATG---TTT"))
  expect_identical(AAsite_by_DNA(bss, 1, 9)$seq_aa, c("M-P", "M-F"))
})
