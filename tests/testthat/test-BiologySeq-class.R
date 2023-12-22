test_that("BiologySeq", {
  expect_error(BiologySeq("ATGAA"))
  expect_error(BiologySeq("ATGAAA--T"))

  bs <- BiologySeq("ATGAAA---")
  expect_identical(DNA(bs), Biostrings::DNAString("ATGAAA---"))
  expect_identical(AA(bs), Biostrings::AAString("MK-"))
})

test_that("DNA<-, BiologySeq", {
  bs <- BiologySeq("ATGAAA---")

  DNA(bs) <- "ATGAAA"
  expect_identical(AA(bs), Biostrings::AAString("MK"))

  DNA(bs)[1:3] <- "A"
  expect_identical(AA(bs), Biostrings::AAString("KK"))
})

test_that("DNAsite_by_AA, BiologySeq", {
  bs <- BiologySeq("ATGAAA---")
  expect_identical(DNAsite_by_AA(bs, 2, 3)$seq_nt, "AAA---")
})

test_that("AAsite_by_DNA, BiologySeq", {
  bs <- BiologySeq("ATGAAA---")
  expect_identical(AAsite_by_DNA(bs, 2, 4)$seq_aa, "MK")
})
