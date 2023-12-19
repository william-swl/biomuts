test_that("BiologySeq", {
  expect_error(BiologySeq("ATGAA"))
  expect_error(BiologySeq("ATGAAA--T"))

  bs <- BiologySeq("ATGAAA---")
  expect_identical(DNA(bs), Biostrings::DNAString("ATGAAA---"))
  expect_identical(AA(bs), Biostrings::AAString("MK-"))
})
