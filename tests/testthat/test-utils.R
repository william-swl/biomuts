test_that("DNA_gaps_corr", {
  expect_identical(DNA_gaps_corr("AT---G"), "ATG---")
  expect_identical(DNA_gaps_corr("A---TG"), "---ATG")
  expect_identical(DNA_gaps_corr("---ATG"), "---ATG")
  expect_error(DNA_gaps_corr("AT--GG"))
  expect_identical(DNA_gaps_corr("AT--GG", dismiss_no3gaps = TRUE), "AT--GG")
})


test_that("dna2aa", {
  expect_identical(dna2aa("ATG---AAA"), "M-K")
  expect_identical(dna2aa(Biostrings::DNAString("ATG---")), "M-")
  expect_error(dna2aa("ATG--"))
})
