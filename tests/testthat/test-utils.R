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


test_that("dna2aa", {
  expect_snapshot(call_mut(BiologySeq("ATGCCCTTT"), BiologySeq("CTTCCCTTC")))
})


test_that("call_mut", {
  ref <- BiologySeq("CTTCCCTTC")
  query <- BiologySeq("ATGCCCTTT")

  expect_snapshot(call_mut(query, ref))
})




test_that("count_muts", {
  mut_list <- list(
    mut1 = c("D123G", "D6F", "C878C", "C878T"),
    mut2 = c("D123G", "D6F"),
    mut3 = c("C878C", "D123G", "C878T")
  )

  muts <- BiologyAAmutSet(mut_list)

  expect_snapshot(count_muts(muts))
})
