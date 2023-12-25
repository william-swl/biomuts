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

test_that("call_AAmut", {
  query <- Biostrings::AAString("MQVNPTE")
  ref <- Biostrings::AAString("MQVCTTE")

  expect_identical(call_AAmut(query, ref), c("C4N", "T5P"))
})



test_that("call_AAmutSet", {
  bss <- BiologySeqSet(
    c("ATGCAGGTAAACCCTACTGAG", "ATGCAGGTTACTGAG", "ATGCAGGTAACTGTG")
  )
  alnbs <- BiologySeqMSA(bss)

  expect_snapshot(call_AAmutSet(alnbs))
  expect_snapshot(call_AAmutSet(alnbs, ref = "3"))
})


test_that("count_muts", {
  mut_list <- list(
    mut1 = c("D123G", "D6F", "C878C", "C878T"),
    mut2 = c("D123G", "D6F"),
    mut3 = c("C878C", "D123G", "C878T")
  )

  muts <- BiologyAAmutSet(mut_list)
  numbering(muts) <- c("6" = "6", "123" = "123")

  expect_snapshot(count_muts(muts))

  expect_snapshot(count_muts(muts, use_numbering = TRUE))
})


test_that("aa_info", {
  expect_snapshot(aa_info("C"))
})


test_that("compare_aa", {
  expect_snapshot(compare_aa(c("A", "T"), "C"))
})


test_that("number_mut", {
  expect_identical(
    number_mut(c("A12T", "C18P"), c("12" = "12", "18" = "18")),
    c("A[12]T", "C[18]P")
  )
})
