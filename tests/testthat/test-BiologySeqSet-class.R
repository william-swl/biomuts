test_that("BiologySeqSet", {
  expect_error(BiologySeqSet(c("ATG--CCC", "ATG--TTT")))
  bss <- BiologySeqSet(c("ATG---CCC", "ATG---TTT"))
  expect_identical(
    AA(bss),
    Biostrings::AAStringSet(c("1" = "M-P", "2" = "M-F"))
  )

  expect_true(is(bss[[1]], 'BiologySeq'))

})
