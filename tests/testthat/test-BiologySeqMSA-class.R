test_that("BiologySeqMSA", {
  bss <- BiologySeqSet(
    c("ATGCAGGTAAACCCTACTGAG", "ATGCAGGTTACTGAG", "ATGCAGGTAACTGTG")
  )
  alnbs <- BiologySeqMSA(bss)

  expect_identical(
    AA(alnbs),
    Biostrings::AAStringSet(
      c("1" = "MQVNPTE", "2" = "MQV--TE", "3" = "MQV--TV")
    )
  )

  expect_identical(
    consDNA(alnbs),
    Biostrings::DNAString("ATGCAGGTA------ACTGAG")
  )

  expect_identical(
    consAA(alnbs),
    Biostrings::AAString("MQV--TE")
  )

  expect_identical(
    consDNA(alnbs),
    Biostrings::DNAString("ATGCAGGTA------ACTGAG")
  )


  expect_identical(
    consDNAfreq(alnbs),
    c(
      1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
      1.000, 1.000, 0.667, 0.667, 0.667, 0.667,
      0.667, 0.667, 0.667, 1.000, 1.000, 1.000,
      1.000, 0.667, 1.000
    )
  )

  expect_identical(
    consAAfreq(alnbs),
    c(1.000, 1.000, 1.000, 0.667, 0.667, 1.000, 0.667)
  )
})
