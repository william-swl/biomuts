test_that("BiologySeqMSA-1", {
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
    Biostrings::DNAString("ATGCAGGTAACT------GAG")
  )

  expect_identical(
    consAA(alnbs),
    Biostrings::AAString("MQVT--E")
  )

  expect_identical(
    consDNAfreq(alnbs),
    c(
      1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
      1.000, 0.667, 1.000, 0.667, 0.667, 0.667, 0.667,
      0.667, 0.667, 0.667, 0.667, 1.000, 0.667, 1.000
    )
  )

  expect_identical(
    consAAfreq(alnbs),
    c(1.000, 1.000, 1.000, 0.667, 0.667, 0.667, 0.667)
  )

  expect_true(is(alnbs[[3]], 'BiologySeq'))

  expect_identical(length(alnbs), as.integer(3))

})

test_that("BiologySeqMSA-2", {
  bss <- BiologySeqSet(
    c(
      "ATGCAGGTTACTGTG", "ATGCAGGTAACTGTG",
      "ATGCAGGTAAACCCTACTGAG", "ATGCAGGTAACTGAC"
    )
  )
  alnbs <- BiologySeqMSA(bss, order = "aligned")

  expect_identical(
    AA(alnbs),
    Biostrings::AAStringSet(
      c(
        "3" = "MQVNPTE", "4" = "MQVTD--",
        "1" = "MQVTV--", "2" = "MQVTV--"
      )
    )
  )
})
