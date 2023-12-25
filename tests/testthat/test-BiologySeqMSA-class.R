test_that("BiologySeqMSA-1", {
  bss <- BiologySeqSet(
    c("ATGCAGGTAAACCCTACTGAG", "ATGCAGGTTACTGAG", "ATGCAGGTAACTGTG")
  )
  alnbs <- BiologySeqMSA(bss)

  expect_error(BiologySeqMSA(bss, skip_align = TRUE))

  expect_identical(
    as.character(AA(alnbs)),
    c("1" = "MQVNPTE", "2" = "MQV--TE", "3" = "MQV--TV")
  )

  expect_identical(
    toString(consDNA(alnbs)),
    "ATGCAGGTA------ACTGAG"
  )

  expect_identical(toString(consAA(alnbs)), "MQV--TE")

  expect_identical(
    consDNAfreq(alnbs),
    c(
      1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
      1.000, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667,
      0.667, 1.000, 1.000, 1.000, 1.000, 0.667, 1.000
    )
  )

  expect_identical(
    consAAfreq(alnbs),
    c(1.000, 1.000, 1.000, 0.667, 0.667, 1.000, 0.667)
  )

  expect_true(is(alnbs[[3]], "BiologySeq"))

  expect_identical(length(alnbs), as.integer(3))

  expect_error(AAnumbering(alnbs) <- c("A", "B", "C"))

  AAnumbering(alnbs) <- c("A", "B", "C", "D", "E", "F", "G")

  expect_identical(
    AAnumbering(alnbs),
    c(
      "1" = "A", "2" = "B", "3" = "C", "4" = "D",
      "5" = "E", "6" = "F", "7" = "G"
    )
  )

  expect_identical(
    numbering(call_AAmutSet(alnbs)),
    c(
      "1" = "A", "2" = "B", "3" = "C", "4" = "D",
      "5" = "E", "6" = "F", "7" = "G"
    )
  )
})

test_that("BiologySeqMSA-2", {
  bss <- BiologySeqSet(
    c(
      "ATGCAGGTTACTGTG", "ATGCAGGTAACTGTG",
      "ATGCAGGTAAACCCTACTGAG", "ATGCAGGTAACTGAC"
    )
  )
  alnbs <- BiologySeqMSA(bss, order = "aligned", method = "ClustalOmega")

  expect_identical(
    as.character(AA(alnbs)),
    c(
      "3" = "MQVNPTE", "4" = "MQVTD--",
      "1" = "MQVTV--", "2" = "MQVTV--"
    )
  )
})
