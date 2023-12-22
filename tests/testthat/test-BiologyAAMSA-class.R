test_that("BiologyAAMSA", {
  AA <- c("MQVNPTE", "MQVTE", "MQVTV")
  aln <- BiologyAAMSA(AA)

  expect_error(BiologyAAMSA(AA, skip_align = TRUE))

  expect_identical(
    as.character(AA(aln)),
    c("1" = "MQVNPTE", "2" = "MQV--TE", "3" = "MQV--TV")
  )

  expect_identical(toString(consAA(aln)), "MQV--TE")

  expect_identical(
    consAAfreq(aln),
    c(1.000, 1.000, 1.000, 0.667, 0.667, 1.000, 0.667)
  )

  expect_true(is(aln[[2]], "AAString"))

  expect_identical(length(aln), as.integer(3))

  expect_snapshot(call_AAmutSet(aln))
})
