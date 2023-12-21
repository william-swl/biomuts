test_that("BiologyAAmut", {
  expect_error(BiologyAAmut(c("C12A", "T18FC")))
  expect_error(BiologyAAmut(c("AC12A", "T18FC")))

  mut <- BiologyAAmut(c("D123G", "D6F", "C878C", "D123G", "C878T"))

  expect_identical(
    unique(sort(mut)),
    BiologyAAmut(c("D6F", "D123G", "C878C", "C878T"))
  )
  expect_identical(
    unique(sort(mut), bysite = TRUE),
    BiologyAAmut(c("D6F", "D123G", "C878C"))
  )

  expect_identical(
    select_mut(mut, 1, 150),
    BiologyAAmut(c("D123G", "D6F", "D123G"))
  )
})
