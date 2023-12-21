test_that("BiologyAAmutSet", {
  mut_list <- list(
    mut1 = c("D123G", "D6F", "C878C", "C878T"),
    mut2 = c("D123G", "D6F"),
    mut3 = c("C878C", "D123G", "C878T")
  )

  muts <- BiologyAAmutSet(mut_list)


  expect_true(is(muts[3], "BiologyAAmutSet"))
  expect_true(is(muts[[3]], "BiologyAAmut"))

  expect_identical(
    sort(muts),
    BiologyAAmutSet(
      list(
        mut1 = c("D6F", "D123G", "C878C", "C878T"),
        mut2 = c("D6F", "D123G"),
        mut3 = c("D123G", "C878C", "C878T")
      )
    )
  )

  expect_identical(
    select_mut(muts, 1, 130),
    BiologyAAmutSet(list(
      mut1 = c("D123G", "D6F"),
      mut2 = c("D123G", "D6F"),
      mut3 = c("D123G")
    ))
  )

  expect_identical(
    unique(muts, method = "site"),
    BiologyAAmutSet(
      list(
        mut1 = c("D123G", "D6F", "C878C"),
        mut2 = c("D123G", "D6F"),
        mut3 = c("C878C", "D123G")
      )
    )
  )

  expect_identical(
    unique(BiologyAAmutSet(list("D12C", "D12C"))),
    BiologyAAmutSet(list("D12C"))
  )

  names(muts)[3] <- "3"
  expect_identical(names(muts), c("mut1", "mut2", "3"))
})