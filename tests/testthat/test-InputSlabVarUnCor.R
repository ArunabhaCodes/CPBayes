context("Input UpdateSlabVar, MinSlabVar and MaxSlabVar of function for uncorrelated version")

test_that("Throws warning if UpdateSlabVar parameter is not logical scalar", {
  skip_on_cran()
  expect_warning(cpbayes_uncor(1:10, 1:10, UpdateSlabVar = 1), "UpdateSlabVar not provided as logical*")
  expect_warning(cpbayes_uncor(1:10, 1:10, UpdateSlabVar = "TRUE"), "UpdateSlabVar not provided as logical*")
  expect_warning(cpbayes_uncor(1:10, 1:10, UpdateSlabVar = 1:10), "UpdateSlabVar is not a scalar*")
  expect_warning(cpbayes_uncor(1:10, 1:10, UpdateSlabVar = matrix(1,1,1)), "UpdateSlabVar is not a scalar*")
  expect_warning(cpbayes_uncor(1:10, 1:10, UpdateSlabVar = c(TRUE, FALSE)), "UpdateSlabVar is not a scalar*")
})

test_that("Throws warning if MinSlabVar parameter is not a numeric scalar", {
  skip_on_cran()
  expect_warning(cpbayes_uncor(1:10, 1:10, MinSlabVar = "A"), "MinSlabVar is not numeric*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MinSlabVar = TRUE), "MinSlabVar is not numeric*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MinSlabVar = 1:10), "MinSlabVar is not a scalar*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MinSlabVar = matrix(1,1,1)), "MinSlabVar is not a scalar*")
})

test_that("Throws warning if MaxSlabVar parameter is not a numeric scalar", {
  skip_on_cran()
  expect_warning(cpbayes_uncor(1:10, 1:10, MaxSlabVar = "A"), "MaxSlabVar is not numeric*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MaxSlabVar = TRUE), "MaxSlabVar is not numeric*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MaxSlabVar = 1:10), "MaxSlabVar is not a scalar*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MaxSlabVar = matrix(1,1,1)), "MaxSlabVar is not a scalar*")
})
