context("Primary Inputs of locFDR BF theoretic function cor")

test_that("Absent of BetaHat and SE will throw error", {
  expect_error(analytic_locFDR_BF_cor(1:10), "BetaHat or SE vector is missing!")
  expect_error(analytic_locFDR_BF_cor(SpikeVar=0.3), "BetaHat or SE vector is missing!")
  expect_error(analytic_locFDR_BF_cor(1:10, 1:10), "Correlation matrix is missing!")
})

test_that("Throws error if BetaHat is not a numeric vector of length more than 1 and no missing value", {
  expect_error(analytic_locFDR_BF_cor(1, 1:10, ExampleDataCor$cor), "Number of elements in the BetaHat vector must be more than 1!")
  expect_error(analytic_locFDR_BF_cor(NA, 1:10, ExampleDataCor$cor), "BetaHat must be a numeric vector.")
  expect_error(analytic_locFDR_BF_cor(c(NA,1), 1:10, ExampleDataCor$cor), "BetaHat for one or more phenotypes are missing!")
  expect_error(analytic_locFDR_BF_cor("AB", 1:10, ExampleDataCor$cor), "BetaHat must be a numeric vector.")
  expect_error(analytic_locFDR_BF_cor(c("A", 1), 1:10, ExampleDataCor$cor), "BetaHat must be a numeric vector.")
  expect_error(analytic_locFDR_BF_cor(matrix(0,2,2), 1:10, ExampleDataCor$cor), "BetaHat must be a vector.")
  expect_error(analytic_locFDR_BF_cor(data.frame(rep(1:10)), 1:10, ExampleDataCor$cor), "BetaHat must be a vector.")
  expect_warning(analytic_locFDR_BF_cor(matrix(1:10, 10, 1), 1:10, ExampleDataCor$cor), "BetaHat is a matrix!")
})

test_that("Throws error if SE is not a positive numeric vector of length more than 1 and no missing value", {
  expect_error(analytic_locFDR_BF_cor(1:10, 1, ExampleDataCor$cor), "Number of elements in the SE vector must be more than 1!")
  expect_error(analytic_locFDR_BF_cor(1:10, NA, ExampleDataCor$cor), "SE must be a numeric vector.")
  expect_error(analytic_locFDR_BF_cor(1:10, c(NA,1), ExampleDataCor$cor), "SE for one or more phenotypes are missing!")
  expect_error(analytic_locFDR_BF_cor(1:10, "AB", ExampleDataCor$cor), "SE must be a numeric vector.")
  expect_error(analytic_locFDR_BF_cor(1:10, c("A", 1), ExampleDataCor$cor), "SE must be a numeric vector.")
  expect_error(analytic_locFDR_BF_cor(1:10, matrix(0,2,2), ExampleDataCor$cor), "SE must be a vector.")
  expect_error(analytic_locFDR_BF_cor(1:10, data.frame(rep(1:10)), ExampleDataCor$cor), "SE must be a vector.")
  expect_error(analytic_locFDR_BF_cor(1:10, c(-1, 2,3), ExampleDataCor$cor), "One or more elements in the SE vector are not positive!")
  expect_warning(analytic_locFDR_BF_cor(1:10, matrix(1:10, 10, 1), ExampleDataCor$cor), "SE is a matrix!")
})

test_that("Throws error if Beta and SE are not of same length", {
  expect_error(analytic_locFDR_BF_cor(1:10, 1:15, ExampleDataCor$cor), "BetaHat and SE vectors must have the same number of elements!")
  expect_error(analytic_locFDR_BF_cor(1:100, 1:15, ExampleDataCor$cor), "BetaHat and SE vectors must have the same number of elements!")
})

