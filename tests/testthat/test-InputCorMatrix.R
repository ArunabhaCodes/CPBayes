context("Input Correlation matrix of function for correlated version")

test_that("Throws error if the Corln is not a square matrix or Corln matrix has any missing entry", {
  expect_error(cpbayes_cor(1:10, 1:10, as.data.frame(ExampleDataCor$cor)), "Corln must be a matrix not a data.frame.*")
  expect_error(cpbayes_cor(1:10, 1:10, 1:10), "Corln must be a matrix.")
  expect_error(cpbayes_cor(1:10, 1:10, as.matrix(c("a", "b", "c", "d"))), "Corln must be a numeric matrix.")
  expect_error(cpbayes_cor(1:10, 1:10, matrix(c(1:99, NA), 10, 10)), "One or more entries of Corln are missing!")
  expect_error(cpbayes_cor(1:10, 1:10, matrix(1:100, 25, 4)), "Number of rows and columns of Corln are different!")
})

test_that("Throws error or warning if the Corln is not a symmetric positive definite matrix", {
  expect_error(cpbayes_cor(1:10, 1:10, matrix(1:100, 10, 10)), "Corln is not symmetric!")
  expect_error(cpbayes_cor(1:10, 1:10, matrix(1:25, 5, 5)), "Number of rows of Corln and length of BetaHat do not match!")
  expect_error(cpbayes_cor(1:10, 1:10, diag(c(-1, 2:10))), "Corln is negative definite!")
  expect_error(cpbayes_cor(1:10, 1:10, diag(-c(1:10))), "Diagonal elements of Corln are not 1!")
  expect_warning(cpbayes_cor(1:10, 1:10, matrix(rep(1, 100), 10, 10)), "Corln is a singular matrix!")
})