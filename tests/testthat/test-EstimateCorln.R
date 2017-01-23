context("Inputs of estimate_corln function")

test_that("Throws error if inputs are missing", {
  expect_error(estimate_corln(), "n11, n00 or n10 matrix is missing!")
  expect_error(estimate_corln(SampleOverlapMatrix$n11), "n11, n00 or n10 matrix is missing!")
  expect_error(estimate_corln(SampleOverlapMatrix$n11, SampleOverlapMatrix$n00), "n11, n00 or n10 matrix is missing!")
})

test_that("Throws error if inputs are not matrices", {
  expect_error(estimate_corln(as.data.frame(SampleOverlapMatrix$n11), SampleOverlapMatrix$n00, SampleOverlapMatrix$n10),
    "n11 must be a matrix not a data.frame.*")
  expect_error(estimate_corln(SampleOverlapMatrix$n11, as.data.frame(SampleOverlapMatrix$n00), SampleOverlapMatrix$n10),
    "n00 must be a matrix not a data.frame.*")
  expect_error(estimate_corln(SampleOverlapMatrix$n11, SampleOverlapMatrix$n00, as.data.frame(SampleOverlapMatrix$n10)),
    "n10 must be a matrix not a data.frame.*")
  expect_error(estimate_corln(1:10, 1, 1), "n11 must be a matrix.")
  expect_error(estimate_corln(SampleOverlapMatrix$n11, 1, 1), "n00 must be a matrix.")
  expect_error(estimate_corln(SampleOverlapMatrix$n11, SampleOverlapMatrix$n00, 1), "n10 must be a matrix.")
})

test_that("Throws error if inputs are not square matrices", {
  expect_error(estimate_corln(matrix(0, 2,3), SampleOverlapMatrix$n11, SampleOverlapMatrix$n00), 
    "Number of rows and columns of n11 are different!")
  expect_error(estimate_corln(SampleOverlapMatrix$n11, matrix(0, 2,3), SampleOverlapMatrix$n00), 
    "Number of rows and columns of n00 are different!")
  expect_error(estimate_corln(SampleOverlapMatrix$n11, SampleOverlapMatrix$n00, matrix(0, 2,3)), 
    "Number of rows and columns of n10 are different!")
})

test_that("Throws error if inputs are not numeric matrices with no missing entries", {
  expect_error(estimate_corln(matrix("a", 5, 5), SampleOverlapMatrix$n00, SampleOverlapMatrix$n10), 
    "n11 must be a numeric matrix.")
  expect_error(estimate_corln(SampleOverlapMatrix$n00, matrix("a", 5, 5), SampleOverlapMatrix$n10), 
    "n00 must be a numeric matrix.")
  expect_error(estimate_corln(SampleOverlapMatrix$n00, SampleOverlapMatrix$n10, matrix("a", 5, 5)), 
    "n10 must be a numeric matrix.")
  expect_error(estimate_corln(matrix(NA, 5, 5), SampleOverlapMatrix$n00, SampleOverlapMatrix$n10), 
    "One or more entries of n11 are missing!")
  expect_error(estimate_corln(SampleOverlapMatrix$n00, matrix(NA, 5, 5), SampleOverlapMatrix$n10), 
    "One or more entries of n00 are missing!")
  expect_error(estimate_corln(SampleOverlapMatrix$n00, SampleOverlapMatrix$n10, matrix(c(1:24, NA), 5, 5)), 
    "One or more entries of n10 are missing!")
})

test_that("Throws error if inputs are not non-negative integer matrices", {
  expect_error(estimate_corln(SampleOverlapMatrix$n00, SampleOverlapMatrix$n10, matrix(c(1:24, 1.5), 5, 5)),
    "Every element of n10 must be a non-negative integer.")
  expect_error(estimate_corln(SampleOverlapMatrix$n00, SampleOverlapMatrix$n10, matrix(c(1:24, -1), 5, 5)),
    "Every element of n10 must be a non-negative integer.")
  expect_error(estimate_corln(SampleOverlapMatrix$n00, matrix(c(1:24, 1.5), 5, 5), SampleOverlapMatrix$n10),
    "Every element of n00 must be a non-negative integer.")
  expect_error(estimate_corln(SampleOverlapMatrix$n00, matrix(c(1:24, -1), 5, 5), SampleOverlapMatrix$n10),
    "Every element of n00 must be a non-negative integer.")
  expect_error(estimate_corln(matrix(c(1:24, 1.5), 5, 5), SampleOverlapMatrix$n00, SampleOverlapMatrix$n10),
    "Every element of n11 must be a non-negative integer.")
  expect_error(estimate_corln(matrix(c(1:24, -1), 5, 5), SampleOverlapMatrix$n00, SampleOverlapMatrix$n10),
    "Every element of n11 must be a non-negative integer.")
})

test_that("Throws error if input matrices are not all of same size", {
  expect_error(estimate_corln(SampleOverlapMatrix$n00, SampleOverlapMatrix$n10, matrix(1:16,4,4)),
    "n11, n00, n10 matrices must have same dimension.")
  expect_error(estimate_corln(SampleOverlapMatrix$n00, matrix(1:16,4,4), SampleOverlapMatrix$n10),
    "n11, n00, n10 matrices must have same dimension.")
  expect_error(estimate_corln(matrix(1:16,4,4), SampleOverlapMatrix$n00, SampleOverlapMatrix$n10),
    "n11, n00, n10 matrices must have same dimension.")
})

test_that("Throws error if n11 and n00 are not symmetric matrices", {
  expect_error(estimate_corln(SampleOverlapMatrix$n00, SampleOverlapMatrix$n10, SampleOverlapMatrix$n00), 
    "n00 must be symmetric!")
  expect_error(estimate_corln(SampleOverlapMatrix$n10, SampleOverlapMatrix$n00, SampleOverlapMatrix$n00), 
    "n11 must be symmetric!")
})

test_that("Throws error if any diagonal entry of n00 and n11 matrices is zero", {
  expect_error(estimate_corln(SampleOverlapMatrix$n00, matrix(1, 5, 5)-diag(c(1,rep(0, 4))), SampleOverlapMatrix$n10),
    "Diagonal elements of n00 must be positive integer.")
  expect_error(estimate_corln(matrix(1, 5, 5)-diag(c(1,rep(0, 4))), SampleOverlapMatrix$n00, SampleOverlapMatrix$n10),
    "Diagonal elements of n11 must be positive integer.")
})

test_that("Throws error if any diagonal entry of n10 matrices is not zero", {
  expect_error(estimate_corln(SampleOverlapMatrix$n00, SampleOverlapMatrix$n11, matrix(1, 5, 5)-diag(c(1,rep(0, 4)))),
    "Diagonal elements of n10 must be zero.")
  expect_error(estimate_corln(SampleOverlapMatrix$n00, SampleOverlapMatrix$n11, matrix(1, 5, 5)-diag(c(0,rep(1, 4)))),
    "Diagonal elements of n10 must be zero.")
})

