context("Input GenVar and Traits of function for Uncorrelated version")

test_that("Throws error if the argument Phenotypes is present but not a character vector of same length as BetaHat with no duplicate entry", {
  expect_error(cpbayes_uncor(1:10, 1:10, as.matrix(c("phen1", "phen2"))), "Phenotypes must be a vector.")
  expect_error(cpbayes_uncor(1:10, 1:10, as.data.frame(c("phen1", "phen2"))), "Phenotypes must be a vector.")
  expect_error(cpbayes_uncor(1:10, 1:10, 1:10), "Phenotypes must be a character vector.")
  expect_error(cpbayes_uncor(1:10, 1:10, rep("a",10)), "Two or more phenotypes have the same name!")
  expect_error(cpbayes_uncor(1:10, 1:10, c("A", "B", "C")), "BetaHat and Phenotypes vectors must have the same number of elements!")
})

test_that("Throws error if the argument Variant is present but not a character vector of 1", {
  expect_error(cpbayes_uncor(1:10, 1:10, Variant =as.matrix("A")), "Variant must be a vector.")
  expect_error(cpbayes_uncor(1:10, 1:10, Variant =as.data.frame("A")), "Variant must be a vector.")
  expect_warning(cpbayes_uncor(1:10, 1:10, Variant =1), "Variant is not a character vactor!")
  expect_error(cpbayes_uncor(1:10, 1:10, Variant =1:2), "Variant must be a vector of length 1.")
  expect_warning(cpbayes_uncor(1:10, 1:10, Variant =NA), "Variant is NA!")
})