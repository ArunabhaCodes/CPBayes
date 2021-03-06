context("Input MCMCiter and Burnin of function for uncorrelated version")

test_that("Throws warning if MCMCiter parameter is not a integer greater than a cutoff", {
  skip_on_cran()
  expect_warning(cpbayes_uncor(1:10, 1:10, MCMCiter = 1), "MCMCiter should be at least 2200*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MCMCiter = "123"), "MCMCiter not provided as integer*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MCMCiter = 123.5), "MCMCiter not provided as integer*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MCMCiter = -1), "MCMCiter should be at least 2200*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MCMCiter = c(10000, 20000)), "MCMCiter is not a scalar*")
  expect_warning(cpbayes_uncor(1:10, 1:10, MCMCiter = as.matrix(c(10000, 20000))), "MCMCiter is not a scalar*")
})

test_that("Throws warning if Burnin parameter is not a integer greater than a cutoff", {
  skip_on_cran()
  expect_warning(cpbayes_uncor(1:10, 1:10, Burnin = 1), "Burnin should be at least 200*")
  expect_warning(cpbayes_uncor(1:10, 1:10, Burnin = "123"), "Burnin not provided as integer*")
  expect_warning(cpbayes_uncor(1:10, 1:10, Burnin = 123.5), "Burnin not provided as integer*")
  expect_warning(cpbayes_uncor(1:10, 1:10, Burnin = -1), "Burnin should be at least 200*")
  expect_warning(cpbayes_uncor(1:10, 1:10, Burnin = c(10000, 20000)), "Burnin is not a scalar*")
  expect_warning(cpbayes_uncor(1:10, 1:10, Burnin = as.matrix(c(10000, 20000))), "Burnin is not a scalar*")
})

test_that("Throws warning if MCMC sample size is less than 5000", {
  skip_on_cran()
  expect_warning(cpbayes_uncor(1:10, 1:10, Burnin = 10000, MCMCiter = 10000), "*provided less than 2000*")
  expect_warning(cpbayes_uncor(1:10, 1:10, Burnin = 12000, MCMCiter = 10000), "*provided less than 2000*")
  expect_warning(cpbayes_uncor(1:10, 1:10, Burnin = 10000, MCMCiter = 11000), "*provided less than 2000*")
})

