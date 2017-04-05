context("Input MCMCiter and Burnin of function for correlated version")

test_that("Throws warning if MCMCiter parameter is not a integer greater than a cutoff", {
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MCMCiter = 1), "MCMCiter should be at least 10000*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MCMCiter = "123"), "MCMCiter not provided as integer*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MCMCiter = 123.5), "MCMCiter not provided as integer*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MCMCiter = -1), "MCMCiter should be at least 10000*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MCMCiter = c(10000, 20000)), "MCMCiter is not a vector of length 1*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MCMCiter = as.matrix(c(10000, 20000))), "MCMCiter is not a vector of length 1*")
})

test_that("Throws warning if Burnin parameter is not a integer greater than a cutoff", {
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, Burnin = 1), "Burnin should be at least 5000*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, Burnin = "123"), "Burnin not provided as integer*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, Burnin = 123.5), "Burnin not provided as integer*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, Burnin = -1), "Burnin should be at least 5000*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, Burnin = c(10000, 20000)), "Burnin is not a vector of length 1*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, Burnin = as.matrix(c(10000, 20000))), "Burnin is not a vector of length 1*")
})

test_that("Throws warning if MCMC sample size is less than 5000", {
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, Burnin = 10000, MCMCiter = 10000), "*provided less than 5000*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, Burnin = 12000, MCMCiter = 10000), "*provided less than 5000*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, Burnin = 10000, MCMCiter = 12000), "*provided less than 5000*")
})

