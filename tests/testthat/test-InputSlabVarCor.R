context("Input UpdateSlabVar, MinSlabVar and MaxSlabVar of function for correlated version")

test_that("Throws warning if UpdateSlabVar parameter is not logical vector of length 1", {
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, UpdateSlabVar = 1), "UpdateSlabVar not provided as logical*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, UpdateSlabVar = "TRUE"), "UpdateSlabVar not provided as logical*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, UpdateSlabVar = 1:10), "UpdateSlabVar is not a vector of length 1*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, UpdateSlabVar = matrix(1,1,1)), "UpdateSlabVar is not a vector of length 1*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, UpdateSlabVar = c(TRUE, FALSE)), "UpdateSlabVar is not a vector of length 1*")
})

test_that("Throws warning if MinSlabVar parameter is not a numeric vector of length 1", {
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MinSlabVar = "A"), "MinSlabVar is not numeric*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MinSlabVar = TRUE), "MinSlabVar is not numeric*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MinSlabVar = 1:10), "MinSlabVar is not a vector of length 1*") 
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MinSlabVar = matrix(1,1,1)), "MinSlabVar is not a vector of length 1*")
})

test_that("Throws warning if MaxSlabVar parameter is not a numeric vector of length 1", {
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MaxSlabVar = "A"), "MaxSlabVar is not numeric*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MaxSlabVar = TRUE), "MaxSlabVar is not numeric*")
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MaxSlabVar = 1:10), "MaxSlabVar is not a vector of length 1*") 
  expect_warning(cpbayes_cor(1:10, 1:10, ExampleDataCor$cor, MaxSlabVar = matrix(1,1,1)), "MaxSlabVar is not a vector of length 1*")
})