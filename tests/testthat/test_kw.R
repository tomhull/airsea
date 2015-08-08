library(airsea)
context("Kw")

test_that("kw returns expected values", {
              expect_equal(kw('O2', 10, 9, 35, method = 'WA09'), 4.123147e-05, tolerance = 1E-6)
              expect_equal(kw('CO2', 4, 10, 10, method = 'WA09'), 3.829843e-05, tolerance = 1E-6)
              expect_equal(kw('DMS', 7, 7, 34, method = 'WA14'), 1.857818e-05, tolerance = 1E-6)
})
