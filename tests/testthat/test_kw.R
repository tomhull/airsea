library(airsea)
context("Kw")

test_that("kw returns expected values", {
              expect_equal(kw('O2', 10, 9, 35, method = 'WA09'), 3.795919e-05, tolerance = 1E-6)
              expect_equal(kw('CO2', 4, 10, 10, method = 'WA14'), 4.492012e-05, tolerance = 1E-6)
              expect_equal(kw('DMS', 7, 7, 34, method = 'NG00'), 1.902423e-05, tolerance = 1E-6)
})
