library(airsea)
context("Kw")

test_that("kw returns expected values", {
              expect_equal(kw('O2', 10, 9, 35, method = 'WA09', schmidt_method = "JS"), 3.795919e-05, tolerance = 1E-5)
              expect_equal(kw('CO2', 4, 10, 10, method = 'WA14', schmidt_method = "JS"), 4.492012e-05, tolerance = 1E-5)
              expect_equal(kw('DMS', 7, 7, 34, method = 'NG00', schmidt_method = "JS"), 1.902423e-05, tolerance = 1E-5)
})
