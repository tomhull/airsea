library(airsea)
context("Schmidt")

test_that("Sch returns expected values", {
              expect_equal(Sch('O2', 10, 35), 876.2457, tolerance = 0.00001)
              expect_equal(Sch('O2', 15, 5), 614.8143, tolerance = 0.00001)
              expect_equal(Sch('CO2', 20, 15), 638.3747, tolerance = 0.00001)
              expect_equal(Sch('DMS', 10, 25), 1799.161, tolerance = 0.00001)
              expect_equal(Sch('acetone', 10, 35), 1865.888, tolerance = 0.00001)
              expect_equal(Sch('Ar', 10, 35), 948.7035, tolerance = 0.00001)
})

test_that("Sch method switching works", {
              expect_equal(Sch('CO2', 10, 35, method = 'mean'), 1172.965, tolerance = 0.0001)
              expect_equal(Sch('CO2', 10, 35, method = 'HM'), 1241.108, tolerance = 0.0001)
              expect_equal(Sch('CO2', 10, 35, method = 'WC'), 1111.916, tolerance = 0.0001)
})
