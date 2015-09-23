library(airsea)
context("schmidt")

test_that("Sch returns expected values", {
              expect_equal(Sch('O2', 10, 35), 1033.831, tolerance = 0.00001)
              expect_equal(Sch('O2', 15, 5), 718.2619, tolerance = 0.00001)
              expect_equal(Sch('CO2', 20, 15), 638.3747, tolerance = 0.00001)
              expect_equal(Sch('CO2', 20, 35), 679.8753, tolerance = 0.00001)
              expect_equal(Sch('DMS', 10, 25), 1799.161, tolerance = 0.00001)
              expect_equal(Sch('acetone', 10, 35), 1865.888, tolerance = 0.00001)
              expect_equal(Sch('Ar', 10, 35), 948.7035, tolerance = 0.00001)
})

test_that("Sch method switching works", {
              expect_equal(Sch('CO2', 20, 35, method = 'mean'), 679.8753, tolerance = 1.0001)
              expect_equal(Sch('CO2', 20, 35, method = 'HM'), 727.8766, tolerance = 0.0001)
              expect_equal(Sch('CO2', 20, 35, method = 'WC'), 637.8135, tolerance = 0.0001)
})
