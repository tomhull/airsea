library(airsea)
context("schmidt")

test_that("Calculation of molar volume at boiling point is correct", {
  expect_equal(Vb("O2"), 21)
  expect_equal(Vb("CO2"), 35)
  expect_equal(Vb("MeOH"), 42.7)
})

test_that("Seawater viscosity", {
  expect_equal(n_sw(20, 35), 1.071764, tolerance = 0.00001)
  expect_equal(n_sw(10, 15), 1.340564, tolerance = 0.00001)
  expect_equal(v_sw(20, 35), 0.01046451, tolerance = 0.00001)
})

test_that("Sch returns expected values for various compounds", {
              expect_equal(Sch('O2', 10, 35, method = "JS"), 876.2457, tolerance = 0.00001)
              expect_equal(Sch('O2', 15, 5, method = "JS"), 614.8143, tolerance = 0.00001)
              expect_equal(Sch('CO2', 20, 15, method = "JS"), 638.3747, tolerance = 0.00001)
              expect_equal(Sch('CO2', 20, 35, method = "JS"), 679.8753, tolerance = 0.00001)
              expect_equal(Sch('DMS', 10, 25, method = "JS"), 1799.161, tolerance = 0.00001)
              expect_equal(Sch('Ar', 10, 35, method = "JS"), 948.7035, tolerance = 0.00001)
})

test_that("Wanninkof check values", {
  expect_equal(Sch("Ar", 20, 35, method = "WA"), 615, tolerance = 0.1)
  expect_equal(Sch("O2", 20, 35, method = "WA"), 568, tolerance = 0.1)
  expect_equal(Sch("DMS", 20, 35, method = "WA"), 941, tolerance = 0.1)
})

test_that("Sch method switching works", {
              expect_equal(Sch('CO2', 20, 35, method = 'JS'), 679.8753, tolerance = 0.0001)
              expect_equal(Sch('CO2', 20, 35, method = 'HM'), 727.8766, tolerance = 0.0001)
              expect_equal(Sch('CO2', 20, 35, method = 'WC'), 637.8135, tolerance = 0.0001)
              expect_equal(Sch('CO2', 20, 35, method = 'HL'), 705.9441, tolerance = 0.0001)
})

# Vb for O2?
# 28.04
# 21 (predicted from kcalcs)