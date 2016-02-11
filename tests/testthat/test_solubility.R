library(airsea)
context("solubility")

test_that("Csat gives expected values", {
              expect_equal(Csat(10, 35)/44.6608, 6.315, tolerance = 0.0001)
})