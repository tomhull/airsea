library(airsea)
context("tools")

test_that("convert_u10 works", {
              expect_equal(convert_u10(10, 15, method = 'liutang'), 9.477369, tolerance = 0.0000001)
              expect_error(convert_u10(10, 9), "formulations only work for z > 10")
})
