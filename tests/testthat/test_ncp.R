library(airsea)
context("ncp")

dat = data.frame(timePeriod = 90000, C0 = 282, C1 = 292, T0 = 10, T1 = 10.5, S0 = 35, S1 = 34.5, h0 = 25, h1 = 45, Pslp0 = 1000, Pslp1 = 998, u0 = 7, u1 = 5, Cb0 = 282, Cb1 = 282, entrainment = T)

test_that("O2NCP.mean returns expected values", {
              expect_equal(O2NCP.mean(dat), 479.1553, tolerance = 0.00001)
})