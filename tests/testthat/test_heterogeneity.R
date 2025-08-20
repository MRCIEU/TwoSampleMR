context("heterogeneity")

# dat2 <- make_dat()
load(system.file("extdata", "test_commondata.RData", package = "TwoSampleMR"))

test_that("heterogeneity", {
  w <- mr_heterogeneity(dat2)
  expect_true(nrow(w) == 8)
})

test_that("egger intercept", {
  w <- mr_pleiotropy_test(dat2)
  expect_true(nrow(w) == 4)
})
