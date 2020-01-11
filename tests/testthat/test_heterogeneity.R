context("heterogeneity")
library(TwoSampleMR)


test_that("heterogeneity", {
	dat <- make_dat()
	w <- mr_heterogeneity(dat)
	expect_true(nrow(w) == 8)
})

test_that("egger intercept", {
	dat <- make_dat()
	w <- mr_pleiotropy_test(dat)
	expect_true(nrow(w) == 4)
})

