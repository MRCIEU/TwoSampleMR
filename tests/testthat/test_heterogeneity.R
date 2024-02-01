context("heterogeneity")

dat <- make_dat()

test_that("heterogeneity", {
	w <- mr_heterogeneity(dat)
	expect_true(nrow(w) == 8)
})

test_that("egger intercept", {
	w <- mr_pleiotropy_test(dat)
	expect_true(nrow(w) == 4)
})
