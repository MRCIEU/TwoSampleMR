context("eve")
library(TwoSampleMR)

dat <- make_dat(2, 7)

test_that("wrapper", {
	w <- mr_wrapper(dat)
	expect_true(length(w) == 1)
	expect_true(length(names(w[[1]])) == 5)
})


