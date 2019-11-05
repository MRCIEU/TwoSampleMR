context("steiger")
library(TwoSampleMR)

test_that("directionality", {
	w <- make_dat(2, 7)
	o <- directionality_test(w)
	expect_true(nrow(o) == 1)
})


test_that("directionality cc", {
	w <- make_dat(2, 7)
	w$r.outcome <- get_r_from_lor(w$beta.outcome, w$eaf.outcome, w$samplesize.outcome/2, w$samplesize.outcome/2, 0.1)
	o <- directionality_test(w)
	expect_true(nrow(o) == 1)
})


test_that("steiger filtering", {
	w <- make_dat(2, 7)
	w <- steiger_filtering(w)
	expect_true("steiger_pval" %in% names(w))
})
