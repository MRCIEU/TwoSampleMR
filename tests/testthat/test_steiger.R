context("steiger")

w <- make_dat(2, 7)

test_that("directionality", {
	o <- directionality_test(w)
	expect_true(nrow(o) == 1)
})


test_that("directionality cc", {
	w$r.outcome <- get_r_from_lor(w$beta.outcome, w$eaf.outcome, w$samplesize.outcome/2, w$samplesize.outcome/2, 0.1)
	o <- directionality_test(w)
	expect_true(nrow(o) == 1)
})


test_that("steiger filtering", {
	w <- steiger_filtering(w)
	expect_true("steiger_pval" %in% names(w))
})
