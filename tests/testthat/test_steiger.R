context("steiger")

# dat <- make_dat("ieu-a-2", "ieu-a-7")
load(system.file("extdata", "test_commondata.RData", package="TwoSampleMR"))

test_that("directionality", {
	o <- directionality_test(dat)
	expect_true(nrow(o) == 1)
})

test_that("directionality cc", {
	dat$r.outcome <- get_r_from_lor(dat$beta.outcome, dat$eaf.outcome, dat$samplesize.outcome/2, dat$samplesize.outcome/2, 0.1)
	o <- directionality_test(dat)
	expect_true(nrow(o) == 1)
})

test_that("steiger filtering", {
	expect_warning(dat <- steiger_filtering(dat))
	expect_true("steiger_pval" %in% names(dat))
})


test_that("steiger sensitivity uconf", {
	r <- steiger_sensitivity_conf(bxy=0.1, bxyo=0.2, bgx=0.5, vx=1, vy=1, vg=0.5, beta_a=1, beta_b=10)
	expect_equal(r$prop, 1)
})
