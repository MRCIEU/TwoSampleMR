context("other formats")
library(TwoSampleMR)

dat <- make_dat("ieu-a-2", "ieu-a-7")

test_that("MRInput", {
	w <- dat_to_MRInput(dat, get_correlations=FALSE)
	expect_true(length(w) == 1)	
	expect_true(class(w) == "list")
	expect_true(class(w[[1]]) == "MRInput")
})

test_that("MRInput with cor", {
	expect_warning(w <- dat_to_MRInput(dat, get_correlations=TRUE)[[1]])
	expect_true(nrow(w@correlation) == length(w@betaX))
})


test_that("mrpresso", {
	expect_warning(w <- run_mr_presso(dat))
	expect_true(length(w) == 1)	
	expect_true(class(w) == "list")
	expect_true(class(w[[1]]) == "list")
	expect_true(length(w[[1]]) == 2)
})


test_that("radial MR dat", {
	w <- dat_to_RadialMR(dat)
	expect_true(class(w) == "list")
	expect_true(nrow(w[[1]]) == nrow(dat))
})


test_that("radial MR", {
	w <- dat_to_RadialMR(dat)
	expect_warning(o <- RadialMR::ivw_radial(w[[1]]))
	expect_true(class(o) == "IVW")
})


test_that("radial MR", {
	w <- mr_ivw_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
	expect_true(class(w) == "list")
})

