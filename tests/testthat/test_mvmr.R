context("mvmr")

skip_on_ci()
skip_on_cran()
skip_if_offline()
skip_if_offline(host = "api.opengwas.io")

test_that("control", {
	lipids <- mv_extract_exposures(c("ieu-a-299","ieu-a-300","ieu-a-302"))
	chd <- extract_outcome_data(lipids$SNP, "ieu-a-7")
	control <- mv_harmonise_data(lipids, chd)
	expect_output(print(mv_residual(control, intercept=TRUE, instrument_specific=TRUE)))
	expect_output(print(mv_residual(control, intercept=FALSE, instrument_specific=TRUE)))
	expect_output(print(mv_residual(control, intercept=TRUE, instrument_specific=FALSE)))
	expect_output(print(mv_residual(control, intercept=FALSE, instrument_specific=FALSE)))
	expect_output(print(mv_multiple(control, intercept=TRUE, instrument_specific=TRUE)))
	expect_output(print(mv_multiple(control, intercept=FALSE, instrument_specific=TRUE)))
	expect_output(print(mv_multiple(control, intercept=TRUE, instrument_specific=FALSE)))
	expect_output(print(mv_multiple(control, intercept=FALSE, instrument_specific=FALSE)))
})


test_that("dat", {
	a <- mv_extract_exposures(c("ukb-b-5238", "ieu-a-1001"))
	b <- extract_outcome_data(a$SNP, "ieu-a-297")
	dat <- mv_harmonise_data(a, b)
	expect_output(print(mv_residual(dat, intercept=TRUE, instrument_specific=TRUE)))
	expect_output(print(mv_residual(dat, intercept=FALSE, instrument_specific=TRUE)))
	expect_output(print(mv_residual(dat, intercept=TRUE, instrument_specific=FALSE)))
	expect_output(print(mv_residual(dat, intercept=FALSE, instrument_specific=FALSE)))
	expect_output(print(mv_multiple(dat, intercept=TRUE, instrument_specific=TRUE)))
	expect_output(print(mv_multiple(dat, intercept=FALSE, instrument_specific=TRUE)))
	expect_output(print(mv_multiple(dat, intercept=TRUE, instrument_specific=FALSE)))
	expect_output(print(mv_multiple(dat, intercept=FALSE, instrument_specific=FALSE)))
	expect_output(print(mv_ivw(dat)))
	expect_output(print(mv_basic(dat)))
})


test_that("ordering 1", {
	lipids1 <- mv_extract_exposures(c("ieu-a-299","ieu-a-300","ieu-a-302"))
	chd1 <- extract_outcome_data(lipids1$SNP, "ieu-a-7")
	control1 <- mv_harmonise_data(lipids1, chd1)
	expect_output(print(mv_multiple(control1)))
})

test_that("ordering 2", {
	lipids2 <- mv_extract_exposures(c("ieu-a-302","ieu-a-300","ieu-a-299"))
	chd2 <- extract_outcome_data(lipids2$SNP, "ieu-a-7")
	control2 <- mv_harmonise_data(lipids2, chd2)
	expect_output(print(mv_multiple(control2)))
})
