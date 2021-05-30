context("mvmr")
library(TwoSampleMR)

test_that("control", {
	lipids <- mv_extract_exposures(c("ieu-a-299","ieu-a-300","ieu-a-302"))
	chd <- extract_outcome_data(lipids$SNP, "ieu-a-7")
	control <- mv_harmonise_data(lipids, chd)
	mv_residual(control, intercept=TRUE, instrument_specific=TRUE)$result
	mv_residual(control, intercept=FALSE, instrument_specific=TRUE)$result
	mv_residual(control, intercept=TRUE, instrument_specific=FALSE)$result
	mv_residual(control, intercept=FALSE, instrument_specific=FALSE)$result
	mv_multiple(control, intercept=TRUE, instrument_specific=TRUE)$result
	mv_multiple(control, intercept=FALSE, instrument_specific=TRUE)$result
	mv_multiple(control, intercept=TRUE, instrument_specific=FALSE)$result
	mv_multiple(control, intercept=FALSE, instrument_specific=FALSE)$result
})


test_that("dat", {
	a <- mv_extract_exposures(c("ukb-b-5238", "ieu-a-1001"))
	b <- extract_outcome_data(a$SNP, "ieu-a-297")
	dat <- mv_harmonise_data(a, b)
	mv_residual(dat, intercept=TRUE, instrument_specific=TRUE)$result
	mv_residual(dat, intercept=FALSE, instrument_specific=TRUE)$result
	mv_residual(dat, intercept=TRUE, instrument_specific=FALSE)$result
	mv_residual(dat, intercept=FALSE, instrument_specific=FALSE)$result
	mv_multiple(dat, intercept=TRUE, instrument_specific=TRUE)$result
	mv_multiple(dat, intercept=FALSE, instrument_specific=TRUE)$result
	mv_multiple(dat, intercept=TRUE, instrument_specific=FALSE)$result
	mv_multiple(dat, intercept=FALSE, instrument_specific=FALSE)$result
	mv_ivw(dat)$result
	mv_basic(dat)$result
})


test_that("ordering", {
	lipids1 <- mv_extract_exposures(c("ieu-a-299","ieu-a-300","ieu-a-302"))
	chd1 <- extract_outcome_data(lipids1$SNP, "ieu-a-7")
	control1 <- mv_harmonise_data(lipids1, chd1)
	mv_multiple(control1)

	lipids2 <- mv_extract_exposures(c("ieu-a-302","ieu-a-300","ieu-a-299"))
	chd2 <- extract_outcome_data(lipids2$SNP, "ieu-a-7")
	control2 <- mv_harmonise_data(lipids2, chd2)
	mv_multiple(control2)

})

