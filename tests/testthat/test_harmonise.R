context("harmonise")
library(TwoSampleMR)

exp_dat <- extract_instruments("ieu-a-2")
out_dat <- extract_outcome_data(exp_dat$SNP, "ieu-a-7")

test_that("check columns before harmonising", {
	expect_null(check_required_columns(exp_dat, "exposure"))
	expect_null(check_required_columns(out_dat, "outcome"))
	expect_error(check_required_columns(subset(exp_dat, select=-id.exposure), "exposure"))
})

