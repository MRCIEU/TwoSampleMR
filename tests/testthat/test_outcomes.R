context("outcome")
library(TwoSampleMR)

test_that("outcomes", {
	
	expect_warning(a <- extract_instruments(7))
	expect_warning(b <- extract_outcome_data(a$SNP, 2, proxies=FALSE))
	expect_true(nrow(b) < 30 & nrow(b) > 15)

	expect_warning(b <- extract_outcome_data(a$SNP, 2, proxies=TRUE))
	expect_true(nrow(b) > 30 & nrow(b) < nrow(a))

	expect_warning(b <- extract_outcome_data(a$SNP, c(2, "a"), proxies=FALSE))
	expect_true(nrow(b) < 30 & nrow(b) > 15)

	expect_warning(b <- extract_outcome_data(a$SNP, c(2, "a"), proxies=TRUE))
	expect_true(nrow(b) > 30 & nrow(b) < nrow(a))

	expect_warning(b <- extract_outcome_data(a$SNP, c(2, 7), proxies=FALSE))
	expect_true(nrow(b) > 60)

	expect_warning(b <- extract_outcome_data(a$SNP, c(2, 7), proxies=TRUE))
	expect_true(nrow(b) > 70)

})
