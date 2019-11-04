context("outcome")
library(TwoSampleMR)

test_that("outcomes", {
	
	a <- extract_instruments(7)
	b <- extract_outcome_data(a$SNP, 2, proxies=FALSE)
	expect_true(nrow(b) < 30 & nrow(b) > 15)

	b <- extract_outcome_data(a$SNP, 2, proxies=TRUE)
	expect_true(nrow(b) > 30 & nrow(b) < nrow(a))

	b <- extract_outcome_data(a$SNP, c(2, "a"), proxies=FALSE)
	expect_true(nrow(b) < 30 & nrow(b) > 15)

	b <- extract_outcome_data(a$SNP, c(2, "a"), proxies=TRUE)
	expect_true(nrow(b) > 30 & nrow(b) < nrow(a))

	b <- extract_outcome_data(a$SNP, c(2, 7), proxies=FALSE)
	expect_true(nrow(b) > 60)

	b <- extract_outcome_data(a$SNP, c(2, 7), proxies=TRUE)
	expect_true(nrow(b) > 70)

})
