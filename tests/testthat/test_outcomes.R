context("outcome")

test_that("outcomes", {

  skip("Skip unless you have good access to the API.")
  skip_on_ci()
  skip_on_cran()

	a <- extract_instruments("ieu-a-7")
	b <- extract_outcome_data(a$SNP, "ieu-a-2", proxies=FALSE)
	expect_true(nrow(b) < 30 & nrow(b) > 15)

	b <- extract_outcome_data(a$SNP, "ieu-a-2", proxies=TRUE)
	expect_true(nrow(b) > 30 & nrow(b) < nrow(a))

	b <- extract_outcome_data(a$SNP, c("ieu-a-2", "a"), proxies=FALSE)
	expect_true(nrow(b) < 30 & nrow(b) > 15)

	b <- extract_outcome_data(a$SNP, c("ieu-a-2", "a"), proxies=TRUE)
	expect_true(nrow(b) > 30 & nrow(b) < nrow(a))

	b <- extract_outcome_data(a$SNP, c("ieu-a-2", "ieu-a-7"), proxies=FALSE)
	expect_true(nrow(b) > 60)

	b <- extract_outcome_data(a$SNP, c("ieu-a-2", "ieu-a-7"), proxies=TRUE)
	expect_true(nrow(b) > 70)

})
