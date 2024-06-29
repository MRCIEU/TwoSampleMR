context("outcome")

test_that("outcomes", {

  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
  skip_on_cran()
  skip_on_ci()

	a <- extract_instruments("ieu-a-7")
	if(inherits(a, "response")) skip("Server issues")
	b <- extract_outcome_data(a$SNP, "ieu-a-2", proxies=FALSE)
	if(inherits(b, "response")) skip("Server issues")
	expect_true(nrow(b) < 30 & nrow(b) > 15)

	b <- extract_outcome_data(a$SNP, "ieu-a-2", proxies=TRUE)
	if(inherits(b, "response")) skip("Server issues")
	expect_true(nrow(b) > 30 & nrow(b) < nrow(a))

	b <- extract_outcome_data(a$SNP, c("ieu-a-2", "a"), proxies=FALSE)
	if(inherits(b, "response")) skip("Server issues")
	expect_true(nrow(b) < 30 & nrow(b) > 15)

	b <- extract_outcome_data(a$SNP, c("ieu-a-2", "a"), proxies=TRUE)
	if(inherits(b, "response")) skip("Server issues")
	expect_true(nrow(b) > 30 & nrow(b) < nrow(a))

	b <- extract_outcome_data(a$SNP, c("ieu-a-2", "ieu-a-7"), proxies=FALSE)
	if(inherits(b, "response")) skip("Server issues")
	expect_true(nrow(b) > 60)

	b <- extract_outcome_data(a$SNP, c("ieu-a-2", "ieu-a-7"), proxies=TRUE)
	if(inherits(b, "response")) skip("Server issues")
	expect_true(nrow(b) > 70)

})
