context("outcome")

skip_if_opengwas_tests_disabled()
skip_if_offline()
skip_if_offline(host = "api.opengwas.io")
skip_on_cran()

a <- try(extract_instruments("ieu-a-7"))
if (inherits(a, "try-error")) skip("Server issues")

test_that("outcomes 1", {
	b <- try(extract_outcome_data(a$SNP, "ieu-a-2", proxies=FALSE))
	if (inherits(b, "try-error")) skip("Server issues")
	expect_true(nrow(b) < 30 & nrow(b) > 15)
})

test_that("outcomes 2", {
	b <- try(extract_outcome_data(a$SNP, "ieu-a-2", proxies=TRUE))
	if (inherits(b, "try-error")) skip("Server issues")
	expect_true(nrow(b) > 30 & nrow(b) < nrow(a))
})

test_that("outcomes 3", {
	b <- try(extract_outcome_data(a$SNP, c("ieu-a-2", "a"), proxies=FALSE))
	if (inherits(b, "try-error")) skip("Server issues")
	expect_true(nrow(b) < 30 & nrow(b) > 15)
})

test_that("outcomes 4", {
	b <- try(extract_outcome_data(a$SNP, c("ieu-a-2", "a"), proxies=TRUE))
	if (inherits(b, "try-error")) skip("Server issues")
	expect_true(nrow(b) > 30 & nrow(b) < nrow(a))
})

test_that("outcomes 5", {
	b <- try(extract_outcome_data(a$SNP, c("ieu-a-2", "ieu-a-7"), proxies=FALSE))
	if (inherits(b, "try-error")) skip("Server issues")
	expect_true(nrow(b) > 60)
})

test_that("outcomes 6", {
	b <- try(extract_outcome_data(a$SNP, c("ieu-a-2", "ieu-a-7"), proxies=TRUE))
	if (inherits(b, "try-error")) skip("Server issues")
	expect_true(nrow(b) > 70)
})
