context("ld")

test_that("extract some data", {
  skip("Skip unless you have good access to the API.")
  skip_on_ci()
  skip_on_cran()
  a <- extract_instruments(2, clump=FALSE)
  out <- clump_data(a)
})

test_that("clump", {
	skip_on_ci()
	skip_on_cran()
	expect_equal(ncol(a), ncol(out))
	expect_true(nrow(a) > nrow(out))
	expect_true(nrow(out) > 0)
})


test_that("matrix", {
	skip_on_ci()
	skip_on_cran()
	b <- ld_matrix(out$SNP)
	expect_equal(nrow(b), nrow(out))
	expect_equal(ncol(b), nrow(out))
})


test_that("clump multiple", {
	skip_on_ci()
	skip_on_cran()
	a <- extract_instruments(c("ieu-a-2", "ieu-a-1001"), clump=FALSE)
	out <- clump_data(a)
	expect_equal(length(unique(a$id.exposure)), length(unique(out$id.exposure)))
})

test_that("clump local", {
	skip_on_ci()
	skip_on_cran()
	skip_if_not(file.exists("/Users/gh13047/repo/opengwas-api-internal/opengwas-api/app/ld_files/EUR.bim"))
	aclump <- clump_data(a, bfile="/Users/gh13047/repo/opengwas-api-internal/opengwas-api/app/ld_files/EUR", plink_bin="plink")
})
