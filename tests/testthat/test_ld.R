context("ld")

skip_if_offline()
skip_if_offline(host = "api.opengwas.io")
skip_on_cran()
skip_on_ci()

test_that("extract some data", {
  a <- extract_instruments(2, clump=FALSE)
  if(length(a) == 0) skip("Server issues")
  out <- clump_data(a)
})

test_that("clump", {
  skip_if_not(exists('a'), "a not created in test above")
  skip_if_not(exists('out'), "out not created in test above")
	expect_equal(ncol(a), ncol(out))
	expect_true(nrow(a) > nrow(out))
	expect_true(nrow(out) > 0)
})


test_that("matrix", {
  skip_if_not(exists('out'), "out not created in test above")
	b <- ld_matrix(out$SNP)
	expect_equal(nrow(b), nrow(out))
	expect_equal(ncol(b), nrow(out))
})


test_that("clump multiple", {
	a <- extract_instruments(c("ieu-a-2", "ieu-a-1001"), clump=FALSE)
	out <- clump_data(a)
	expect_equal(length(unique(a$id.exposure)), length(unique(out$id.exposure)))
})

test_that("clump local", {
  skip_if_not(exists('a'), "a not created in test above")
	skip_if_not(file.exists("/Users/gh13047/repo/opengwas-api-internal/opengwas-api/app/ld_files/EUR.bim"))
	aclump <- clump_data(a, bfile="/Users/gh13047/repo/opengwas-api-internal/opengwas-api/app/ld_files/EUR", plink_bin="plink")
})
