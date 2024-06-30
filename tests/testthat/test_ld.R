context("ld")

skip_if_offline()
skip_if_offline(host = "api.opengwas.io")
skip_on_cran()

# extract some data
a <- try(extract_instruments("ieu-a-2", clump=FALSE))
if (class(a) == "try-error") skip("Server issues")
out <- try(clump_data(a))
if (class(out) == "try-error") skip("Server issues")


test_that("clump", {
  skip_if_not(exists('a'), "a not created in test above")
  skip_if_not(exists('out'), "out not created in test above")
	expect_equal(ncol(a), ncol(out))
	expect_true(nrow(a) > nrow(out))
	expect_true(nrow(out) > 0)
})


test_that("matrix", {
  skip_if_not(exists('out'), "out not created in test above")
	b <- try(ld_matrix(out$SNP))
	if (inherits(b, "try-error")) skip("Server issues")
	expect_equal(nrow(b), nrow(out))
	expect_equal(ncol(b), nrow(out))
})


test_that("clump multiple", {
	a <- try(extract_instruments(c("ieu-a-2", "ieu-a-1001"), clump=FALSE))
	if (class(a) == "try-error") skip("Server issues")
	out <- clump_data(a)
	expect_equal(length(unique(a$id.exposure)), length(unique(out$id.exposure)))
})

test_that("clump local", {
  skip_if_not(exists('a'), "a not created in test above")
	skip_if_not(file.exists("/Users/gh13047/repo/opengwas-api-internal/opengwas-api/app/ld_files/EUR.bim"))
	aclump <- clump_data(a, bfile="/Users/gh13047/repo/opengwas-api-internal/opengwas-api/app/ld_files/EUR", plink_bin="plink")
})
