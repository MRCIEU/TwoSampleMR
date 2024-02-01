context("ld")


a <- extract_instruments(2, clump=FALSE)
out <- clump_data(a)


test_that("clump", {
	expect_equal(ncol(a), ncol(out))
	expect_true(nrow(a) > nrow(out))
	expect_true(nrow(out) > 0)
})

b <- ld_matrix(out$SNP)

test_that("matrix", {
	expect_equal(nrow(b), nrow(out))
	expect_equal(ncol(b), nrow(out))
})



a <- extract_instruments(c("ieu-a-2", "ieu-a-1001"), clump=FALSE)
out <- clump_data(a)

test_that("clump multiple", {
	expect_equal(length(unique(a$id.exposure)), length(unique(out$id.exposure)))
})
