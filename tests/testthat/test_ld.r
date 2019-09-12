context("LD functions")
library(TwoSampleMR)


a <- extract_instruments(2)

test_that("ld ref", {
	expect_equal(
		nrow(a), nrow(clump_data(a))
	)
})


test_that("ld matrix", {
	expect_equal(
		length(unique(a$SNP)), nrow(ld_matrix(a$SNP))
	)
})

