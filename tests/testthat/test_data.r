context("Associations")
library(TwoSampleMR)
library(dplyr)
library(magrittr)

a <- extract_instruments(c(300, 301))
b <- extract_outcome_data(a$SNP, c(2,7))

test_that("extract_outcome_data is working", {
	expect_gt(nrow(b), 0)
})

ab <- harmonise_data(a,b)

test_that("harmonisation", {
	expect_equal(
		unique(paste(ab$id.exposure, ab$id.outcome)) %>% length,
		4
	)
})

