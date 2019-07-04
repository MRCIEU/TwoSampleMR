context("Retrieving and organising data")
library(TwoSampleMR)
library(dplyr)
library(magrittr)

test_that("extract_instruments is working", {
	expect_gt(nrow(extract_instruments(2)), 0)
	expect_equal(nrow(extract_instruments(2)), nrow(extract_instruments(2, force_server=TRUE)))
	expect_equal(
		extract_instruments(c(300, 301, 302)) %>%
			subset(id.exposure == 300) %>% nrow(), 
		extract_instruments(300, force_server=TRUE) %>% nrow()
	)
	expect_equal(
		extract_instruments(c(300, 301)) %$% unique(id.exposure) %>% length(),
		2
	)
})

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
