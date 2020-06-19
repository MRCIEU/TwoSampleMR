context("add rsq")
library(TwoSampleMR)

test_that("exposure data 1", {
	d <- extract_instruments('ieu-a-2') %>% add_rsq()
	expect_true("rsq.exposure" %in% names(d))
})

test_that("exposure data 2", {
	d <- extract_instruments(c('ieu-a-2', 'ieu-a-7')) %>% add_rsq()
	expect_true("rsq.exposure" %in% names(d))
})

exposure <- extract_instruments("ieu-a-2")[1:5,]

test_that("outcome data 1", {
	d <- extract_outcome_data(exposure$SNP, 'ieu-a-2') %>% add_rsq()
	expect_true("rsq.outcome" %in% names(d))
})

test_that("outcome data 2", {
	d <- extract_outcome_data(exposure$SNP, c('ieu-a-2', 'ieu-a-7')) %>% add_rsq()
	expect_true("rsq.outcome" %in% names(d))
})

test_that("dat 2", {
	d <- make_dat(proxies=FALSE) %>% add_rsq()
	expect_true("rsq.outcome" %in% names(d) & "rsq.exposure" %in% names(d))
})

test_that("dat ukb-d", {
	d <- make_dat(exposure="ukb-d-30710_irnt", proxies=FALSE) %>% add_rsq()
	expect_true("rsq.outcome" %in% names(d) & "rsq.exposure" %in% names(d))
})

