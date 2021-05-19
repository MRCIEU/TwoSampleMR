context("add metadata")
library(TwoSampleMR)

test_that("exposure data 1", {
	d <- extract_instruments('ieu-a-2') %>% add_metadata()
	expect_true("units.exposure" %in% names(d))
})

test_that("exposure data 2", {
	d <- extract_instruments(c('ieu-a-2', 'ieu-a-7')) %>% add_metadata()
	expect_true("units.exposure" %in% names(d))
})

exposure <- extract_instruments("ieu-a-2")[1:5,]

test_that("outcome data 1", {
	d <- extract_outcome_data(exposure$SNP, 'ieu-a-2') %>% add_metadata()
	expect_true("units.outcome" %in% names(d))
})

test_that("outcome data 2", {
	d <- extract_outcome_data(exposure$SNP, c('ieu-a-2', 'ieu-a-7')) %>% add_metadata()
	expect_true("units.outcome" %in% names(d))
})

test_that("dat 2", {
	d <- make_dat(proxies=FALSE) %>% add_metadata()
	expect_true("units.outcome" %in% names(d) & "units.exposure" %in% names(d))
})

test_that("no id1", {
	d <- extract_instruments("ieu-a-2")[1:5,]
	d$id.exposure <- "not a real id"
	d <- add_metadata(d)
	expect_true(!"units.exposure" %in% names(d))
})

test_that("no id2", {
	d <- extract_outcome_data(exposure$SNP, 'ieu-a-2')
	d$id.outcome <- "not a real id"
	d <- add_metadata(d)
	expect_true(!"units.outcome" %in% names(d))
})

test_that("ukb-d", {
	d <- extract_outcome_data(exposure$SNP, 'ukb-d-30710_irnt')
	d <- add_metadata(d)
	expect_true("units.outcome" %in% names(d))
})

test_that("bbj-a-1", {
	d <- extract_instruments("bbj-a-1")	%>% add_metadata()
	expect_true("samplesize.exposure" %in% names(d))
	expect_true(all(!is.na(d$samplesize.exposure)))
})

test_that("ieu-b-109", {
	expect_warning(d <- extract_instruments("ieu-b-109") %>% add_metadata())
	expect_true("samplesize.exposure" %in% names(d))
	expect_true(all(!is.na(d$samplesize.exposure)))
})

