context("add metadata")

# get required data
# d1 <- extract_instruments('ieu-a-2')
# d2 <- extract_instruments(c('ieu-a-2', 'ieu-a-7'))
# exposure <- extract_instruments("ieu-a-2")[1:5,]
# d3 <- extract_outcome_data(exposure$SNP, 'ieu-a-2')
# d4 <- extract_outcome_data(exposure$SNP, c('ieu-a-2', 'ieu-a-7'))
# d5 <- make_dat(proxies=FALSE)
# d6 <- extract_instruments("ieu-a-2")[1:5,]
# d7 <- extract_outcome_data(exposure$SNP, 'ieu-a-2')
# d8 <- extract_outcome_data(exposure$SNP, 'ukb-d-30710_irnt')
# d9 <- extract_instruments("bbj-a-1") 
# d10 <- extract_instruments("ieu-b-109")

# save(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, file="inst/extdata/test_add_metadata.RData")

load(system.file("extdata", "test_add_metadata.RData", package="TwoSampleMR"))

test_that("exposure data 1", {
	d1 <- d1 %>% add_metadata()
	expect_true("units.exposure" %in% names(d1))
})

test_that("exposure data 2", {
	d2 <- d2 %>% add_metadata()
	expect_true("units.exposure" %in% names(d2))
})

test_that("outcome data 1", {
	d3 <- d3 %>% add_metadata()
	expect_true("units.outcome" %in% names(d3))
})

test_that("outcome data 2", {
	d4 <- d4 %>% add_metadata()
	expect_true("units.outcome" %in% names(d4))
})

test_that("dat 2", {
	d5 <- d5 %>% add_metadata()
	expect_true("units.outcome" %in% names(d5) & "units.exposure" %in% names(d5))
})

test_that("no id1", {
	d6$id.exposure <- "not a real id"
	d6 <- add_metadata(d6)
	expect_true(!"units.exposure" %in% names(d6))
})

test_that("no id2", {
	d7$id.outcome <- "not a real id"
	d7 <- add_metadata(d7)
	expect_true(!"units.outcome" %in% names(d7))
})

test_that("ukb-d", {
	d8 <- add_metadata(d8)
	expect_true("units.outcome" %in% names(d8))
})

test_that("bbj-a-1", {
	d9 <- d9 %>% add_metadata()
	expect_true("samplesize.exposure" %in% names(d9))
	expect_true(all(!is.na(d9$samplesize.exposure)))
})

test_that("ieu-b-109", {
	d10 <- d10 %>% add_metadata()
	expect_true("samplesize.exposure" %in% names(d10))
	expect_true(all(!is.na(d10$samplesize.exposure)))
})

