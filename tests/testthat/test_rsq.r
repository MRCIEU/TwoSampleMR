context("add rsq")

load(system.file("extdata", "test_commondata.RData", package="TwoSampleMR"))

test_that("exposure data 1", {
	# exp_dat <- extract_instruments("ieu-a-2")
	d <- exp_dat %>% add_rsq()
	expect_true("rsq.exposure" %in% names(d))
	expect_true("effective_n.exposure" %in% names(d))
})

test_that("exposure data 2", {
  skip("Skip unless you have good access to the API.")
  d <- extract_instruments(c('ieu-a-2', 'ieu-a-7')) %>% add_rsq()
	expect_true("rsq.exposure" %in% names(d))
	expect_true("effective_n.exposure" %in% names(d))
})

exposure <- exp_dat[1:5,]

test_that("outcome data 1", {
  skip("Skip unless you have good access to the API.")
  d <- extract_outcome_data(exposure$SNP, 'ieu-a-2') %>% add_rsq()
	expect_true("rsq.outcome" %in% names(d))
	expect_true("effective_n.outcome" %in% names(d))
})

test_that("outcome data 2", {
  skip("Skip unless you have good access to the API")
  d <- extract_outcome_data(exposure$SNP, c('ieu-a-2', 'ieu-a-7')) %>% add_rsq()
	expect_true("rsq.outcome" %in% names(d))
	expect_true("effective_n.outcome" %in% names(d))
})

test_that("dat 2", {
  skip("Skip unless you have good access to the API.")
  d <- make_dat(proxies=FALSE) %>% add_rsq()
	expect_true("rsq.outcome" %in% names(d) & "rsq.exposure" %in% names(d))
	expect_true("effective_n.outcome" %in% names(d) & "effective_n.exposure" %in% names(d))
})

test_that("dat ukb-d", {
  skip("Skip unless you have good access to the API.")
  skip_on_ci()
  skip_on_cran()
  d <- make_dat(exposure="ukb-d-30710_irnt", proxies=FALSE) %>% add_rsq()
	expect_true("rsq.outcome" %in% names(d) & "rsq.exposure" %in% names(d))
})

test_that("effective n", {
	effn <- effective_n(c(1000, 20000), c(49000, 30000))
	expect_true(
		effn[1] < effn[2]
	)
})

test_that("get_population_allele_frequency", {
  skip("Skip unless you have good access to the API.")  
	d <- extract_instruments("ieu-a-7")
	d <- add_metadata(d)
	d$eaf.exposure.controls <- get_population_allele_frequency(
		af = d$eaf.exposure,
		prop = d$ncase.exposure / (d$ncase.exposure + d$ncontrol.exposure),
		odds_ratio = exp(d$beta.exposure),
		prevalence = 0.2
	)
	expect_equal(cor(d$eaf.exposure, d$eaf.exposure.controls), 1, tolerance = 0.1)
})

test_that("bbj-a-1", {
  skip("Skip unless you have good access to the API.")
  skip_on_ci()
  skip_on_cran()
  d <- extract_instruments('bbj-a-1') %>% add_metadata() %>% add_rsq()
	expect_true(all(!is.na(d$rsq.exposure)))
})

test_that("bsen vs pn", {
	skip("Skip unless you have good access to the API.")
  d <- extract_instruments("ieu-a-2")
	r1 <- get_r_from_bsen(d$beta.exposure, d$se.exposure, d$samplesize.exposure)
	r2 <- get_r_from_pn(d$pval.exposure, d$samplesize.exposure)
	expect_true(cor(abs(r1), r2) > 0.99)
})
