context("add rsq")

load(system.file("extdata", "test_commondata.RData", package="TwoSampleMR"))

test_that("exposure data 1", {
	# exp_dat <- extract_instruments("ieu-a-2")
	d <- exp_dat %>% add_rsq()
	expect_true("rsq.exposure" %in% names(d))
	expect_true("effective_n.exposure" %in% names(d))
})

test_that("exposure data 2", {
  skip_if_opengwas_tests_disabled()
  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
  skip_on_cran()
  d <- try(extract_instruments(c('ieu-a-2', 'ieu-a-7')))
  if(inherits(d, "try-error")) skip("Server issues")
  d <- d %>% add_rsq()
	expect_true("rsq.exposure" %in% names(d))
	expect_true("effective_n.exposure" %in% names(d))
})

exposure <- exp_dat[1:5,]

test_that("outcome data 1", {
  skip_if_opengwas_tests_disabled()
  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
  skip_on_cran()
  skip_on_ci()
  d <- try(extract_outcome_data(exposure$SNP, 'ieu-a-2'))
  if(inherits(d, "try-error")) skip("Server issues")
  d <- try(d %>% add_rsq())
  if (inherits(d, "try-error")) skip("Server issues")
	expect_true("rsq.outcome" %in% names(d))
	expect_true("effective_n.outcome" %in% names(d))
})

test_that("outcome data 2", {
  skip_if_opengwas_tests_disabled()
  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
  skip_on_cran()
  skip_on_ci()
  d <- try(extract_outcome_data(exposure$SNP, c('ieu-a-2', 'ieu-a-7')))
  if(inherits(d, "try-error")) skip("Server issues")
  d <- try(d %>% add_rsq())
  if(inherits(d, "try-error")) skip("Server issues")
	expect_true("rsq.outcome" %in% names(d))
	expect_true("effective_n.outcome" %in% names(d))
})

test_that("dat 2", {
  skip_if_opengwas_tests_disabled()
  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
  skip_on_cran()
  d <- try(make_dat(proxies=FALSE))
  if(inherits(d, "try-error")) skip("Server issues")
  d <- d %>% add_rsq()
	expect_true("rsq.outcome" %in% names(d) & "rsq.exposure" %in% names(d))
	expect_true("effective_n.outcome" %in% names(d) & "effective_n.exposure" %in% names(d))
})

test_that("dat ukb-d", {
  skip_if_opengwas_tests_disabled()
  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
  skip_on_cran()
  d <- try(make_dat(exposure="ukb-d-30710_irnt", proxies=FALSE))
  if(inherits(d, "try-error")) skip("Server issues")
  d <- d %>% add_rsq()
	expect_true("rsq.outcome" %in% names(d) & "rsq.exposure" %in% names(d))
})

test_that("effective n", {
	effn <- effective_n(c(1000, 20000), c(49000, 30000))
	expect_true(
		effn[1] < effn[2]
	)
})

test_that("get_population_allele_frequency", {
  skip_if_opengwas_tests_disabled()
  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
  skip_on_cran()
	d <- try(extract_instruments("ieu-a-7"))
	if(inherits(d, "try-error")) skip("Server issues")
	d <- try(add_metadata(d))
	if(inherits(d, "try-error")) skip("Server issues")
	d$eaf.exposure.controls <- get_population_allele_frequency(
		af = d$eaf.exposure,
		prop = d$ncase.exposure / (d$ncase.exposure + d$ncontrol.exposure),
		odds_ratio = exp(d$beta.exposure),
		prevalence = 0.2
	)
	expect_equal(cor(d$eaf.exposure, d$eaf.exposure.controls), 1, tolerance = 0.1)
})

test_that("bbj-a-1", {
  skip_if_opengwas_tests_disabled()
  skip_on_cran()
  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
  d <- try(extract_instruments('bbj-a-1'))
  if(inherits(d, "try-error")) skip("Server issues")
  d <- try(d %>% add_metadata() %>% add_rsq())
  if(inherits(d, "try-error")) skip("Server issues")
	expect_true(all(!is.na(d$rsq.exposure)))
})

test_that("bsen vs pn", {
  skip_if_opengwas_tests_disabled()
  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
  skip_on_cran()
  d <- try(extract_instruments("ieu-a-2"))
  if(inherits(d, "try-error")) skip("Server issues")
	r1 <- get_r_from_bsen(d$beta.exposure, d$se.exposure, d$samplesize.exposure)
	r2 <- get_r_from_pn(d$pval.exposure, d$samplesize.exposure)
	expect_true(cor(abs(r1), r2) > 0.99)
})
