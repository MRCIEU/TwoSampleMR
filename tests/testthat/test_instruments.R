context("Instruments")

skip_if_offline()
skip_if_offline(host = "api.opengwas.io")
skip_on_cran()
skip_on_ci()

test_that("server and mrinstruments 1", {
	# no no
  exp_dat <- try(extract_instruments(outcomes=c("ieu-a-1032")))
  if (class(exp_dat) == "try-error") skip("Server issues")
	expect_true(length(unique(exp_dat$id)) == 0)
})


test_that("server and mrinstruments 2", {
	# no yes
	exp_dat <- try(extract_instruments(outcomes=c("ebi-a-GCST004634")))
	if (class(exp_dat) == "try-error") skip("Server issues")
	expect_true(length(unique(exp_dat$id)) == 1)
})
	
test_that("server and mrinstruments 3", {
	# yes no
	exp_dat <- try(extract_instruments(outcomes=c("ieu-a-2", "ieu-a-1032")))
	if (class(exp_dat) == "try-error") skip("Server issues")
	expect_true(length(unique(exp_dat$id)) == 1)
})

test_that("server and mrinstruments 4", {
	# yes yes
	exp_dat <- try(extract_instruments(outcomes=c("ieu-a-2", "ebi-a-GCST004634")))
	if (class(exp_dat) == "try-error") skip("Server issues")
	expect_true(length(unique(exp_dat$id)) == 2)
})

test_that("server and mrinstruments 5", {
	exp_dat <- try(extract_instruments(outcomes=c("ieu-a-1032", "ebi-a-GCST004634")))
	if (class(exp_dat) == "try-error") skip("Server issues")
	expect_true(length(unique(exp_dat$id)) == 1)
})

test_that("server and mrinstruments 6", {
	exp_dat <- try(extract_instruments(outcomes=c(2,100,"ieu-a-1032",104,72,999)))
	if (class(exp_dat) == "try-error") skip("Server issues")
	expect_true(length(unique(exp_dat$id)) == 5)
})

test_that("server and mrinstruments 7", {
	exp_dat <- try(extract_instruments(outcomes=c(2,100,"ieu-a-1032",104,72,999, "ebi-a-GCST004634")))
	if (class(exp_dat) == "try-error") skip("Server issues")
	expect_true(length(unique(exp_dat$id)) == 6)
})

load(system.file("extdata", "test_commondata.RData", package="TwoSampleMR"))

test_that("read data", {
	exp_dat <- try(extract_instruments("ieu-a-2"))
	if (class(exp_dat) == "try-error") skip("Server issues")
	names(exp_dat) <- gsub(".exposure", "", names(exp_dat))
	fn <- tempfile()
	write.table(exp_dat, file=fn, row=FALSE, col=TRUE, qu=FALSE, sep="\t")

	a <- read_exposure_data(fn, sep="\t")
	b <- read_outcome_data(fn, sep="\t")

	expect_true("chr.outcome" %in% names(b))
	expect_true("chr.exposure" %in% names(a))
})
