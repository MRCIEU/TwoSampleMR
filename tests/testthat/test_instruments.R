context("Instruments")
library(TwoSampleMR)

test_that("server and mrinstruments", {

	# no no
	exp_dat <- extract_instruments(outcomes=c("ieu-a-1032"))
	expect_true(length(unique(exp_dat$id)) == 0)

	# no yes
	exp_dat <- extract_instruments(outcomes=c("ebi-a-GCST004634"))
	expect_true(length(unique(exp_dat$id)) == 1)

	# yes no
	exp_dat <- extract_instruments(outcomes=c("ieu-a-2", "ieu-a-1032"))
	expect_true(length(unique(exp_dat$id)) == 1)

	# yes yes
	exp_dat <- extract_instruments(outcomes=c("ieu-a-2", "ebi-a-GCST004634"))
	expect_true(length(unique(exp_dat$id)) == 2)

	exp_dat <- extract_instruments(outcomes=c("ieu-a-1032", "ebi-a-GCST004634"))
	expect_true(length(unique(exp_dat$id)) == 1)

	exp_dat <- extract_instruments(outcomes=c(2,100,"ieu-a-1032",104,72,999))
	expect_true(length(unique(exp_dat$id)) == 5)

	exp_dat <- extract_instruments(outcomes=c(2,100,"ieu-a-1032",104,72,999, "ebi-a-GCST004634"))
	expect_true(length(unique(exp_dat$id)) == 6)
})


test_that("read data", {
	exp_dat <- extract_instruments("ieu-a-2")
	names(exp_dat) <- gsub(".exposure", "", names(exp_dat))
	fn <- tempfile()
	write.table(exp_dat, file=fn, row=F, col=T, qu=F, sep="\t")

	a <- read_exposure_data(fn, sep="\t")
	b <- read_outcome_data(fn, sep="\t")

	expect_true("chr.outcome" %in% names(b))
	expect_true("chr.exposure" %in% names(a))
})