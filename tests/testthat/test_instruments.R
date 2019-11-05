context("Instruments")
library(TwoSampleMR)

test_that("server and mrinstruments", {

	# no no
	exp_dat <- extract_instruments(outcomes=c(1032))
	expect_true(length(unique(exp_dat$id)) == 0)

	# no yes
	exp_dat <- extract_instruments(outcomes=c(1282))
	expect_true(length(unique(exp_dat$id)) == 1)

	# yes no
	exp_dat <- extract_instruments(outcomes=c(2, 1032))
	expect_true(length(unique(exp_dat$id)) == 1)

	# yes yes
	exp_dat <- extract_instruments(outcomes=c(2, 1282))
	expect_true(length(unique(exp_dat$id)) == 2)

	exp_dat <- extract_instruments(outcomes=c(1032, 1282))
	expect_true(length(unique(exp_dat$id)) == 1)

	exp_dat <- extract_instruments(outcomes=c(2,100,1032,104,1,72,999))
	expect_true(length(unique(exp_dat$id)) == 6)

	exp_dat <- extract_instruments(outcomes=c(2,100,1032,104,1,72,999, 1282))
	expect_true(length(unique(exp_dat$id)) == 7)

})
