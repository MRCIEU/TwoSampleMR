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
