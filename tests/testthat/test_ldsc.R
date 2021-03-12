context("LDSC")
library(TwoSampleMR)


test_that("get afl2", {
	skip("Very slow")
	hm3info <- ieugwasr::afl2_list("hapmap3")
	s <- ieugwasr::afl2_list()
	a <- ldsc_h2("ieu-a-2", snpinfo=hm3info)
	b <- ldsc_rg("ieu-a-2", "ieu-a-2", snpinfo=hm3info)
	c <- ldsc_rg("ukb-a-248", "ukb-b-19953", snpinfo=hm3info)
	height <- ldsc_h2("ukb-b-10787", snpinfo=hm3info)
})

