context("Tophits")
library(TwoSampleMR)
library(dplyr)
library(magrittr)


test_that("extract_instruments is working", {
	expect_gt(nrow(extract_instruments(2)), 0)
	expect_equal(nrow(extract_instruments(2)), nrow(extract_instruments("IEU-a-2", force_server=TRUE)))
	expect_equal(
		extract_instruments(c(300, 301, 302)) %>%
			subset(id.exposure == 301) %>% nrow(), 
		extract_instruments("IEU-a-301", force_server=TRUE) %>% nrow()
	)
	expect_equal(
		extract_instruments(c(300, 301)) %$% unique(id.exposure) %>% length(),
		2
	)
})

