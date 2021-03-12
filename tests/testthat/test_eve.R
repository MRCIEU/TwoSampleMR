context("eve")
library(TwoSampleMR)

dat <- make_dat("ieu-a-2", "ieu-a-7") %>% add_metadata()

test_that("wrapper", {
	expect_warning(w <- mr_wrapper(dat))
	expect_true(length(w) == 1)
	expect_true(length(names(w[[1]])) == 5)
})


