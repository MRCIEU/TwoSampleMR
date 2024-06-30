context("plots")

load(system.file("extdata", "test_commondata.RData", package="TwoSampleMR"))

test_that("scatter plot", {
	# dat <- make_dat("ieu-a-2", "ieu-a-7")
	m <- mr(dat, method_list="mr_ivw")
	p <- mr_scatter_plot(m, dat)
	expect_true(is.list(p))
	expect_true(length(p) == 1)
})
