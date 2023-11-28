context("plots")

test_that("scatter plot", {
	dat <- make_dat(2,7)
	m <- mr(dat, method_list="mr_ivw")
	p <- mr_scatter_plot(m, dat)
	expect_true(is.list(p))
	expect_true(length(p) == 1)
})
