context("plots")

load(system.file("extdata", "test_commondata.RData", package="TwoSampleMR"))

test_that("MR scatter plot for mr_ivw", {
	# dat <- make_dat("ieu-a-2", "ieu-a-7")
	m <- mr(dat, method_list="mr_ivw")
	p <- mr_scatter_plot(m, dat)
	expect_true(is.list(p))
	expect_true(length(p) == 1L)
})

test_that("Scatter plot for default set of estimates", {
  # dat <- make_dat("ieu-a-2", "ieu-a-7")
  m2 <- mr(dat)
  p2 <- mr_scatter_plot(m2, dat)
  expect_true(is.list(p2))
  expect_true(length(p2) == 1L)
})
