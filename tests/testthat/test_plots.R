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

test_that("Scatter plot for mr_grip", {
  m3 <- mr(dat, method_list = "mr_grip")
  p3 <- mr_scatter_plot(m3, dat)
  expect_true(is.list(p3))
  expect_true(length(p3) == 1L)
})

test_that("A second scatter plot for mr_grip", {
  m4 <- mr(dat2, method_list = "mr_grip")
  p4 <- mr_scatter_plot(m4, dat2)
  expect_true(is.list(p4))
  expect_true(length(p4) == 4L)
})
