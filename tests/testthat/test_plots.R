context("plots")

load(system.file("extdata", "test_commondata.RData", package="TwoSampleMR"))

test_that("MR scatter plot for mr_ivw", {
	# dat <- make_dat("ieu-a-2", "ieu-a-7")
	m <- mr(dat, method_list="mr_ivw")
	expect_no_error(p <- mr_scatter_plot(m, dat))
	expect_true(is.list(p))
	expect_true(length(p) == 1L)
})

test_that("Scatter plot for default set of estimates", {
  # dat <- make_dat("ieu-a-2", "ieu-a-7")
  m2 <- mr(dat)
  expect_no_error(p2 <- mr_scatter_plot(m2, dat))
  expect_true(is.list(p2))
  expect_true(length(p2) == 1L)
})

test_that("Scatter plot for mr_grip", {
  m3 <- mr(dat, method_list = "mr_grip")
  expect_no_error(p3 <- mr_scatter_plot(m3, dat))
  expect_true(is.list(p3))
  expect_true(length(p3) == 1L)
})

test_that("A second scatter plot for mr_grip", {
  m4 <- mr(dat2, method_list = "mr_grip")
  expect_no_error(p4 <- mr_scatter_plot(m4, dat2))
  expect_true(is.list(p4))
  expect_true(length(p4) == 4L)
})

res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)

test_that("Forest plot", {
  expect_no_error(p5 <- mr_forest_plot(res_single))
  expect_true(length(p5) == 1L)
})

test_that("Leave one out plot", {
  expect_no_error(p6 <- mr_leaveoneout_plot(res_loo))
  expect_true(length(p6) == 1L)
})

test_that("Funnel plot", {
  expect_no_error(p7 <- mr_funnel_plot(res_single))
  expect_true(length(p7) == 1L)
})
