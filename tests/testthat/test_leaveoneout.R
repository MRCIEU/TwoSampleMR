context("leaveoneout")

# dat <- make_dat("ieu-a-2", "ieu-a-7")
load(system.file("extdata", "test_commondata.RData", package = "TwoSampleMR"))

test_that("leaveoneout", {
  w <- mr_leaveoneout(dat)
  expect_true(nrow(w) == sum(dat$mr_keep) + 1)
})

test_that("leaveoneout_plot", {
  w <- mr_leaveoneout(dat)
  p <- mr_leaveoneout_plot(w)
  expect_true(length(p) == 1)
})
