# Test mr()

# dat <- make_dat("ieu-a-2", "ieu-a-7")
load(system.file("extdata", "test_commondata.RData", package="TwoSampleMR"))

test_that("Test mr(): MR Egger, Weighted median, Inverse variance weighted, Simple mode, Weighted mode", {
  res <- mr(dat)
  expect_equal(nrow(res), 5L)
  expect_equal(ncol(res), 9L)
  expect_equal(res[1, "b"], 0.5025, tolerance = 1e-3)
  expect_equal(res[2, "b"], 0.3870, tolerance = 1e-3)
  expect_equal(res[3, "b"], 0.4459, tolerance = 1e-3)
  expect_equal(res[4, "b"], 0.3402, tolerance = 1e-3)
  expect_equal(res[5, "b"], 0.3791, tolerance = 1e-1)
})

test_that("mr.raps", {
  skip_if_not_installed("mr.raps")
  res2 <- suppressWarnings(mr(dat, method_list = "mr_raps"))
  expect_equal(nrow(res2), 1L)
  expect_equal(ncol(res2), 9L)
  expect_equal(res2[1, "b"], 0.4647, tolerance = 1e-3)
})

test_that("mr.raps over.dispersion option", {
  skip_if_not_installed("mr.raps")
  params <- default_parameters()
  params$over.dispersion <- FALSE
  res3 <- suppressWarnings(mr(dat, method_list = "mr_raps", parameters = params))
  expect_equal(nrow(res3), 1L)
  expect_equal(ncol(res3), 9L)
  expect_equal(res3[1, "b"], 0.4682, tolerance = 1e-3)
})

test_that("mr.raps loss.function option", {
  skip_if_not_installed("mr.raps")
  params <- default_parameters()
  params$loss.function <- "tukey"
  res4 <- suppressWarnings(mr(dat, method_list = "mr_raps", parameters = params))
  expect_equal(nrow(res4), 1L)
  expect_equal(ncol(res4), 9L)
  expect_equal(res4[1, "b"], 0.4788, tolerance = 1e-3)
})

test_that("mr.raps shrinkage option", {
  skip_if_not_installed("mr.raps")
  params <- default_parameters()
  params$shrinkage <- TRUE
  res5 <- suppressWarnings(mr(dat, method_list = "mr_raps", parameters = params))
  expect_equal(nrow(res5), 1L)
  expect_equal(ncol(res5), 9L)
  expect_equal(res5[1, "b"], 0.4647, tolerance = 1e-3)
})

test_that("mr_grip()", {
  res6 <- suppressWarnings(mr(dat, method_list = "mr_grip"))
  expect_equal(nrow(res6), 1L)
  expect_equal(ncol(res6), 9L)
  expect_equal(res6[1, "b"], 0.490, tolerance = 1e-3)
})

test_that("mr_grip() not from mr()", {
  tst <- mr_grip(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
  expect_equal(tst$b, 0.49, tol = 1e-4)
  expect_equal(tst$se, 0.0947, tol = 1e-4)
  expect_equal(tst$b_i, -3.38e-5, tol = 1e-6)
  expect_equal(tst$se_i, 5.66e-5, tol = 1e-6)
  expect_equal(tst$b.adj, 0.44736, tol = 1e-4)
  expect_equal(tst$se.adj, 0.16389, tol = 1e-4)
  expect_equal(tst$nsnp, 79L)
})
