# Test mr()

# dat <- make_dat("ieu-a-2", "ieu-a-7")
load(system.file("extdata", "test_commondata.RData", package="TwoSampleMR"))

test_that("Test mr(): MR Egger, Weighted median, Inverse variance weighted, Simple mode, Weighted mode", {
  res <- mr(dat)
  expect_equal(nrow(res), 5L)
  expect_equal(ncol(res), 9L)
  expect_equal(res[1, "b"], 0.5025, tolerance = 1e-4)
  expect_equal(res[2, "b"], 0.3870, tolerance = 1e-4)
  expect_equal(res[3, "b"], 0.4459, tolerance = 1e-4)
  expect_equal(res[4, "b"], 0.3402, tolerance = 1e-4)
  expect_equal(res[5, "b"], 0.3791, tolerance = 1e-4)
})

test_that("mr.raps", {
  res2 <- suppressWarnings(mr(dat, method_list = "mr_raps"))
  expect_equal(nrow(res2), 1L)
  expect_equal(ncol(res2), 9L)
  expect_equal(res2[1, "b"], 0.4647, tolerance = 1e-4)
})
