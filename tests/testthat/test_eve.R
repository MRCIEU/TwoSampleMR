context("eve")

# dat <- make_dat("ieu-a-2", "ieu-a-7") %>% add_metadata()
load(system.file("extdata", "test_commondata.RData", package = "TwoSampleMR"))

test_that("wrapper", {
  skip_if_not_installed("car")
  expect_warning(w <- mr_wrapper(dat))
  expect_true(length(w) == 1)
  expect_true(length(names(w[[1]])) == 5)
})
