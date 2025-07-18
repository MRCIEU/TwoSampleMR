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

# forest plot 1 to many
rm(list = ls())
load(system.file("extdata", "forestplot_1_to_many_data.RData", package="TwoSampleMR"))

test_that("Forest plot 1 to many", {
  expect_warning(p8 <- forest_plot_1_to_many(
    res,
    b = "b",
    se = "se",
    exponentiate = TRUE,
    ao_slc = FALSE,
    lo = 0.3,
    up = 2.5,
    TraitM = "exposure",
    col1_width = 2,
    by = NULL,
    trans = "log2",
    xlab = "OR for CHD per SD increase in risk factor (95% confidence interval)",
    weight = "weight"
  ), regexp = "Removed 6 rows")
})

test_that("Forest plot 1 to many test 2", {
  res$pval<-formatC(res$pval, format = "e", digits = 2)
  expect_warning(p9 <- forest_plot_1_to_many(
    res,
    b = "b",
    se = "se",
    exponentiate = TRUE,
    ao_slc = FALSE,
    lo = 0.3,
    up = 2.5,
    TraitM = "exposure",
    by = NULL,
    trans = "log2",
    xlab = "OR for CHD per SD increase in risk factor (95% CI)",
    weight = "weight",
    subheading_size = 11,
    col1_title = "Risk factor",
    col1_width = 2.5,
    col_text_size = 4,
    addcols = c("nsnp", "pval"),
    addcol_widths = c(1.0, 1.0),
    addcol_titles = c("No. SNPs", "P-val")
  ), regexp = "Removed 6 rows")
})

test_that("Forest plot 1 to many test 3 - with subcategory in by argument", {
  res <- mr(dat2)
  res <- split_exposure(res)
  res <- subset_on_method(res)
  res$subcategory[res$exposure %in% c("Adiponectin", "Hip circumference", "Waist circumference")] <- "Group 1"
  res$subcategory[is.na(res$subcategory)] <- "Group 2"
  res$weight <- 1/res$se
  res <- sort_1_to_many(res, sort_action = 1, group = "subcategory")
  expect_warning(p10 <- forest_plot_1_to_many(
    res,
    b = "b",
    se = "se",
    exponentiate = TRUE,
    trans = "log2",
    ao_slc = FALSE,
    lo = 0.3,
    up = 2.5,
    TraitM = "exposure",
    col_text_size = 4,
    col1_width = 1.5,
    by = "subcategory",
    xlab = "OR for CHD per SD increase in risk factor (95% confidence interval)",
    subheading_size = 14,
    weight = "weight"
  ), regexp = "Removed 3 rows")
})
