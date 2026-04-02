context("singlesnp")

# dat <- make_dat("ieu-a-2", "ieu-a-7")
load(system.file("extdata", "test_commondata.RData", package = "TwoSampleMR"))

test_that("singlesnp", {
  w <- mr_singlesnp(dat)
  expect_true(nrow(w) == sum(dat$mr_keep) + 2)
})

test_that("singlesnp_plot", {
  w <- mr_singlesnp(dat)
  p <- mr_forest_plot(w)
  expect_true(length(p) == 1)
})

test_that("density plot", {
  w <- mr_singlesnp(dat)
  m <- mr(dat, method_list = "mr_ivw")
  p <- mr_density_plot(w, m)
  expect_true(length(p) == 1)
})

test_that("singlesnp_plot", {
  w <- mr_singlesnp(dat)
  p <- mr_funnel_plot(w)
  expect_true(length(p) == 1)
})

test_that("forest_plot_with_category", {
  w <- mr_singlesnp(dat)
  snp_rows <- !grepl("^All", w$SNP)
  w$category <- NA_character_
  w$category[snp_rows] <- sample(
    c("Cluster 1", "Cluster 2", "Junk"),
    size = sum(snp_rows),
    replace = TRUE
  )
  p <- mr_forest_plot(w)
  expect_true(length(p) == 1)
  expect_s3_class(p[[1]], "ggplot")
})

rm(list = ls())

test_that("RadialMR outliers on forest plot with categories", {
  load(system.file("extdata", "test_commondata.RData", package = "TwoSampleMR"))

  # Subset to 20 SNPs for a cleaner plot
  set.seed(42)
  snps_keep <- sample(unique(dat$SNP), 20)
  dat20 <- dat[dat$SNP %in% snps_keep, ]

  res_single <- mr_singlesnp(dat20)

  # Run RadialMR to detect outliers
  radial_dat <- dat_to_RadialMR(dat20)
  invisible(utils::capture.output(
    radial_res <- suppressWarnings(RadialMR::ivw_radial(radial_dat[[1]], alpha = 0.05, weights = 3))
  ))
  outlier_snps <- radial_res$outliers$SNP

  # Label SNPs as Outlier or Main
  snp_rows <- !grepl("^All", res_single$SNP)
  res_single$category <- NA_character_
  res_single$category[snp_rows] <- ifelse(
    res_single$SNP[snp_rows] %in% outlier_snps, "Outlier", "Main"
  )

  p <- mr_forest_plot(res_single)
  expect_true(length(p) == 1)
  expect_s3_class(p[[1]], "ggplot")
})


