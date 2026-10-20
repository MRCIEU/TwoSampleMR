context("other formats")

# dat <- make_dat("ieu-a-2", "ieu-a-7")
load(system.file("extdata", "test_commondata.RData", package = "TwoSampleMR"))


test_that("MRInput", {
  skip_if(getRversion() < '4.4.0') # because MendelianRandomization will not be installed
  skip_if_not_installed("MendelianRandomization")
  w <- dat_to_MRInput(dat, get_correlations = FALSE)
  expect_true(length(w) == 1)
  expect_true(class(w) == "list")
  expect_true(class(w[[1]]) == "MRInput")
})

test_that("MRInput with cor", {
  skip_if_opengwas_tests_disabled()
  skip_on_cran()
  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
  skip_if(getRversion() < '4.4.0') # because MendelianRandomization will not be installed
  skip_if_not_installed("MendelianRandomization")
  w <- tryCatch(
    dat_to_MRInput(dat, get_correlations = TRUE)[[1]],
    error = skip("Server issues"),
    warning = expect_warning()
  )
  expect_true(nrow(w@correlation) == length(w@betaX))
})


test_that("harmonise_ld_dat handles LD matrix in a different SNP order", {
  x <- data.frame(
    SNP = c("rs1", "rs2", "rs3"),
    effect_allele.exposure = c("A", "C", "G"),
    other_allele.exposure = c("G", "T", "A"),
    beta.exposure = c(0.1, 0.2, 0.3),
    beta.outcome = c(0.01, 0.02, 0.03)
  )
  # LD matrix in a different order to x, with rs1 on the opposite effect allele
  ldnames <- c("rs3_G_A", "rs1_G_A", "rs2_C_T")
  ld <- matrix(
    c(1, 0.13, 0.23,
      0.13, 1, 0.12,
      0.23, 0.12, 1),
    nrow = 3,
    dimnames = list(ldnames, ldnames)
  )
  out <- harmonise_ld_dat(x, ld)
  expect_equal(out$x$SNP, c("rs1", "rs2", "rs3"))
  expect_equal(rownames(out$ld), c("rs1", "rs2", "rs3"))
  expect_equal(colnames(out$ld), c("rs1", "rs2", "rs3"))
  expect_equal(out$ld["rs1", "rs2"], 0.12)
  expect_equal(out$ld["rs1", "rs3"], 0.13)
  expect_equal(out$ld["rs2", "rs3"], 0.23)
  expect_equal(out$x$beta.exposure, c(-0.1, 0.2, 0.3))
  expect_equal(out$x$beta.outcome, c(-0.01, 0.02, 0.03))
  expect_equal(out$x$effect_allele.exposure, c("G", "C", "G"))
  expect_equal(out$x$other_allele.exposure, c("A", "T", "A"))
})


test_that("run_mrmix removes SNPs with eaf.exposure of 0 or 1", {
  skip_if_not_installed("MRMix")
  d <- dat
  d$eaf.exposure[1] <- 1
  d$eaf.exposure[2] <- 0
  expect_warning(
    expect_message(w <- run_mrmix(d)),
    "Removing 2 SNP\\(s\\)"
  )
  expect_message(w0 <- run_mrmix(dat[-(1:2), ]))
  expect_equal(w, w0)
})

test_that("run_mrmix returns NULL when no SNPs have a valid eaf.exposure", {
  skip_if_not_installed("MRMix")
  d <- dat
  d$eaf.exposure <- 1
  warnings <- character()
  expect_message(
    w <- withCallingHandlers(
      run_mrmix(d),
      warning = function(cnd) {
        warnings <<- c(warnings, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      }
    )
  )
  expect_match(warnings[1], "Removing 79 SNP\\(s\\)")
  expect_match(warnings[2], "No SNPs remaining")
  expect_null(w[[1]])
})


test_that("mrpresso", {
  w <- expect_warning(run_mr_presso(dat))
  expect_true(length(w) == 1)
  expect_true(class(w) == "list")
  expect_true(class(w[[1]]) == "list")
  expect_true(length(w[[1]]) == 2)
})


test_that("radial MR dat", {
  w <- dat_to_RadialMR(dat)
  expect_true(class(w) == "list")
  expect_true(nrow(w[[1]]) == nrow(dat))
})


test_that("radial MR", {
  w <- dat_to_RadialMR(dat)
  utils::capture.output(expect_warning(o <- RadialMR::ivw_radial(w[[1]])))
  expect_true(class(o) == "IVW")
})


test_that("radial MR", {
  utils::capture.output(
    w <- mr_ivw_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
  )
  expect_true(class(w) == "list")
})
