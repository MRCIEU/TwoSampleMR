make_mvdat <- function() {
  set.seed(1234)
  nsnp <- 40
  exposure_beta <- matrix(stats::rnorm(nsnp * 2, sd = 0.1), ncol = 2)
  colnames(exposure_beta) <- c("exp1", "exp2")
  # 30 SNPs pass 5e-8 for exp1 and 10 for exp2
  exposure_pval <- cbind(
    exp1 = c(rep(1e-10, 30), rep(0.01, 10)),
    exp2 = c(rep(0.01, 30), rep(1e-10, 10))
  )
  outcome_beta <- as.vector(exposure_beta %*% c(0.5, -0.3)) + stats::rnorm(nsnp, sd = 0.01)
  list(
    exposure_beta = exposure_beta,
    exposure_pval = exposure_pval,
    exposure_se = matrix(0.01, nsnp, 2, dimnames = list(NULL, c("exp1", "exp2"))),
    outcome_beta = outcome_beta,
    outcome_se = rep(0.01, nsnp),
    outcome_pval = rep(0.5, nsnp),
    expname = data.frame(id.exposure = c("exp1", "exp2"), exposure = c("Exposure 1", "Exposure 2")),
    outname = data.frame(id.outcome = "out", outcome = "Outcome")
  )
}

test_that("mv_multiple() nsnp is all SNPs when instrument_specific = FALSE (issue #78)", {
  mvdat <- make_mvdat()
  res <- mv_multiple(mvdat)$result
  expect_equal(res$nsnp, c(40, 40))
  # pval_threshold does not affect the estimates when instrument_specific = FALSE
  res1 <- mv_multiple(mvdat, pval_threshold = 1)$result
  expect_equal(res, res1)
  res_is <- mv_multiple(mvdat, instrument_specific = TRUE)$result
  expect_equal(res_is$nsnp, c(30, 10))
})

test_that("mv_residual() nsnp is all SNPs when instrument_specific = FALSE (issue #78)", {
  mvdat <- make_mvdat()
  res <- mv_residual(mvdat)$result
  expect_equal(res$nsnp, c(40, 40))
  res1 <- mv_residual(mvdat, pval_threshold = 1)$result
  expect_equal(res, res1)
  res_is <- mv_residual(mvdat, instrument_specific = TRUE)$result
  expect_equal(res_is$nsnp, c(30, 10))
})

test_that("mv_residual() estimates are not set to NA when few SNPs pass the threshold and instrument_specific = FALSE (issue #78)", {
  mvdat <- make_mvdat()
  # Only 1 SNP passes the threshold for exp2
  mvdat$exposure_pval[, "exp2"] <- c(1e-10, rep(0.01, 39))
  res <- mv_residual(mvdat)$result
  expect_false(anyNA(res$b))
  res_is <- mv_residual(mvdat, instrument_specific = TRUE)$result
  expect_true(is.na(res_is$b[res_is$id.exposure == "exp2"]))
})
