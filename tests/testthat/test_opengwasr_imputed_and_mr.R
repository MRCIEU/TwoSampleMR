# BMI (per SD) instruments and their coronary artery disease (log odds) Associations, 4 of them
# Imputed, put through harmonise_data() and mr().
bmi_cad <- function() {
  local_fake_opengwasr(env = parent.frame())
  exp_dat <- extract_instruments("ukb-b-19953")
  out_dat <- suppressMessages(extract_outcome_data(exp_dat$SNP, "ebi-a-GCST003116"))
  list(exp = exp_dat, out = out_dat)
}

test_that("harmonise_data() carries each outcome Association's status and says how many were Imputed", {
  d <- bmi_cad()
  expect_message(
    h <- harmonise_data(d$exp, d$out),
    "[0-9]+ of [0-9]+ outcome Associations used for MR \\(Coronary artery disease.*\\) were Imputed"
  )
  expect_true("association_status.outcome" %in% names(h))
  expect_equal(sum(h$association_status.outcome == "imputed"), 4)
  expect_equal(sum(h$mr_keep & h$association_status.outcome == "imputed"), 4)
  expect_message(
    harmonise_data(d$exp, d$out),
    "4 of 15 outcome Associations used for MR"
  )
})

test_that("harmonise_data() says nothing about Imputed Associations for data that have no status", {
  d <- bmi_cad()
  d$out$association_status.outcome <- NULL
  msgs <- testthat::capture_messages(h <- harmonise_data(d$exp, d$out))
  expect_false(any(grepl("Imputed", msgs)))
  expect_false("association_status.outcome" %in% names(h))
})

test_that("mr() carries the number of Imputed outcome Associations per exposure-outcome pair", {
  d <- bmi_cad()
  h <- suppressMessages(harmonise_data(d$exp, d$out))
  res <- suppressMessages(mr(h, method_list = c("mr_ivw", "mr_egger_regression")))
  expect_equal(names(res)[1:9], c(
    "id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval"
  ))
  expect_equal(names(res)[10], "nsnp_imputed")
  expect_equal(res$nsnp_imputed, c(4, 4))

  # Without the status the output is what it always was.
  h$association_status.outcome <- NULL
  res <- suppressMessages(mr(h, method_list = "mr_ivw"))
  expect_false("nsnp_imputed" %in% names(res))
})

test_that("mr() on opengwasr data gives the per-SD estimate worked out by hand", {
  d <- bmi_cad()
  h <- suppressMessages(harmonise_data(d$exp, d$out))
  expect_equal(unique(h$effect_scale.exposure), "sd")
  expect_equal(unique(h$effect_scale.outcome), "log_or")

  # One instrument: the Wald ratio of its CAD log odds ratio to its BMI effect (SD per allele).
  # rs56094641: BMI beta -0.07357775 (se 0.002016341), CAD beta -0.029787272 (se 0.009613037).
  wald <- suppressMessages(mr(h[h$SNP == "rs56094641", ], method_list = "mr_wald_ratio"))
  expect_equal(wald$b, -0.029787272 / -0.07357775, tolerance = 1e-6)
  expect_equal(wald$se, 0.009613037 / 0.07357775, tolerance = 1e-6)

  # All 15 that survive harmonisation, inverse variance weighted, from the formula
  # b = sum(bx * by / se_y^2) / sum(bx^2 / se_y^2), with the standard error of the fixed effect
  # model times the residual standard error where that exceeds 1 (multiplicative random effects).
  kept <- h[h$mr_keep, ]
  bx <- kept$beta.exposure
  by <- kept$beta.outcome
  sy <- kept$se.outcome
  ivw <- suppressMessages(mr(h, method_list = "mr_ivw"))
  b <- sum(bx * by / sy^2) / sum(bx^2 / sy^2)
  residual_se <- sqrt(sum((by - b * bx)^2 / sy^2) / (length(bx) - 1))
  expect_equal(ivw$nsnp, 15)
  expect_equal(ivw$b, b, tolerance = 1e-9)
  expect_equal(ivw$se, 1 / sqrt(sum(bx^2 / sy^2)) * max(1, residual_se), tolerance = 1e-9)
  # The same numbers, as mr() gave them before the switch to opengwasr on these data.
  expect_equal(ivw$b, 0.5359352949, tolerance = 1e-8)
  expect_equal(ivw$se, 0.08094567703, tolerance = 1e-8)
})
