# Opt-in: these run against a real OpenGWAS Data API, only with credentials (`OPENGWAS_URL` and
# `OPENGWAS_DEV_KEY`, or the Gateway's `OPENGWAS_JWT` once it is live), and never on CRAN.
skip_if_no_opengwasr_credentials <- function() {
  skip_on_cran()
  skip_if_not_installed("opengwasr")
  if (Sys.getenv("OPENGWAS_URL") == "" || Sys.getenv("OPENGWAS_DEV_KEY") == "") {
    skip("No OpenGWAS credentials (OPENGWAS_URL, OPENGWAS_DEV_KEY).")
  }
}

test_that("live: instruments, outcome data, harmonisation and MR through the real API", {
  skip_if_no_opengwasr_credentials()
  src <- opengwasr::og_api()

  exp_dat <- suppressWarnings(extract_instruments("ukb-b-19953", src = src))
  expect_equal(unique(exp_dat$id.exposure), "OGS-00009:ukb-b-19953")
  expect_true(nrow(exp_dat) > 500)
  expect_equal(unique(exp_dat$effect_scale.exposure), "sd")
  expect_equal(unique(exp_dat$clump_ld_choice.exposure), "matched")
  expect_equal(unique(exp_dat$clump_ld.exposure), "hgdp1kgp-v2-EUR-k10")

  exp_dat <- exp_dat[1:50, ]
  out_dat <- extract_outcome_data(exp_dat$SNP, "ebi-a-GCST003116", src = src)
  expect_equal(unique(out_dat$id.outcome), "OGS-00011:GCST003116")
  expect_true(nrow(out_dat) > 40)
  expect_true(all(out_dat$association_status.outcome %in% c("observed", "imputed")))

  h <- suppressMessages(harmonise_data(exp_dat, out_dat))
  res <- suppressMessages(mr(h, method_list = "mr_ivw"))
  expect_true(res$b > 0)
  expect_true(res$nsnp_imputed >= 0)
})

test_that("live: available_outcomes() and add_metadata()", {
  skip_if_no_opengwasr_credentials()
  src <- opengwasr::og_api()
  ao <- available_outcomes(src = src)
  expect_true(all(c("OGS-00009:ukb-b-19953", "OGS-00011:GCST003116") %in% ao$id))

  exp_dat <- suppressWarnings(extract_instruments("ukb-b-19953", src = src))[1:5, ]
  exp_dat$samplesize.exposure <- NULL
  out <- add_metadata(exp_dat, src = src)
  expect_equal(unique(out$samplesize.exposure), 461460)
})

test_that("live: an unknown Analysis is an error naming it", {
  skip_if_no_opengwasr_credentials()
  expect_error(
    extract_instruments("not-a-real-analysis"),
    class = "opengwasr_analysis_not_found"
  )
})
