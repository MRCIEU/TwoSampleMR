test_that("available_outcomes() lists the Analyses opengwasr holds", {
  state <- local_fake_opengwasr()
  ao <- available_outcomes()
  expect_equal(ao$id, c(bmi_id, cad_id))
  expect_equal(ao$trait, c("Body mass index (BMI)", "Coronary artery disease"))
  expect_equal(ao$analysis_label, ao$trait)
  expect_length(calls_to(state, "og_analyses"), 1)
  expect_null(calls_to(state, "og_analyses")[[1]]$args$ids)
})

test_that("available_outcomes() passes the Source, and says opengwas_jwt is retired", {
  state <- local_fake_opengwasr()
  my_src <- structure(list(), class = "my_source")
  available_outcomes(src = my_src)
  expect_identical(calls_to(state, "og_analyses")[[1]]$args$src, my_src)
  expect_message(
    available_outcomes(opengwas_jwt = "x"),
    "opengwas_jwt",
    class = "twosamplemr_deprecated_argument"
  )
})

test_that("add_metadata() adds sample sizes, case counts and the effect scale to exposure data", {
  local_fake_opengwasr()
  exp_dat <- extract_instruments("ukb-b-19953")
  exp_dat$samplesize.exposure <- NULL
  exp_dat$effect_scale.exposure <- NULL
  out <- add_metadata(exp_dat)
  expect_equal(unique(out$samplesize.exposure), 461460)
  expect_equal(unique(out$effect_scale.exposure), "sd")
  expect_true(all(is.na(out$ncase.exposure)))
  expect_true(all(is.na(out$sd.exposure)))
  expect_equal(unique(out$units.exposure), "SD")
  expect_equal(out$SNP, exp_dat$SNP)
})

test_that("add_metadata() reports case and control counts for outcome data", {
  local_fake_opengwasr()
  dat <- extract_outcome_data("rs6265", "ebi-a-GCST003116")
  dat$samplesize.outcome <- NULL
  dat$ncase.outcome <- NULL
  dat$ncontrol.outcome <- NULL
  out <- add_metadata(dat)
  expect_equal(out$samplesize.outcome, 187599)
  expect_equal(out$ncase.outcome, 61289)
  expect_equal(out$ncontrol.outcome, 126310)
  expect_equal(out$effect_scale.outcome, "log_or")
})

test_that("add_metadata() accepts the retired `unit` field and reports effect_scale instead", {
  local_fake_opengwasr()
  exp_dat <- extract_instruments("ukb-b-19953")
  exp_dat$effect_scale.exposure <- NULL
  expect_message(
    out <- add_metadata(exp_dat, cols = c("sample_size", "unit")),
    "effect_scale",
    class = "twosamplemr_deprecated_argument"
  )
  expect_equal(unique(out$effect_scale.exposure), "sd")
  expect_equal(unique(out$units.exposure), "SD")
})

test_that("add_metadata() copes with ids that name no Analysis", {
  local_fake_opengwasr()
  exp_dat <- extract_instruments("ukb-b-19953")
  exp_dat$id.exposure <- "not a real id"
  exp_dat$effect_scale.exposure <- NULL
  exp_dat$samplesize.exposure <- NULL
  expect_message(out <- add_metadata(exp_dat), "none of the IDs found")
  expect_false("effect_scale.exposure" %in% names(out))
  expect_false("samplesize.exposure" %in% names(out))

  mixed <- rbind(exp_dat[1:2, ], transform(exp_dat[3:4, ], id.exposure = bmi_id))
  out <- add_metadata(mixed)
  expect_equal(out$samplesize.exposure, c(NA, NA, 461460, 461460))
})

test_that("add_metadata() works on legacy ids, which opengwasr resolves", {
  local_fake_opengwasr()
  exp_dat <- extract_instruments("ukb-b-19953")
  exp_dat$id.exposure <- "ukb-b-19953"
  exp_dat$samplesize.exposure <- NULL
  out <- add_metadata(exp_dat)
  expect_equal(unique(out$samplesize.exposure), 461460)
})
