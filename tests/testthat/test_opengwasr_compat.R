harmonised <- function() {
  local_fake_opengwasr(env = parent.frame())
  exp_dat <- extract_instruments("ukb-b-19953")
  out_dat <- suppressMessages(extract_outcome_data(exp_dat$SNP, "ebi-a-GCST003116"))
  suppressMessages(harmonise_data(exp_dat, out_dat))
}

# Units: the readers of `units.exposure` and `units.outcome` still work -----------------------------

test_that("add_metadata() adds units.* from effect_scale, as it always has", {
  dat <- harmonised()
  dat$units.exposure <- NULL
  dat$units.outcome <- NULL
  out <- add_metadata(dat)
  expect_equal(unique(out$units.exposure), "SD")
  expect_equal(unique(out$units.outcome), "log odds")
  expect_equal(nrow(out), nrow(dat))
})

test_that("an effect scale with no units label is NA, with a message", {
  expect_message(units <- units_from_effect_scale(c("sd", "log_or", "original")), "original")
  expect_equal(units, c("SD", "log odds", NA))
})

test_that("standardise_units() runs on data from add_metadata() and leaves log odds alone", {
  dat <- add_metadata(harmonised())
  out <- standardise_units(dat)
  expect_equal(nrow(out), nrow(dat))
  expect_equal(unique(out$units.outcome), "log odds")
  expect_equal(out$beta.outcome, dat$beta.outcome)
})

test_that("add_rsq() takes the binary path for a log odds outcome", {
  dat <- add_metadata(harmonised())
  out <- suppressWarnings(add_rsq(dat))
  expected <- get_r_from_lor(
    out$beta.outcome,
    out$eaf.outcome,
    out$ncase.outcome,
    out$ncontrol.outcome,
    0.1
  )^2
  expect_equal(out$rsq.outcome, expected)
  expect_equal(
    out$rsq.exposure,
    2 * out$beta.exposure^2 * out$eaf.exposure * (1 - out$eaf.exposure)
  )
})

test_that("steiger_filtering() and directionality_test() run on data from add_metadata()", {
  dat <- add_metadata(harmonised())
  st <- suppressWarnings(steiger_filtering(dat))
  expect_true(all(c("steiger_dir", "steiger_pval") %in% names(st)))
  dir <- suppressWarnings(directionality_test(dat))
  expect_equal(nrow(dir), 1)
})

test_that("get_rsq() and run_mrmix() read units.* from add_metadata()", {
  dat <- add_metadata(harmonised())
  dat <- dat[dat$mr_keep, ]
  r <- suppressWarnings(get_rsq(dat))
  expect_equal(
    r$rsq.outcome,
    get_r_from_lor(r$beta.outcome, r$eaf.outcome, r$ncase.outcome, r$ncontrol.outcome, 0.1)^2
  )
  expect_no_error(suppressWarnings(suppressMessages(run_mrmix(dat))))
})

# ao_slc = TRUE keeps working with the new catalogue -----------------------------------------------

test_that("format_mr_results() works at its default ao_slc = TRUE, and says what it skips", {
  dat <- harmonised()
  res <- suppressMessages(mr(dat))
  expect_message(out <- format_mr_results(res), "subcategory")
  expect_s3_class(out, "data.frame")
  expect_equal(unique(out$outcome), "Coronary artery disease")
  expect_equal(unique(out$sample_size), 187599)
  expect_equal(out$effect, res$b[res$method == "Inverse variance weighted"])
})

test_that("combine_all_mrresults() and the forest plots work at ao_slc = TRUE", {
  dat <- harmonised()
  res <- suppressMessages(mr(dat))
  het <- suppressMessages(mr_heterogeneity(dat))
  plt <- suppressMessages(mr_pleiotropy_test(dat))
  sin <- suppressMessages(mr_singlesnp(dat))
  all <- suppressMessages(combine_all_mrresults(res, het, plt, sin))
  expect_s3_class(all, "data.frame")
  expect_equal(unique(all$id.outcome), cad_id)

  p <- suppressWarnings(suppressMessages(forest_plot_1_to_many(res, ao_slc = TRUE)))
  expect_s3_class(p, "ggplot")
})

# Numeric legacy ids --------------------------------------------------------------------------------

test_that("numeric ids are mapped as the old OpenGWAS client mapped them, before opengwasr sees them", {
  state <- local_fake_opengwasr()
  suppressWarnings(try(extract_instruments(c(2, 300, 303, "UKB-b:19953")), silent = TRUE))
  asked <- calls_to(state, "og_analyses")[[1]]$args$ids
  expect_equal(asked, c("ieu-a-2", "ieu-a-300", "met-a-303", "ukb-b-19953"))
})

# Several ids, some unknown or without instruments --------------------------------------------------

test_that("extract_instruments() copes with some Analyses having no instruments", {
  local_fake_opengwasr(empty = cad_id)
  exp_dat <- extract_instruments(c("ukb-b-19953", "ebi-a-GCST003116"))
  expect_equal(unique(exp_dat$id.exposure), bmi_id)

  local_fake_opengwasr()
  both <- extract_instruments(c("ukb-b-19953", "ebi-a-GCST003116"))
  expect_equal(unique(both$id.exposure), c(bmi_id, cad_id))
  expect_equal(nrow(both), 17 + 10)
})

test_that("an id that names no Analysis is a warning when others are known, an error when none is", {
  local_fake_opengwasr()
  expect_warning(
    exp_dat <- extract_instruments(c("ukb-b-19953", "not-an-id")),
    "not-an-id",
    class = "twosamplemr_unknown_analyses"
  )
  expect_equal(unique(exp_dat$id.exposure), bmi_id)

  expect_warning(
    out <- suppressMessages(extract_outcome_data("rs6265", c("ebi-a-GCST003116", "a"))),
    "a",
    class = "twosamplemr_unknown_analyses"
  )
  expect_equal(unique(out$id.outcome), cad_id)

  expect_error(extract_instruments("not-an-id"), class = "opengwasr_analysis_not_found")
  expect_error(extract_outcome_data("rs6265", "a"), class = "opengwasr_analysis_not_found")
})

test_that("extract_outcome_data() takes several outcomes", {
  local_fake_opengwasr()
  snps <- c("rs539515", "rs6265")
  out <- extract_outcome_data(snps, c("ebi-a-GCST003116", "ukb-b-19953"))
  expect_equal(sort(unique(out$id.outcome)), sort(c(bmi_id, cad_id)))
  expect_equal(nrow(out), 4)
  expect_equal(out$ncase.outcome[out$id.outcome == bmi_id], c(NA_real_, NA_real_))
})

test_that("extract_outcome_data() says that splitsize is no longer needed", {
  local_fake_opengwasr()
  expect_message(
    extract_outcome_data("rs6265", "ebi-a-GCST003116", splitsize = 100),
    "splitsize",
    class = "twosamplemr_deprecated_argument"
  )
})

# A binary exposure, and the round trip through a file ----------------------------------------------

test_that("extract_instruments() maps case and control counts for a binary exposure", {
  local_fake_opengwasr()
  exp_dat <- extract_instruments("ebi-a-GCST003116")
  expect_equal(nrow(exp_dat), 10)
  expect_equal(unique(exp_dat$ncase.exposure), 61289)
  expect_equal(unique(exp_dat$ncontrol.exposure), 126310)
  expect_equal(unique(exp_dat$samplesize.exposure), 187599)
  expect_equal(unique(exp_dat$effect_scale.exposure), "log_or")
})

test_that("instrument data written to a file read back as exposure and outcome data", {
  local_fake_opengwasr()
  exp_dat <- extract_instruments("ukb-b-19953")
  names(exp_dat) <- gsub(".exposure", "", names(exp_dat))
  fn <- tempfile()
  on.exit(unlink(fn))
  utils::write.table(
    exp_dat,
    file = fn,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
  )

  a <- read_exposure_data(fn, sep = "\t")
  b <- read_outcome_data(fn, sep = "\t")
  expect_true("chr.exposure" %in% names(a))
  expect_true("chr.outcome" %in% names(b))
  expect_equal(nrow(a), 17)
  expect_no_warning(read_exposure_data(fn, sep = "\t"))
})

test_that("add_metadata() works on harmonised data and keeps its rows in order", {
  dat <- harmonised()
  dat$samplesize.exposure <- NULL
  dat$samplesize.outcome <- NULL
  out <- add_metadata(dat)
  expect_equal(out$SNP, dat$SNP)
  expect_equal(unique(out$samplesize.exposure), 461460)
  expect_equal(unique(out$samplesize.outcome), 187599)
})

# Clumping provenance in mr() ----------------------------------------------------------------------

test_that("mr() carries the clumping provenance of each exposure, only when the data have it", {
  dat <- harmonised()
  res <- suppressMessages(mr(dat, method_list = "mr_ivw"))
  expect_equal(res$clump_ld.exposure, "hgdp1kgp-v2-EUR-k10")
  expect_equal(res$clump_ld_choice.exposure, "matched")
  expect_equal(res$clump_r2.exposure, 0.03186544, tolerance = 1e-6)
  expect_equal(res$clump_kb.exposure, 10000)
  expect_equal(names(res)[10], "nsnp_imputed")

  dat <- dat[, !grepl("^clump_", names(dat))]
  res <- suppressMessages(mr(dat, method_list = "mr_ivw"))
  expect_false(any(grepl("^clump_", names(res))))
})
