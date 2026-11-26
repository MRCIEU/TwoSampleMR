# Two exposures, BMI and coronary artery disease, as `mv_extract_exposures()` takes them. The fake
# serves BMI's recorded top hits, and for coronary artery disease the first five Associations
# recorded for it, given the clumping columns an opengwasr top hit has.
local_fake_two_exposures <- function(env = parent.frame()) {
  state <- local_fake_opengwasr(env = env)
  bmi <- fixture("tophits-bmi")
  cad <- fixture("associations-cad")[1:5, ]
  clump_cols <- c("clump_ld", "clump_ld_choice", "clump_r2", "clump_kb", "clump_size")
  cad[, clump_cols] <- bmi[1:5, clump_cols]
  testthat::local_mocked_bindings(
    og_tophits = function(
      src,
      analyses,
      threshold = 5e-8,
      limit = 1000,
      observed_only = FALSE,
      clump = TRUE,
      r2 = NULL,
      kb = NULL,
      original_units = FALSE,
      releases = NULL
    ) {
      state$calls[[length(state$calls) + 1]] <- list(
        fn = "og_tophits",
        args = list(analyses = analyses, r2 = r2, kb = kb)
      )
      rbind(
        if (bmi_id %in% analyses || "ukb-b-19953" %in% analyses) bmi,
        if (cad_id %in% analyses || "ebi-a-GCST003116" %in% analyses) cad
      )
    },
    .package = "opengwasr",
    .env = env
  )
  # The final clumping, through opengwasr, keeps every variant.
  state$ld <- local_fake_ld(keep = "all", env = env)
  state
}

test_that("mv_extract_exposures() extracts instruments and their effects on every exposure through opengwasr", {
  state <- local_fake_two_exposures()
  mv <- mv_extract_exposures(c("ukb-b-19953", "ebi-a-GCST003116"))

  expect_setequal(unique(mv$id.exposure), c(bmi_id, cad_id))
  # Every SNP kept has an effect on each exposure, and the palindromic A/T variants and the ones
  # that fail harmonisation are gone: 11 of the 17 BMI instruments (and CAD's own, all among them).
  expect_equal(nrow(mv), 22)
  expect_true(all(table(mv$SNP) == 2))
  expect_false("rs10203386" %in% mv$SNP)
  expect_true(all(mv$SNP %in% fixture("tophits-bmi")$rsid))

  asked <- calls_to(state, "og_tophits")[[1]]$args
  expect_equal(asked$analyses, c(bmi_id, cad_id))
  expect_null(asked$r2)
  expect_null(asked$kb)
})

test_that("mv_extract_exposures() sends clump_r2 and clump_kb to the LD Store only when they are set", {
  state <- local_fake_two_exposures()
  suppressMessages(mv_extract_exposures(
    c("ukb-b-19953", "ebi-a-GCST003116"),
    clump_r2 = 0.01,
    clump_kb = 3000
  ))
  asked <- calls_to(state, "og_tophits")[[1]]$args
  expect_equal(asked$r2, 0.01)
  expect_equal(asked$kb, 3000)
})

test_that("mv_extract_exposures() says which arguments are retired, only when they are passed", {
  local_fake_two_exposures()
  ids <- c("ukb-b-19953", "ebi-a-GCST003116")
  expect_no_message(mv_extract_exposures(ids), class = "twosamplemr_deprecated_argument")
  expect_message(
    mv_extract_exposures(ids, find_proxies = FALSE),
    "find_proxies.*Imputed Associations",
    class = "twosamplemr_deprecated_argument"
  )
  expect_message(
    mv_extract_exposures(ids, opengwas_jwt = "x"),
    "opengwas_jwt",
    class = "twosamplemr_deprecated_argument"
  )
  expect_message(
    mv_extract_exposures(ids, force_server = TRUE),
    "force_server",
    class = "twosamplemr_deprecated_argument"
  )
})

test_that("mv_extract_exposures() clumps across the exposures through og_clump(), with the Analyses' own ids", {
  state <- local_fake_two_exposures()
  ld <- state$ld
  mv <- mv_extract_exposures(c("ukb-b-19953", "ebi-a-GCST003116"))
  asked <- ld$calls[[1]]
  expect_equal(asked$fn, "og_clump")
  expect_true(all(asked$data$id %in% c(bmi_id, cad_id)))
  expect_false("1" %in% asked$data$id)
  expect_null(asked$r2)
  expect_null(asked$kb)
  expect_length(ld$calls, 1)
  expect_gt(nrow(mv), 0)
})

test_that("mv_extract_exposures() passes its Source to the final clump", {
  state <- local_fake_two_exposures()
  src <- structure(list(), class = "fake_source")
  mv_extract_exposures(c("ukb-b-19953", "ebi-a-GCST003116"), src = src)
  expect_identical(state$ld$calls[[1]]$src, src)
})
