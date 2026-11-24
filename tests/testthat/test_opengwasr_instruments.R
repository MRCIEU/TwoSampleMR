test_that("extract_instruments() returns exposure data in TwoSampleMR's format", {
  local_fake_opengwasr()
  exp_dat <- extract_instruments("ukb-b-19953")

  expect_s3_class(exp_dat, "data.frame")
  expect_equal(nrow(exp_dat), 17)
  expect_true(all(
    c(
      "SNP",
      "beta.exposure",
      "se.exposure",
      "pval.exposure",
      "eaf.exposure",
      "effect_allele.exposure",
      "other_allele.exposure",
      "chr.exposure",
      "pos.exposure",
      "samplesize.exposure",
      "exposure",
      "id.exposure",
      "mr_keep.exposure"
    ) %in%
      names(exp_dat)
  ))
  expect_equal(unique(exp_dat$id.exposure), bmi_id)
  expect_equal(unique(exp_dat$exposure), paste0("Body mass index (BMI) || id:", bmi_id))
  expect_equal(unique(exp_dat$samplesize.exposure), 461460)
  expect_equal(exp_dat$SNP[1:3], c("rs56094641", "rs6567160", "rs6744646"))
  expect_equal(exp_dat$beta.exposure[1], -0.07357775, tolerance = 1e-6)
  expect_equal(exp_dat$effect_allele.exposure[1:2], c("A", "C"))
})

test_that("extract_instruments() uses the ALID as SNP where a variant has no rsID", {
  local_fake_opengwasr()
  exp_dat <- extract_instruments("ukb-b-19953")
  expect_equal(
    exp_dat$SNP[13:16],
    c("2:25087361:G:GTA", "5:75701728:A:AAAAT", "7:75543861:A:G", "16:53788264:T:TC")
  )
})

test_that("extract_instruments() reports the effect scale and the clumping that was applied", {
  local_fake_opengwasr()
  exp_dat <- extract_instruments("ukb-b-19953")
  expect_equal(unique(exp_dat$effect_scale.exposure), "sd")
  expect_equal(sum(exp_dat$association_status.exposure == "imputed"), 4)
  expect_equal(sum(exp_dat$association_status.exposure == "observed"), 13)
  expect_equal(unique(exp_dat$clump_ld.exposure), "hgdp1kgp-v2-EUR-k10")
  expect_equal(unique(exp_dat$clump_ld_choice.exposure), "matched")
  expect_equal(unique(exp_dat$clump_r2.exposure), 0.03186544, tolerance = 1e-6)
  expect_equal(unique(exp_dat$clump_kb.exposure), 10000)
})

test_that("extract_instruments() leaves clumping thresholds to the LD Store unless the user sets them", {
  state <- local_fake_opengwasr()
  extract_instruments("ukb-b-19953")
  sent <- calls_to(state, "og_tophits")[[1]]$args
  expect_null(sent$r2)
  expect_null(sent$kb)
  expect_true(sent$clump)
  expect_equal(sent$threshold, 5e-8)

  extract_instruments("ukb-b-19953", p1 = 5e-6, r2 = 0.01, kb = 5000)
  sent <- calls_to(state, "og_tophits")[[2]]$args
  expect_equal(sent$r2, 0.01)
  expect_equal(sent$kb, 5000)
  expect_equal(sent$threshold, 5e-6)
})

test_that("extract_instruments(clump = FALSE) asks for every top hit and carries no clumping columns", {
  state <- local_fake_opengwasr()
  exp_dat <- extract_instruments("ukb-b-19953", clump = FALSE)
  expect_false(calls_to(state, "og_tophits")[[1]]$args$clump)
  expect_false(any(grepl("^clump_", names(exp_dat))))
  expect_error(extract_instruments("ukb-b-19953", clump = 2), "clump argument should be 0 or 1")
})

test_that("extract_instruments() passes the Source it is given to opengwasr", {
  state <- local_fake_opengwasr()
  my_src <- structure(list(), class = "my_source")
  extract_instruments("ukb-b-19953", src = my_src)
  expect_identical(calls_to(state, "og_tophits")[[1]]$args$src, my_src)
})

test_that("extract_instruments() returns NULL when no Analysis has instruments", {
  local_fake_opengwasr(empty = cad_id)
  expect_null(extract_instruments("ebi-a-GCST003116"))
})

test_that("extract_instruments() says that opengwas_jwt and force_server are retired, only when passed", {
  local_fake_opengwasr()
  expect_no_message(extract_instruments("ukb-b-19953"), class = "twosamplemr_deprecated_argument")
  expect_message(
    extract_instruments("ukb-b-19953", opengwas_jwt = "x"),
    "opengwas_jwt.*og_api\\(\\)",
    class = "twosamplemr_deprecated_argument"
  )
  expect_message(
    extract_instruments("ukb-b-19953", force_server = TRUE),
    "force_server",
    class = "twosamplemr_deprecated_argument"
  )
})
