bmi_snps <- function() {
  e <- fixture("tophits-bmi")
  ifelse(is.na(e$rsid), e$alid, e$rsid)
}

test_that("extract_outcome_data() returns outcome data in TwoSampleMR's format", {
  local_fake_opengwasr()
  out_dat <- extract_outcome_data(bmi_snps(), "ebi-a-GCST003116")

  expect_equal(nrow(out_dat), 17)
  expect_true(all(
    c(
      "SNP",
      "chr",
      "pos",
      "beta.outcome",
      "se.outcome",
      "samplesize.outcome",
      "pval.outcome",
      "eaf.outcome",
      "effect_allele.outcome",
      "other_allele.outcome",
      "outcome",
      "id.outcome",
      "mr_keep.outcome",
      "ncase.outcome",
      "ncontrol.outcome"
    ) %in%
      names(out_dat)
  ))
  expect_equal(unique(out_dat$id.outcome), cad_id)
  expect_equal(unique(out_dat$outcome), paste0("Coronary artery disease || id:", cad_id))
  expect_equal(unique(out_dat$samplesize.outcome), 187599)
  expect_equal(unique(out_dat$ncase.outcome), 61289)
  expect_equal(unique(out_dat$ncontrol.outcome), 126310)
  expect_true("rs539515" %in% out_dat$SNP)
  expect_true("2:25087361:G:GTA" %in% out_dat$SNP)
  expect_true(all(out_dat$mr_keep.outcome))
})

test_that("extract_outcome_data() labels the effect scale and which Associations were Imputed", {
  local_fake_opengwasr()
  out_dat <- extract_outcome_data(bmi_snps(), "ebi-a-GCST003116")
  expect_equal(unique(out_dat$effect_scale.outcome), "log_or")
  expect_equal(sum(out_dat$association_status.outcome == "imputed"), 4)
  expect_equal(sum(out_dat$association_status.outcome == "observed"), 13)
})

test_that("extract_outcome_data() asks for each SNP once, in the Analyses it names", {
  state <- local_fake_opengwasr()
  extract_outcome_data(
    c("rs539515", "rs539515", "rs6265"),
    c("ebi-a-GCST003116", "ebi-a-GCST003116")
  )
  sent <- calls_to(state, "og_associations")[[1]]$args
  expect_equal(sent$variants, c("rs539515", "rs6265"))
  expect_equal(sent$analyses, cad_id)
})

test_that("extract_outcome_data() returns NULL when nothing is found", {
  local_fake_opengwasr()
  expect_null(extract_outcome_data("rs0000", "ebi-a-GCST003116"))
})

test_that("extract_outcome_data() passes the Source it is given to opengwasr", {
  state <- local_fake_opengwasr()
  my_src <- structure(list(), class = "my_source")
  extract_outcome_data("rs6265", "ebi-a-GCST003116", src = my_src)
  expect_identical(calls_to(state, "og_associations")[[1]]$args$src, my_src)
})

test_that("extract_outcome_data() explains that proxies are replaced, only when they are passed", {
  local_fake_opengwasr()
  expect_no_message(
    extract_outcome_data("rs6265", "ebi-a-GCST003116"),
    class = "twosamplemr_deprecated_argument"
  )
  for (args in list(
    list(proxies = FALSE),
    list(rsq = 0.9),
    list(align_alleles = 0),
    list(palindromes = 0),
    list(maf_threshold = 0.4),
    list(proxy_splitsize = 100)
  )) {
    expect_message(
      do.call(extract_outcome_data, c(list("rs6265", "ebi-a-GCST003116"), args)),
      "Imputed Associations",
      class = "twosamplemr_deprecated_argument"
    )
  }
})

test_that("extract_outcome_data() says that opengwas_jwt is retired, only when passed", {
  local_fake_opengwasr()
  expect_message(
    extract_outcome_data("rs6265", "ebi-a-GCST003116", opengwas_jwt = "x"),
    "opengwas_jwt",
    class = "twosamplemr_deprecated_argument"
  )
})

test_that("make_dat() does not trigger the proxies message unless asked to", {
  local_fake_opengwasr()
  expect_no_message(
    suppressMessages(make_dat("ukb-b-19953", "ebi-a-GCST003116")),
    class = "twosamplemr_deprecated_argument"
  )
})
