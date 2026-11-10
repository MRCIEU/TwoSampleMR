make_mv_exposure <- function(snps, id, effect_allele = "A", other_allele = "G") {
  data.frame(
    SNP = snps,
    id.exposure = id,
    exposure = paste("Exposure", id),
    effect_allele.exposure = effect_allele,
    other_allele.exposure = other_allele,
    eaf.exposure = 0.3,
    beta.exposure = seq(0.1, by = 0.01, length.out = length(snps)),
    se.exposure = 0.01,
    pval.exposure = 1e-10
  )
}

make_mv_outcome <- function(snps, effect_allele = "A", other_allele = "G") {
  data.frame(
    SNP = snps,
    id.outcome = "out",
    outcome = "Outcome",
    effect_allele.outcome = effect_allele,
    other_allele.outcome = other_allele,
    eaf.outcome = 0.3,
    beta.outcome = 0.05,
    se.outcome = 0.01,
    pval.outcome = 0.01
  )
}

test_that("mv_harmonise_data() keeps only the SNPs present for all exposures", {
  exposure_dat <- rbind(
    make_mv_exposure(c("rs1", "rs2", "rs3", "rs4"), "e1"),
    make_mv_exposure(c("rs1", "rs2", "rs3", "rs5"), "e2")
  )
  outcome_dat <- make_mv_outcome(paste0("rs", 1:5))
  mvdat <- suppressMessages(mv_harmonise_data(exposure_dat, outcome_dat))
  expect_equal(rownames(mvdat$exposure_beta), c("rs1", "rs2", "rs3"))
  expect_equal(colnames(mvdat$exposure_beta), c("e1", "e2"))
  expect_length(mvdat$outcome_beta, 3)
})

test_that("mv_harmonise_data() errors when no SNP is present for all exposures (issue #316)", {
  # e.g. exposure data from extract_instruments() for each exposure, rather than
  # mv_extract_exposures()
  exposure_dat <- rbind(
    make_mv_exposure(c("rs1", "rs2", "rs3"), "e1"),
    make_mv_exposure(c("rs4", "rs5", "rs6"), "e2")
  )
  outcome_dat <- make_mv_outcome(paste0("rs", 1:6))
  expect_error(
    suppressMessages(mv_harmonise_data(exposure_dat, outcome_dat)),
    "None of the SNPs in exposure_dat are present for all 2 exposures",
    fixed = TRUE
  )
  expect_error(
    suppressMessages(mv_harmonise_data(exposure_dat, outcome_dat)),
    "mv_extract_exposures()",
    fixed = TRUE
  )

  # SNPs shared by some, but not all, of the exposures
  exposure_dat <- rbind(
    make_mv_exposure(c("rs1", "rs2"), "e1"),
    make_mv_exposure(c("rs1", "rs2"), "e2"),
    make_mv_exposure(c("rs3", "rs4"), "e3")
  )
  expect_error(
    suppressMessages(mv_harmonise_data(exposure_dat, outcome_dat)),
    "None of the SNPs in exposure_dat are present for all 3 exposures",
    fixed = TRUE
  )
})

test_that("mv_harmonise_data() errors when no shared SNPs are in outcome_dat (issue #316)", {
  exposure_dat <- rbind(
    make_mv_exposure(c("rs1", "rs2", "rs3"), "e1"),
    make_mv_exposure(c("rs1", "rs2", "rs3"), "e2")
  )
  outcome_dat <- make_mv_outcome(c("rs10", "rs11", "rs12"))
  expect_error(
    suppressMessages(mv_harmonise_data(exposure_dat, outcome_dat)),
    "None of the 3 SNPs present for all exposures were found in outcome_dat",
    fixed = TRUE
  )
})

test_that("mv_harmonise_data() errors when no SNPs remain after harmonising (issue #316)", {
  # Palindromic SNPs are removed with harmonise_strictness = 3
  exposure_dat <- rbind(
    make_mv_exposure(c("rs1", "rs2", "rs3"), "e1", "A", "T"),
    make_mv_exposure(c("rs1", "rs2", "rs3"), "e2", "A", "T")
  )
  outcome_dat <- make_mv_outcome(c("rs1", "rs2", "rs3"), "A", "T")
  expect_error(
    suppressMessages(mv_harmonise_data(exposure_dat, outcome_dat, harmonise_strictness = 3)),
    "None of the 3 SNPs present for all exposures and in outcome_dat were kept after harmonising",
    fixed = TRUE
  )

  # Incompatible alleles
  exposure_dat <- rbind(
    make_mv_exposure(c("rs1", "rs2", "rs3"), "e1"),
    make_mv_exposure(c("rs1", "rs2", "rs3"), "e2")
  )
  outcome_dat <- make_mv_outcome(c("rs1", "rs2", "rs3"), "A", "C")
  expect_error(
    suppressMessages(mv_harmonise_data(exposure_dat, outcome_dat)),
    "harmonise_strictness",
    fixed = TRUE
  )
})

test_that("mv_harmonise_data() handles duplicated SNPs within an exposure", {
  # rs1 is duplicated in e1 and missing from e3, so a naive count of rows per
  # SNP (3) would wrongly treat it as present for all 3 exposures
  exposure_dat <- rbind(
    make_mv_exposure(c("rs1", "rs1", "rs2"), "e1"),
    make_mv_exposure(c("rs1", "rs2"), "e2"),
    make_mv_exposure("rs2", "e3")
  )
  outcome_dat <- make_mv_outcome(c("rs1", "rs2"))
  expect_warning(
    mvdat <- suppressMessages(mv_harmonise_data(exposure_dat, outcome_dat)),
    "1 duplicated SNP-exposure rows found in exposure_dat",
    fixed = TRUE
  )
  expect_equal(rownames(mvdat$exposure_beta), "rs2")
  expect_equal(colnames(mvdat$exposure_beta), c("e1", "e2", "e3"))
  expect_false(is.list(mvdat$exposure_beta))
})
