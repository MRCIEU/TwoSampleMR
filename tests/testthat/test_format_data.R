df <- data.frame(
  "SNP" = c("9_69001927_C_T", "9_69459263_A_G", "9_69508544_G_A"),
  "effect_allele" = c("T", "G", "A"),
  "other_allele" = c("C", "A", "G"),
  "se" = c(0.01664, 0.02038, 0.04585),
  "beta" = c(0.367464, -0.265656, 0.254032),
  "eaf" = c(0.319894, 0.39234, 0.343751),
  "pval" = c(0.250677, 0.498338, 0.459907),
  "n" = c(30100, 30100, 30100),
  "pheno_id" = c("traita", "traita", "traita")
)

df <- data.table::data.table(df)

test_that("format_data() should error on a dat object of class data.table", {
  expect_error(format_data(
    df,
    type = "exposure",
    snp_col = "SNP",
    pval_col = "pval",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    phenotype_col = "pheno_id",
    samplesize_col = "n"
  ))
})

df <- data.frame(df)

test_that("format_data() should not error after having its data.table class removed", {
  expect_no_error(format_data(
    df,
    type = "exposure",
    snp_col = "SNP",
    pval_col = "pval",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    phenotype_col = "pheno_id",
    samplesize_col = "n"
  ))
})

test_that("format_data() should not cause a stack overflow", {
  a <- data.table::data.table(x = sample(1:1e6))
  expect_error(do.call(format_data, list(a)))
})

test_that("Inferred p-value calculation is correct for outcome", {
  df2 <- data.frame(
    "SNP" = c("9_69001927_C_T", "9_69459263_A_G"),
    "effect_allele" = c("T", "G"),
    "other_allele" = c("C", "A"),
    "se" = c(0.6, 0.3),
    "beta" = c(1.1, 0.5),
    "eaf" = c(0.319894, 0.39234),
    "pval" = c(NA, NA),
    "pheno_id" = c("traita", "traita")
  )

  suppressWarnings({
    datformat <- format_data(df2, type = "outcome")
  })
  pvals <- 2 * stats::pnorm(abs(datformat$beta.outcome) / datformat$se.outcome, lower.tail = FALSE)
  expect_equal(datformat$pval.outcome, pvals)
  expect_equal(datformat$pval_origin.outcome, c("inferred", "inferred"))
})

test_that("Inferred p-value calculation is correct for exposure", {
  df3 <- data.frame(
    "SNP" = c("9_69001927_C_T", "9_69459263_A_G"),
    "effect_allele" = c("T", "G"),
    "other_allele" = c("C", "A"),
    "se" = c(0.6, 0.3),
    "beta" = c(1.1, 0.5),
    "eaf" = c(0.319894, 0.39234),
    "pval" = c(NA, NA),
    "pheno_id" = c("traita", "traita")
  )

  suppressWarnings({
    datformat2 <- format_data(df3, type = "exposure")
  })
  pvals2 <- 2 *
    stats::pnorm(abs(datformat2$beta.exposure) / datformat2$se.exposure, lower.tail = FALSE)
  expect_equal(datformat2$pval.exposure, pvals2)
  expect_equal(datformat2$pval_origin.exposure, c("inferred", "inferred"))
})
