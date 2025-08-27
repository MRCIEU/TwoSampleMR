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
