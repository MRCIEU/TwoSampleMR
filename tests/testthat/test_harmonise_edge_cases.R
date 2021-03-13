context("harmonise")

set.seed(1)

# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)


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

df_exp <- format_data(
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
)

df_out <- format_data(
    df,
    type = "outcome",
    snp_col = "SNP",
    pval_col = "pval",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    phenotype_col = "pheno_id",
    samplesize_col = "n"
)


test_that("harmonise_data works when exposure and outcome df are 1 row.", {
    for (i in seq(1,3)) {
        result <- harmonise_data(
            exposure_dat = df_exp[1,],
            outcome_dat = df_out[1,],
            action = 1
        )
        expect_equal(nrow(result), 1)
    }
})


test_that("harmonise_data works when there are no matching SNPs.", {
    for (i in seq(1,3)) {
        df_out$SNP <- paste(df_out$SNP, "foo", sep = "_")
        result <- harmonise_data(
            exposure_dat = df_exp,
            outcome_dat = df_out,
            action = 1
        )
        expect_equal(nrow(result), 0)
    }
})


