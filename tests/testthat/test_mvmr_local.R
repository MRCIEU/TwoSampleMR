context("mvmr local")

# a <- ieugwasr::tophits("ieu-a-2")
# b <- ieugwasr::tophits("ieu-a-1001")
# rsid <- unique(c(a$rsid, b$rsid))
# a1 <- ieugwasr::associations(rsid, "ieu-a-2")
# a2 <- ieugwasr::associations(rsid, "ieu-a-1001")
# save(a, b, rsid, a1, a2, file="inst/extdata/test_add_mvmr_local.RData", compress = "xz")

load(system.file("extdata", "test_add_mvmr_local.RData", package="TwoSampleMR"))

test_that("mv exposure local", {
  skip_if_offline()
  skip_if_offline(host = "api.opengwas.io")
	skip_on_cran()
    f1 <- tempfile()
    f2 <- tempfile()
    write.table(a1, file=f1, row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
    write.table(a2, file=f2, row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")

    exposure_dat <- try(mv_extract_exposures_local(
        c(f1, f2),
        sep = "\t",
        snp_col=c("rsid"),
        beta_col=c("beta"),
        se_col=c("se"),
        effect_allele_col=c("ea"),
        other_allele_col=c("nea"),
        pval_col=c("p")
    ))
    if (inherits(exposure_dat, "try-error")) skip("Server issues")
    expect_true(nrow(exposure_dat) > 100)

    exposure_dat2 <- try(mv_extract_exposures_local(
        list(a1, a2),
        sep = "\t",
        snp_col=c("rsid"),
        beta_col=c("beta"),
        se_col=c("se"),
        effect_allele_col=c("ea"),
        other_allele_col=c("nea"),
        pval_col=c("p")
    ))
    if (inherits(exposure_dat2, "try-error")) skip("Server issues")
    expect_true(all.equal(exposure_dat, exposure_dat2))
})
