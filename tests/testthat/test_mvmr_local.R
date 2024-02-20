context("mvmr local")

test_that("mv exposure local", {
	skip_on_ci()
	skip_on_cran()

    a <- ieugwasr::tophits("ieu-a-2")
    b <- ieugwasr::tophits("ieu-a-1001")
    rsid <- unique(c(a$rsid, b$rsid))
    a1 <- ieugwasr::associations(rsid, "ieu-a-2")
    a2 <- ieugwasr::associations(rsid, "ieu-a-1001")

    f1 <- tempfile()
    f2 <- tempfile()
    write.table(a1, file=f1, row=F, col=T, qu=F, sep="\t")
    write.table(a2, file=f2, row=F, col=T, qu=F, sep="\t")

    exposure_dat <- mv_extract_exposures_local(
        c(f1, f2),
        sep = "\t",
        snp_col=c("rsid"),
        beta_col=c("beta"),
        se_col=c("se"),
        effect_allele_col=c("ea"),
        other_allele_col=c("nea"),
        pval_col=c("p")
    )

    expect_true(nrow(exposure_dat) > 100)
})

