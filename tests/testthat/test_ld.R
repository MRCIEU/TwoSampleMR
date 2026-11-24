context("ld")

# Opt-in: against a real OpenGWAS Data API (see test_opengwasr_live.R).
skip_if_opengwas_tests_disabled()
skip_on_cran()
skip_if_not_installed("opengwasr")
skip_if(
  Sys.getenv("OPENGWAS_URL") == "" || Sys.getenv("OPENGWAS_DEV_KEY") == "",
  "No OpenGWAS credentials (OPENGWAS_URL, OPENGWAS_DEV_KEY)."
)

test_that("clump_data() and ld_matrix() against the LD Stores", {
  src <- opengwasr::og_api()
  a <- suppressWarnings(extract_instruments("ukb-b-19953", clump = FALSE, src = src))
  a <- a[order(a$pval.exposure), ]
  a <- a[seq(1, nrow(a), length.out = 300), ]
  out <- suppressWarnings(clump_data(a, src = src))
  expect_true(nrow(out) > 0)
  expect_true(nrow(a) > nrow(out))
  expect_equal(unique(out$clump_ld_choice.exposure), "matched")
  expect_equal(unique(out$clump_ld.exposure), "hgdp1kgp-v2-EUR-k10")

  b <- suppressWarnings(ld_matrix(head(out$SNP, 5), src = src))
  expect_equal(nrow(b), ncol(b))
  expect_true(all(grepl("^[0-9XYMT]+:[0-9]+:[ACGT]+:[ACGT]+$", rownames(b))))
  expect_equal(unname(diag(b)), rep(1, nrow(b)))
})

test_that("clump_data() of several Analyses of one ancestry", {
  src <- opengwasr::og_api()
  a <- suppressWarnings(
    extract_instruments(c("ukb-b-19953", "ebi-a-GCST003116"), clump = FALSE, src = src)
  )
  a <- a[order(a$pval.exposure), ]
  a <- a[ave(seq_len(nrow(a)), a$id.exposure, FUN = seq_along) <= 100, ]
  out <- suppressWarnings(clump_data(a, src = src))
  expect_equal(length(unique(a$id.exposure)), length(unique(out$id.exposure)))
})

test_that("a refused clumping says so", {
  src <- opengwasr::og_api()
  a <- suppressWarnings(extract_instruments("ukb-b-19953", clump = FALSE, src = src))[1:50, ]
  expect_error(
    clump_data(a, clump_r2 = 0.001, clump_kb = 10000, src = src),
    class = "opengwasr_ld_refusal"
  )
})
