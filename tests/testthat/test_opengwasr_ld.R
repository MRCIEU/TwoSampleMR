# clump_data() and ld_matrix() call opengwasr (og_clump(), og_ld_matrix()), mocked here at that
# boundary; the PLINK tests run PLINK on a small synthetic reference, and skip without a binary.

bmi_dat <- function(state) {
  suppressWarnings(extract_instruments("ukb-b-19953"))
}

test_that("clump_data() clumps through og_clump(), in the LD Store matched from the Analyses, at the LD Store's own thresholds", {
  state <- local_fake_opengwasr()
  ld <- local_fake_ld()
  dat <- bmi_dat()
  out <- clump_data(dat)

  asked <- ld$calls[[1]]
  expect_equal(asked$fn, "og_clump")
  expect_null(asked$ld)
  expect_null(asked$r2)
  expect_null(asked$kb)
  expect_equal(asked$data$alid, alid_of(dat$SNP))
  expect_equal(unique(asked$data$id), bmi_id)
  expect_equal(asked$data$z, abs(qnorm(dat$pval.exposure / 2)))

  expect_equal(out$SNP, dat$SNP[seq(1, nrow(dat), by = 2)])
  expect_equal(names(out)[seq_along(names(dat))], names(dat))
  expect_equal(unique(out$clump_ld.exposure), "hgdp1kgp-v2-EUR-k10")
  expect_equal(unique(out$clump_ld_choice.exposure), "matched")
  expect_equal(unique(out$clump_r2.exposure), 0.0319)
  expect_equal(unique(out$clump_kb.exposure), 10000)
})

test_that("clump_data() sends clump_r2 and clump_kb only when they are set", {
  local_fake_opengwasr()
  ld <- local_fake_ld()
  dat <- bmi_dat()
  clump_data(dat, clump_r2 = 0.1, clump_kb = 500)
  expect_equal(ld$calls[[1]]$r2, 0.1)
  expect_equal(ld$calls[[1]]$kb, 500)
})

test_that("clump_data() leaves out variants above clump_p1", {
  local_fake_opengwasr()
  ld <- local_fake_ld()
  dat <- bmi_dat()
  cut <- sort(dat$pval.exposure)[5]
  clump_data(dat, clump_p1 = cut)
  expect_equal(nrow(ld$calls[[1]]$data), 5)
})

test_that("pop is an ancestry for the LD Store, and 'legacy' says it is gone", {
  local_fake_opengwasr()
  ld <- local_fake_ld()
  dat <- bmi_dat()
  clump_data(dat, pop = "EAS")
  expect_equal(ld$calls[[1]]$ld, "EAS")
  expect_message(clump_data(dat, pop = "legacy"), class = "twosamplemr_deprecated_argument")
  expect_equal(ld$calls[[2]]$ld, "EUR")
})

test_that("data whose ids name no Analysis are refused, asking for pop or bfile; pop may be an LD Store id", {
  local_fake_opengwasr()
  ld <- local_fake_ld()
  dat <- bmi_dat()
  dat$id.exposure <- "from-a-file"
  expect_error(clump_data(dat), class = "twosamplemr_ld_needs_pop")
  expect_error(clump_data(dat), "pop")
  expect_length(ld$calls, 0)
  clump_data(dat, pop = "EUR")
  expect_equal(ld$calls[[1]]$ld, "EUR")
  clump_data(dat, pop = "hgdp1kgp-v2-AFR-k10")
  expect_equal(ld$calls[[2]]$ld, "hgdp1kgp-v2-AFR-k10")
})

test_that("read_exposure_data(clump = TRUE) on a file refuses without pop", {
  local_fake_opengwasr()
  local_fake_ld()
  f <- system.file("extdata", "bmi.txt", package = "TwoSampleMR")
  expect_error(
    suppressMessages(read_exposure_data(f, clump = TRUE)),
    class = "twosamplemr_ld_needs_pop"
  )
})

test_that("complete_ld() sets unmeasured pairs to 0 with a message, except with a PLINK Reference", {
  m <- matrix(c(1, 0.2, NA, 0.2, 1, NA, NA, NA, 1), 3)
  dimnames(m) <- list(letters[1:3], letters[1:3])
  expect_message(out <- complete_ld(m, plink = FALSE), "2 pair.*set to 0")
  expect_equal(sum(is.na(out)), 0)
  expect_equal(out["a", "c"], 0)
  expect_equal(out["a", "b"], 0.2)
  expect_silent(complete_ld(diag(2), plink = FALSE))
  expect_silent(same <- complete_ld(m, plink = TRUE))
  expect_equal(same, m)
})

test_that("variants with no ALID are left out, with a message saying so", {
  local_fake_opengwasr()
  ld <- local_fake_ld()
  dat <- bmi_dat()
  dat$pos.exposure[2] <- NA
  dat$SNP[2] <- "not-an-rsid"
  expect_message(out <- clump_data(dat), "could not be given an ALID")
  expect_equal(nrow(ld$calls[[1]]$data), nrow(dat) - 1)
  expect_false("not-an-rsid" %in% out$SNP)
})

test_that("an rsID without a position is looked up", {
  local_fake_opengwasr()
  ld <- local_fake_ld()
  dat <- bmi_dat()
  looked_up <- fixture("tophits-bmi")[, c("alid", "rsid", "effect_allele", "other_allele")]
  dat$chr.exposure <- NULL
  dat$pos.exposure <- NULL
  testthat::local_mocked_bindings(
    og_variants_rsid = function(src, rsids, releases = NULL) looked_up[looked_up$rsid %in% rsids, ],
    .package = "opengwasr"
  )
  clump_data(dat)
  expect_equal(ld$calls[[1]]$data$alid, alid_of(dat$SNP))
})

test_that("an ALID as SNP is used as it is", {
  local_fake_opengwasr()
  ld <- local_fake_ld()
  dat <- bmi_dat()
  dat$chr.exposure <- NULL
  dat$pos.exposure <- NULL
  dat$SNP <- alid_of(dat$SNP)
  clump_data(dat)
  expect_equal(ld$calls[[1]]$data$alid, dat$SNP)
})

test_that("with bfile the PLINK Reference is the LD, and PLINK's thresholds are TwoSampleMR's old ones", {
  local_fake_opengwasr()
  ld <- local_fake_ld()
  made <- NULL
  testthat::local_mocked_bindings(
    og_plink_reference = function(path, plink_bin = NULL) {
      made <<- list(path = path, plink_bin = plink_bin)
      structure(list(path = path, plink_bin = "plink"), class = "og_plink_reference")
    },
    .package = "opengwasr"
  )
  dat <- bmi_dat()
  clump_data(dat, bfile = "/ref/EUR", plink_bin = "/bin/plink")
  expect_equal(made, list(path = "/ref/EUR", plink_bin = "/bin/plink"))
  expect_s3_class(ld$calls[[1]]$ld, "og_plink_reference")
  expect_null(ld$calls[[1]]$src)
  expect_equal(ld$calls[[1]]$r2, 0.001)
  expect_equal(ld$calls[[1]]$kb, 10000)
  expect_equal(ld$calls[[1]]$data$rsid, ifelse(grepl("^rs", dat$SNP), dat$SNP, NA_character_))
  clump_data(dat, bfile = "/ref/EUR", clump_r2 = 0.01, clump_kb = 1000)
  expect_equal(ld$calls[[2]]$r2, 0.01)
  expect_equal(ld$calls[[2]]$kb, 1000)
})

test_that("a refusal from opengwasr reaches the user", {
  local_fake_opengwasr()
  testthat::local_mocked_bindings(
    og_clump = function(...) {
      stop(structure(class = c("opengwasr_ld_refusal", "error", "condition"), list(message = "refused", call = NULL)))
    },
    .package = "opengwasr"
  )
  expect_error(clump_data(bmi_dat(), clump_r2 = 0.001, clump_kb = 10000), class = "opengwasr_ld_refusal")
})

test_that("ld_matrix() reads an LD Store, in the European one by default, and names the matrix by ALID", {
  ld <- local_fake_ld()
  alids <- c("16:53767042:C:T", "2:24913997:A:T")
  m <- ld_matrix(alids)
  expect_equal(ld$calls[[1]]$ld, "EUR")
  expect_equal(ld$calls[[1]]$variants, alids)
  expect_false(ld$calls[[1]]$fill)
  expect_equal(rownames(m), alids)
  expect_equal(colnames(m), alids)
  expect_false(inherits(m, "og_ld_matrix"))
  expect_true(is.matrix(m))
  ld_matrix(alids, pop = "SAS", fill_below_t_far = TRUE)
  expect_equal(ld$calls[[2]]$ld, "SAS")
  expect_true(ld$calls[[2]]$fill)
})

test_that("ld_matrix() says with_alleles and 'legacy' no longer apply", {
  ld <- local_fake_ld()
  expect_message(ld_matrix("16:53767042:C:T", with_alleles = FALSE), class = "twosamplemr_deprecated_argument")
  expect_message(ld_matrix("16:53767042:C:T", pop = "legacy"), class = "twosamplemr_deprecated_argument")
  expect_equal(ld$calls[[2]]$ld, "EUR")
})

test_that("harmonise_ld_dat() works on a matrix named by ALID, which it names as a PLINK matrix is", {
  local_fake_opengwasr()
  dat <- bmi_dat()[1:4, ]
  dat$beta.outcome <- 0.1
  alid <- alid_of(dat$SNP)
  m <- diag(4)
  dimnames(m) <- list(alid, alid)
  out <- harmonise_ld_dat(dat, m)
  expect_equal(nrow(out$x), 4)
  expect_equal(rownames(out$ld), dat$SNP)
  expect_equal(colnames(out$ld), dat$SNP)
})

test_that("ldsc_h2() and ldsc_rg() stop with a placeholder error naming their opengwasr ticket", {
  expect_error(ldsc_h2("ukb-b-19953"), class = "twosamplemr_not_yet_available")
  expect_error(ldsc_rg("ukb-b-19953", "ebi-a-GCST003116"), "opengwas/opengwasr#26", class = "twosamplemr_not_yet_available")
})

# --- PLINK, on the synthetic reference --------------------------------------------------------------

plink_ref_path <- function() file.path(normalizePath(test_path("fixtures", "plink")), "ref")

skip_without_plink <- function() {
  skip_if_not_installed("opengwasr")
  ok <- tryCatch(
    {
      opengwasr::og_plink_reference(plink_ref_path())
      TRUE
    },
    opengwasr_plink_not_found = function(e) FALSE
  )
  if (!ok) skip("No PLINK binary found (genetics.binaRies, or plink on the PATH).")
}

plink_dat <- function() {
  bim <- utils::read.table(
    paste0(plink_ref_path(), ".bim"),
    col.names = c("chr", "rsid", "cm", "bp", "a1", "a2"),
    stringsAsFactors = FALSE
  )
  chr16 <- bim[bim$chr == 16 & bim$bp < 54000000, ]
  data.frame(
    SNP = chr16$rsid,
    chr.exposure = as.character(chr16$chr),
    pos.exposure = chr16$bp,
    effect_allele.exposure = chr16$a1,
    other_allele.exposure = chr16$a2,
    pval.exposure = c(1e-30, 1e-20, 1e-10)[seq_len(nrow(chr16))],
    id.exposure = "mine",
    stringsAsFactors = FALSE
  )
}

test_that("PLINK: clump_data(bfile =) keeps the strongest variant of variants in LD", {
  skip_without_plink()
  local_no_credentials()
  dat <- plink_dat()
  expect_equal(nrow(dat), 3)
  out <- suppressMessages(clump_data(dat, bfile = plink_ref_path()))
  expect_equal(out$SNP, "rs1421085")
  expect_equal(out$clump_ld_choice.exposure, "user")
  expect_equal(out$clump_r2.exposure, 0.001)
  expect_equal(out$clump_kb.exposure, 10000)
})

test_that("PLINK: ld_matrix(bfile =) is named <rsid>_<a1>_<a2>", {
  skip_without_plink()
  local_no_credentials()
  m <- suppressMessages(ld_matrix(c("rs1421085", "rs1558902"), bfile = plink_ref_path()))
  expect_equal(nrow(m), 2)
  expect_true(all(grepl("^rs[0-9]+_[ACGT]_[ACGT]$", rownames(m))))
  expect_equal(unname(diag(m)), c(1, 1))
  expect_false(is.na(m[1, 2]))
})

test_that("PLINK: every variant of the reference is clumped and matched by rsID, with no API", {
  skip_without_plink()
  local_no_credentials()
  bim <- utils::read.table(paste0(plink_ref_path(), ".bim"), stringsAsFactors = FALSE)
  dat <- data.frame(
    SNP = bim$V2,
    pval.exposure = seq(1e-30, 1e-10, length.out = nrow(bim)),
    id.exposure = "mine",
    stringsAsFactors = FALSE
  )
  dat$SNP[1] <- "rs0000000" # not in the reference
  expect_message(out <- clump_data(dat, bfile = plink_ref_path()), "1")
  expect_false("rs0000000" %in% out$SNP)
  expect_true(nrow(out) > 1)
  expect_true(all(out$SNP %in% bim$V2))
})

test_that("clump_data() clumps each id on its own, as before, each in its own LD Store", {
  local_fake_opengwasr()
  ld <- local_fake_ld(keep = "all")
  dat <- bmi_dat()
  other <- dat[1:5, ]
  other$id.exposure <- cad_id
  both <- rbind(dat, other)
  out <- clump_data(both)
  expect_length(ld$calls, 2)
  expect_setequal(vapply(ld$calls, function(x) unique(x$data$id), ""), c(bmi_id, cad_id))
  expect_equal(nrow(out), nrow(both))
  expect_equal(out$SNP, both$SNP)
})

test_that("a mix of Analyses and own data says which ids have no Analysis", {
  local_fake_opengwasr()
  local_fake_ld()
  dat <- bmi_dat()
  dat$id.exposure[1:2] <- "mine"
  expect_error(clump_data(dat), "mine", class = "twosamplemr_ld_needs_pop")
})

test_that("the PLINK path through dat_to_MRInput() is not filled, the LD Store path is", {
  skip_without_plink()
  local_no_credentials()
  x <- plink_dat()
  x$beta.exposure <- 0.1
  x$beta.outcome <- 0.1
  x$se.outcome <- 0.01
  x$exposure <- "e"
  x$outcome <- "o"
  out <- suppressMessages(ld_for_mr_input(x, pop = "EUR", bfile = plink_ref_path(), plink_bin = NULL, src = NULL))
  expect_false(anyNA(out$ld))
  expect_equal(rownames(out$ld), x$SNP)
})

test_that("the LD Store path through dat_to_MRInput() asks for fill_below_t_far and zeros the rest", {
  local_fake_opengwasr()
  dat <- bmi_dat()[1:3, ]
  seen <- NULL
  testthat::local_mocked_bindings(
    og_ld_matrix = function(src, variants, ld = NULL, fill_below_t_far = FALSE) {
      seen <<- fill_below_t_far
      v <- alid_of(dat$SNP)
      m <- diag(3)
      m[1, 2] <- m[2, 1] <- NA
      dimnames(m) <- list(v, v)
      structure(m, class = c("og_ld_matrix", "matrix", "array"))
    },
    .package = "opengwasr"
  )
  expect_message(out <- ld_for_mr_input(dat, pop = "EUR", bfile = NULL, plink_bin = NULL, src = NULL), "set to 0")
  expect_true(seen)
  expect_false(anyNA(out$ld))
})

test_that("TwoSampleMR's Sources name it as the caller", {
  seen <- NULL
  testthat::local_mocked_bindings(
    og_api = function(..., caller = NULL) {
      seen <<- caller
      structure(list(), class = "fake_source")
    },
    .package = "opengwasr"
  )
  twosamplemr_source()
  expect_match(seen, "^TwoSampleMR/0\\.8\\.0$")
})
