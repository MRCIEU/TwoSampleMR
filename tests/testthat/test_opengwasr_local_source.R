# TwoSampleMR's data functions on an opengwasr Local Source (opengwas/opengwasr#36): opengwasr's
# example Store Releases, read through ogquery, and its synthetic PLINK Reference. Nothing is mocked.
# The expected values come from opengwasr's own verbs on the same Source, mapped to TwoSampleMR's
# columns here. They skip without a Python with ogquery (RETICULATE_PYTHON), and the PLINK ones
# without a PLINK binary.

bmi_local <- "s:bmi_irn"
# Their top hits have no rsID in the example PLINK Reference, so clumping keeps none of them.
trait_a_local <- "hyb:trait_a"
outcomes_local <- c("residual-rag:a", "hyb:trait_b")

# Top hits as TwoSampleMR names them: `SNP` is the rsID, else the ALID, and a variant named twice in
# one Analysis (an rsID that names two variants) is kept once, as format_data() keeps it.
as_exposure <- function(hits) {
  hits$SNP <- ifelse(is.na(hits$rsid), hits$alid, hits$rsid)
  hits <- hits[!duplicated(paste(hits$id, hits$SNP)), , drop = FALSE]
  out <- data.frame(
    SNP = hits$SNP,
    id.exposure = hits$id,
    exposure = paste0(hits$analysis_label, " || id:", hits$id),
    chr.exposure = as.character(hits$chromosome),
    pos.exposure = as.numeric(hits$position),
    beta.exposure = hits$beta,
    se.exposure = hits$se,
    eaf.exposure = hits$eaf,
    pval.exposure = hits$pval,
    effect_allele.exposure = hits$effect_allele,
    other_allele.exposure = hits$other_allele,
    samplesize.exposure = hits$n,
    effect_scale.exposure = hits$effect_scale,
    association_status.exposure = hits$association_status,
    stringsAsFactors = FALSE
  )
  for (column in c("clump_ld", "clump_ld_choice", "clump_r2", "clump_kb")) {
    if (column %in% names(hits)) out[[paste0(column, ".exposure")]] <- hits[[column]]
  }
  sorted_rows(out)
}

sorted_rows <- function(x, id = "id.exposure") {
  x <- x[order(x[[id]], x$SNP), , drop = FALSE]
  rownames(x) <- NULL
  x
}

# An example top hit whose statistics are not finite (rs1000 of s:bmi_irn: beta NaN, se Inf) is kept
# by format_data(), its statistics cleaned to NA and `mr_keep.exposure` FALSE: compared by count.
compare_exposure <- function(got, hits) {
  want <- as_exposure(hits)
  got <- sorted_rows(as.data.frame(got))
  expect_equal(nrow(got), nrow(want))
  finite <- is.finite(want$beta.exposure) & is.finite(want$se.exposure)
  expect_equal(got$mr_keep.exposure, finite)
  got <- got[finite, names(want)]
  rownames(got) <- NULL
  got$pos.exposure <- as.numeric(got$pos.exposure)
  want <- want[finite, ]
  rownames(want) <- NULL
  expect_equal(got, want)
}

# --- Seam 1: instruments --------------------------------------------------------------------------

test_that("extract_instruments() on a Local Source clumps with its PLINK Reference at the r2 and kb given", {
  skip_if_no_local_plink()
  src <- example_local(plink = TRUE)
  ids <- c(bmi_local, trait_a_local)
  got <- suppressWarnings(suppressMessages(
    extract_instruments(ids, src = src, r2 = 0.001, kb = 10000)
  ))
  hits <- suppressMessages(opengwasr::og_tophits(
    src,
    ids,
    threshold = 5e-8,
    limit = 10000,
    clump = TRUE,
    r2 = 0.001,
    kb = 10000
  ))
  expect_gt(nrow(got), 0)
  compare_exposure(got, hits)
  expect_equal(unique(got$clump_ld.exposure), unique(hits$clump_ld))
  expect_equal(unique(got$clump_r2.exposure), 0.001)
  expect_equal(unique(got$clump_kb.exposure), 10000)
  expect_true(all(got$mr_keep.exposure))
})

test_that("extract_instruments(clump = 0) on a Local Source needs no LD resource", {
  skip_if_no_local()
  src <- example_local()
  got <- suppressWarnings(extract_instruments(bmi_local, src = src, clump = 0))
  hits <- opengwasr::og_tophits(src, bmi_local, threshold = 5e-8, limit = 10000, clump = FALSE)
  compare_exposure(got, hits)
  expect_false(any(grepl("^clump_", names(got))))
})

# --- Seam 2: outcome data -------------------------------------------------------------------------

test_that("extract_outcome_data() on a Local Source returns its Associations in TwoSampleMR's format", {
  skip_if_no_local()
  src <- example_local()
  snps <- c("rs1004", "rs1005", "1:1001500:A:G")
  got <- suppressMessages(extract_outcome_data(snps, outcomes_local, src = src))
  assoc <- opengwasr::og_associations(src, snps, outcomes_local)
  assoc$SNP <- ifelse(is.na(assoc$rsid), assoc$alid, assoc$rsid)
  assoc$SNP[assoc$alid %in% snps] <- assoc$alid[assoc$alid %in% snps]
  want <- sorted_rows(
    data.frame(
      SNP = assoc$SNP,
      id.outcome = assoc$id,
      outcome = paste0(assoc$analysis_label, " || id:", assoc$id),
      beta.outcome = assoc$beta,
      se.outcome = assoc$se,
      eaf.outcome = assoc$eaf,
      pval.outcome = assoc$pval,
      effect_allele.outcome = assoc$effect_allele,
      other_allele.outcome = assoc$other_allele,
      samplesize.outcome = assoc$n,
      effect_scale.outcome = assoc$effect_scale,
      association_status.outcome = assoc$association_status,
      stringsAsFactors = FALSE
    ),
    id = "id.outcome"
  )
  expect_gt(nrow(want), 0)
  expect_equal(sorted_rows(got[, names(want)], id = "id.outcome"), want)
  # The ALID asked for comes back as itself, an Imputed Association.
  expect_equal(
    got$association_status.outcome[got$SNP == "1:1001500:A:G" & got$id.outcome == "residual-rag:a"],
    "imputed"
  )
  expect_true(all(got$data_source.outcome == "igd"))
})

test_that("a Local Source runs the two-sample MR pipeline", {
  skip_if_no_local_plink()
  src <- example_local(plink = TRUE)
  exposure <- suppressWarnings(suppressMessages(
    extract_instruments(bmi_local, src = src, r2 = 0.001, kb = 10000)
  ))
  outcome <- suppressWarnings(suppressMessages(
    extract_outcome_data(exposure$SNP, "hyb:trait_b", src = src)
  ))
  dat <- suppressMessages(harmonise_data(exposure, outcome))
  res <- suppressMessages(mr(dat, method_list = "mr_ivw"))
  expect_equal(res$id.exposure, bmi_local)
  expect_equal(res$id.outcome, "hyb:trait_b")
  expect_equal(res$nsnp, sum(dat$mr_keep))
})

# --- Seam 3: clumping needs the Source's LD resource, and its settings ---------------------------

test_that("with a PLINK Reference and no r2 or kb, opengwasr's refusal reaches the user", {
  skip_if_no_local_plink()
  src <- example_local(plink = TRUE)
  expect_error(
    extract_instruments(bmi_local, src = src),
    class = "opengwasr_plink_settings_required"
  )
  expect_error(
    extract_instruments(bmi_local, src = src, r2 = 0.001),
    class = "opengwasr_plink_settings_required"
  )
  dat <- suppressWarnings(extract_instruments(bmi_local, src = src, clump = 0))
  expect_error(clump_data(dat, src = src), class = "opengwasr_plink_settings_required")
})

test_that("with no LD resource, clumping is refused and clump = 0 works", {
  skip_if_no_local()
  src <- example_local()
  expect_error(
    extract_instruments(bmi_local, src = src),
    class = "opengwasr_ld_resource_required"
  )
  expect_error(
    extract_instruments(bmi_local, src = src, r2 = 0.001, kb = 10000),
    class = "opengwasr_ld_resource_required"
  )
  dat <- suppressWarnings(extract_instruments(bmi_local, src = src, clump = 0))
  expect_gt(nrow(dat), 0)
  expect_error(clump_data(dat, src = src), class = "opengwasr_ld_resource_required")
  expect_error(ld_matrix(dat$SNP, src = src), class = "opengwasr_ld_resource_required")
})

test_that("clump_data() on a Local Source clumps with its PLINK Reference, as og_clump() does", {
  skip_if_no_local_plink()
  src <- example_local(plink = TRUE)
  dat <- suppressWarnings(extract_instruments(bmi_local, src = src, clump = 0))
  got <- suppressMessages(clump_data(dat, clump_r2 = 0.001, clump_kb = 10000, src = src))
  hits <- opengwasr::og_tophits(src, bmi_local, threshold = 5e-8, limit = 10000, clump = FALSE)
  # The variants of `dat`: an rsID that names two variants is kept once, as format_data() keeps it.
  name <- ifelse(is.na(hits$rsid), hits$alid, hits$rsid)
  hits <- hits[!duplicated(name), , drop = FALSE]
  leads <- suppressMessages(opengwasr::og_clump(src, hits, r2 = 0.001, kb = 10000))
  expect_setequal(got$SNP, unique(leads$rsid))
  expect_equal(unique(got$clump_ld.exposure), unique(leads$clump_ld))
  expect_equal(unique(got$clump_r2.exposure), 0.001)
})

test_that("ld_matrix() on a Local Source reads its PLINK Reference by default; an LD Store is refused", {
  skip_if_no_local_plink()
  src <- example_local(plink = TRUE)
  snps <- c("rs5000", "rs1004", "rs1005")
  got <- suppressMessages(ld_matrix(snps, src = src))
  want <- unclass(suppressMessages(opengwasr::og_ld_matrix(src, snps)))
  attr(want, "og_provenance") <- NULL
  expect_equal(got, want)
  expect_match(rownames(got), "^rs[0-9]+_[ACGT]+_[ACGT]+$")
  expect_error(ld_matrix(snps, pop = "EUR", src = src), class = "opengwasr_local_unsupported")
})

test_that("mv_extract_exposures() on a Local Source clumps with its PLINK Reference", {
  skip_if_no_local_plink()
  src <- example_local(plink = TRUE)
  ids <- c(bmi_local, "residual-rag:a")
  expect_error(mv_extract_exposures(ids, src = src), class = "opengwasr_plink_settings_required")
  got <- suppressWarnings(suppressMessages(
    mv_extract_exposures(ids, clump_r2 = 0.001, clump_kb = 10000, src = src)
  ))
  expect_setequal(unique(got$id.exposure), ids)
  expect_gt(length(unique(got$SNP)), 0)
  expect_error(
    mv_extract_exposures(ids, src = example_local()),
    class = "opengwasr_ld_resource_required"
  )
})

test_that("dat_to_MRInput()'s LD matrix comes from a Local Source's PLINK Reference, complete", {
  skip_if_no_local_plink()
  src <- example_local(plink = TRUE)
  x <- suppressWarnings(extract_instruments(bmi_local, src = src, clump = 0))
  x <- x[x$SNP %in% c("rs5000", "rs1004", "rs1005", "rs1006"), ]
  x$beta.outcome <- 0.1
  x$se.outcome <- 0.01
  x$outcome <- "o"
  out <- suppressMessages(ld_for_mr_input(x, pop = NULL, bfile = NULL, plink_bin = NULL, src = src))
  expect_false(anyNA(out$ld))
  expect_equal(nrow(out$ld), nrow(out$x))
  expect_equal(sub("_.*", "", rownames(out$ld)), out$x$SNP)
})

# --- Seam 4: legacy ids ---------------------------------------------------------------------------

test_that("a legacy id on a Local Source is refused by opengwasr, never translated through the API", {
  skip_if_no_local_plink()
  for (src in list(example_local(), example_local(plink = TRUE))) {
    # Refused before anything about clumping is checked.
    expect_error(
      extract_instruments("ukb-b-19953", src = src),
      class = "opengwasr_legacy_id_local"
    )
    expect_error(
      extract_outcome_data("rs1004", "ebi-a-GCST003116", src = src),
      class = "opengwasr_legacy_id_local"
    )
    expect_error(
      mv_extract_exposures(c(bmi_local, "ukb-b-19953"), src = src),
      class = "opengwasr_legacy_id_local"
    )
  }
  # A bare number is the legacy `ieu-a-` id it always meant, which the refusal names.
  e <- tryCatch(extract_instruments(2, src = example_local()), error = identity)
  expect_s3_class(e, "opengwasr_legacy_id_local")
  expect_equal(e$ids, "ieu-a-2")
  # Next to a canonical id, too: it is refused, not left out with a warning.
  expect_error(
    extract_instruments(c(bmi_local, "ukb-b-19953"), src = example_local(), clump = 0),
    class = "opengwasr_legacy_id_local"
  )
})

test_that("available_outcomes() and add_metadata() read a Local Source's Analyses", {
  skip_if_no_local()
  src <- example_local()
  ao <- available_outcomes(src = src)
  expect_setequal(ao$id, opengwasr::og_analyses(src)$id)
  expect_equal(ao$trait, ao$analysis_label)
  dat <- suppressWarnings(extract_instruments(bmi_local, src = src, clump = 0))
  dat <- suppressMessages(add_metadata(dat, src = src))
  expect_equal(unique(dat$effect_scale.exposure), "sd")
  expect_equal(unique(dat$units.exposure), "SD")
})
