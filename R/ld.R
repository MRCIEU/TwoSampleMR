#' Perform LD clumping on SNP data
#'
#' Clumps `dat` so that one variant is kept for each independent signal: the strongest variant, and
#' every weaker variant within a window whose LD with it is above a threshold is dropped. The work
#' is done by [opengwasr::og_clump()].
#'
#' @details
#' By default the LD Store that matches the Assigned Ancestry of the Analyses in `dat` (its
#' `id.exposure`) is used, and `clump_r2` and `clump_kb` are the LD Store's own thresholds (r2 about
#' 0.03 within 10,000 kb). They are sent only when you set them: LD Stores refuse r2 = 0.001
#' within 10,000 kb, so scripts that pass those old defaults now get opengwasr's refusal. If the Analyses
#' differ in ancestry, or one is Unassigned, opengwasr refuses unless you choose the LD with `pop`.
#' Data that do not come from OpenGWAS (whose `id.exposure` names no Analysis) have no Assigned Ancestry, so they are refused (error class `twosamplemr_ld_needs_pop`) unless you give `pop` or `bfile`, or the Source carries an LD resource of its own: add `pop = "EUR"` to an old script.
#'
#' Variants are ranked by their p-value (as a |z| score). Each variant is given an ALID
#' (`chromosome:position:allele:allele`, GRCh38, see [opengwasr::og_alid()]): its `SNP`, if that is
#' one; else made from the `chr.exposure` and `pos.exposure` columns (GRCh38 positions) and the
#' alleles; else looked up from its rsID (see [opengwasr::og_variants_rsid()]). Variants with no ALID,
#' or that the LD Store does not hold, cannot be clumped and are left out of the result
#' (opengwasr warns about the latter).
#'
#' With `bfile`, PLINK clumps on your machine using your own reference ([opengwasr::og_plink_reference()]),
#' matching variants by rsID. There the thresholds have no defaults to take, so `clump_r2 = 0.001` and
#' `clump_kb = 10000` (TwoSampleMR's old defaults) are used unless you set them.
#'
#' On a Local Source ([opengwasr::og_local()]) the default is the Source's own LD resource: PLINK
#' clumps with the PLINK Reference it was given as `bfile`, matching variants by rsID. It has no
#' thresholds of its own, so `clump_r2` and `clump_kb` must be set (opengwasr refuses otherwise, with
#' error class `opengwasr_plink_settings_required`). A Local Source with no LD resource is refused
#' (`opengwasr_ld_resource_required`), and so is a `pop`, as a Local Source never reads the API's LD
#' Stores; `bfile` here works as on the API Source.
#'
#' Each `id.exposure` is clumped on its own, in the LD Store of its own ancestry, as the old OpenGWAS client did.
#' To check which of your variants an LD Store holds first, use [opengwasr::og_ld_reflookup()].
#'
#' The result has `clump_ld.exposure`, `clump_ld_choice.exposure`, `clump_r2.exposure` and
#' `clump_kb.exposure` columns recording the LD used (the path of the reference with `bfile`).
#'
#' @param dat Output from [format_data()]. Must have a SNP name column (`SNP`), which is an rsID or an ALID. If `id.exposure` or `pval.exposure` not present they will be generated.
#' @param clump_kb Clumping window in kb. The default, `NULL`, is the LD Store's own (10,000).
#' @param clump_r2 Clumping r2 cutoff. The default, `NULL`, is the LD Store's own.
#' @param clump_p1 Clumping sig level for index SNPs, default is `1`: variants with a larger p-value are left out.
#' @param pop Ancestry (`"EUR"`, `"SAS"`, `"EAS"`, `"AFR"` or `"AMR"`) or LD Store id (see `opengwasr::og_ld_stores()`) of the LD Store to use. The default, `NULL`, is the Source's own LD resource if it carries one (`opengwasr::og_api(ld_store = )`, `opengwasr::og_local(bfile = )`), else matches the Assigned Ancestry of the Analyses in `dat`. `"legacy"` is deprecated and means `"EUR"`.
#' @param bfile Prefix of a PLINK fileset (`.bed`, `.bim`, `.fam`) to clump with locally instead of an LD Store. Default = `NULL`
#' @param plink_bin If `NULL` and `bfile` is not `NULL` then will detect packaged plink binary for specific OS. Otherwise specify path to plink binary. Default = `NULL`
#' @param src A Source from opengwasr: [opengwasr::og_api()], or a Local Source, [opengwasr::og_local()], which reads Store Releases you have downloaded. The default, [opengwasr::og_api()] (named as TwoSampleMR's), is made only when it is needed.
#'
#' @export
#' @return Data frame
clump_data <- function(
  dat,
  clump_kb = NULL,
  clump_r2 = NULL,
  clump_p1 = 1,
  pop = NULL,
  bfile = NULL,
  plink_bin = NULL,
  src = twosamplemr_source()
) {
  clump_data_impl(dat, clump_kb, clump_r2, clump_p1, pop, bfile, plink_bin, src, joint = FALSE)
}

# `clump_data()`, each `id.exposure` on its own (as the old OpenGWAS client did), or with `joint = TRUE` all the ids
# as one set (the clumping across exposures of multivariable MR).
clump_data_impl <- function(dat, clump_kb, clump_r2, clump_p1, pop, bfile, plink_bin, src, joint) {
  pval_column <- "pval.exposure"

  if (!is.data.frame(dat)) {
    stop("Expecting data frame returned from format_data().")
  }

  if ("pval.exposure" %in% names(dat) && "pval.outcome" %in% names(dat)) {
    message("pval.exposure and pval.outcome columns present. Using pval.exposure for clumping.")
  } else if (!"pval.exposure" %in% names(dat) && "pval.outcome" %in% names(dat)) {
    message("pval.exposure column not present, using pval.outcome column for clumping.")
    pval_column <- "pval.outcome"
  } else if (!"pval.exposure" %in% names(dat)) {
    message("pval.exposure not present, setting clumping p-value to 0.99 for all variants")
    dat$pval.exposure <- 0.99
  } else {
    pval_column <- "pval.exposure"
  }

  if (!"id.exposure" %in% names(dat)) {
    dat$id.exposure <- random_string(1)
  }

  sets <- if (joint) {
    list(seq_len(nrow(dat)))
  } else {
    unname(split(seq_len(nrow(dat)), as.character(dat$id.exposure)))
  }
  kept <- lapply(sets, function(rows) {
    set <- dat[rows, , drop = FALSE]
    set$.original_row <- rows
    opengwasr_clump(
      src,
      set,
      pval_column,
      clump_p1,
      ld = ld_choice(src, pop, bfile, plink_bin, set$id.exposure),
      r2 = clump_r2,
      kb = clump_kb
    )
  })
  out <- dplyr::bind_rows(kept)
  out <- out[order(out$.original_row), , drop = FALSE]
  out$.original_row <- NULL
  rownames(out) <- NULL
  out
}


#' Get LD matrix for list of SNPs
#'
#' This function takes a list of SNPs and returns the matrix of LD r values (signed, and not
#' squared) between every pair, read from an OpenGWAS LD Store ([opengwasr::og_ld_matrix()]), or
#' computed by PLINK from your own reference with `bfile`.
#'
#' @details
#' The matrix is named by ALID (`chromosome:position:allele:allele`, GRCh38) from an LD Store, and
#' the sign of r refers to the first allele in each name. From a PLINK reference it is named
#' `<rsid>_<allele1>_<allele2>` with the reference's alleles, as the old OpenGWAS client named a matrix made locally,
#' and the sign of r refers to the first allele. Either way the names carry the alleles, so
#' `with_alleles` has no effect.
#'
#' To check which variants an LD Store holds first, use [opengwasr::og_ld_reflookup()].
#'
#' On a Local Source ([opengwasr::og_local()]) the matrix comes from the Source's own LD resource, the
#' PLINK Reference it was given as `bfile` (named `<rsid>_<allele1>_<allele2>`, the rsIDs of ALIDs
#' looked up in its Store Releases). A Local Source with no LD resource is refused by opengwasr (error
#' class `opengwasr_ld_resource_required`), and so is a `pop`, as a Local Source never reads the API's
#' LD Stores.
#'
#' A variant the LD Store (or the reference) does not hold is left out of the matrix, with a warning
#' (a message for a reference). An LD Store does not hold every pair: a pair it does not hold is `NA`
#' (never 0), unless `fill_below_t_far = TRUE`, which sets to 0 the pairs it holds only as being
#' below its Far Threshold (see [opengwasr::og_ld_matrix()]). Pairs further apart than its far window,
#' or on different chromosomes, remain `NA`.
#'
#' @param snps List of SNPs: rsIDs or ALIDs.
#' @param with_alleles Deprecated and ignored: the names carry the alleles.
#' @param pop Ancestry (`"EUR"`, `"SAS"`, `"EAS"`, `"AFR"` or `"AMR"`) or LD Store id (see `opengwasr::og_ld_stores()`) of the LD Store to use. The default, `NULL`, is the Source's own LD resource if it carries one (`opengwasr::og_api(ld_store = )`, `opengwasr::og_local(bfile = )`), else the European LD Store on the API Source (see Details). `"legacy"` is deprecated and means `"EUR"`. Ignored with `bfile`.
#' @param bfile Prefix of a PLINK fileset (`.bed`, `.bim`, `.fam`) to use instead of an LD Store. Default = `NULL`
#' @param plink_bin If `NULL` and `bfile` is not `NULL` then will detect packaged plink binary for specific OS. Otherwise specify path to plink binary. Default = `NULL`
#' @param fill_below_t_far Set pairs below the LD Store's Far Threshold to 0 (default `FALSE`).
#' @param src A Source from opengwasr: [opengwasr::og_api()], or a Local Source, [opengwasr::og_local()], which reads Store Releases you have downloaded. The default, [opengwasr::og_api()] (named as TwoSampleMR's), is made only when it is needed.
#'
#' @export
#' @return Matrix of LD r values
ld_matrix <- function(
  snps,
  with_alleles = TRUE,
  pop = NULL,
  bfile = NULL,
  plink_bin = NULL,
  fill_below_t_far = FALSE,
  src = twosamplemr_source()
) {
  if (!missing(with_alleles)) {
    deprecated_argument_message(paste0(
      "`with_alleles` is no longer used by `ld_matrix()` and is ignored: ",
      "the names of the matrix always carry the alleles."
    ))
  }
  # The default Source is the API's, with no LD resource of its own, so it is made only when the LD
  # Store is asked; `own_ld` is a promise, which `bfile` and `pop` leave unread.
  src_given <- !missing(src)
  ld <- ld_choice(
    src,
    pop,
    bfile,
    plink_bin,
    ids = NULL,
    own_ld = src_given && defers_to_source_ld(src)
  )
  opengwasr_ld_matrix(src, snps, ld, fill_below_t_far)
}
