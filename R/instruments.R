#' Find instruments for use in MR from the OpenGWAS database
#'
#' This function searches for GWAS significant SNPs (for a given p-value) for a specified set of outcomes.
#' It then performs LD based clumping to return only independent significant associations.
#'
#' The data come from opengwasr (see [opengwasr::og_tophits()]).
#'
#' * Analyses are named `<store_id>:<analysis_id>`, such as `OGS-00009:ukb-b-19953`, and that is
#'   what `id.exposure` holds. Legacy OpenGWAS ids such as `ukb-b-19953` are accepted where they name
#'   one Analysis, and an id that names no Analysis is an error.
#' * Effects of continuous traits are per SD of the trait, and `effect_scale.exposure` says so
#'   (`"sd"` or `"log_or"`). Estimates from [mr()] then are per SD of the exposure.
#' * `SNP` is the rsID, or the ALID (`chromosome:position:effect allele:other allele`, on GRCh38) of a
#'   variant that has no rsID.
#' * Clumping runs on the OpenGWAS server, in the LD Store matching each Analysis's assigned
#'   ancestry, and at that LD Store's own thresholds unless you set `r2` and `kb`. The LD Store used
#'   and the thresholds applied are in `clump_ld.exposure`, `clump_ld_choice.exposure`,
#'   `clump_r2.exposure` and `clump_kb.exposure`. TwoSampleMR's old thresholds, `r2 = 0.001` within
#'   `kb = 10000`, are refused by the LD Stores, so they are not the default.
#' * `association_status.exposure` says whether each Association was observed in the study or
#'   Imputed by Reference Completion.
#'
#' @param outcomes Array of outcome IDs (see [available_outcomes()]).
#' @param p1 Significance threshold. The default is `5e-8`. OpenGWAS indexes `5e-8`, `5e-6` and `5e-4`.
#' @param clump Whether to clump results (`1`) or not (`0`). Default is `1`. (`TRUE` and `FALSE` are also allowed for backwards compatibility.)
#' @param r2 Clumping r2 cut off. The default, `NULL`, is the LD Store's own threshold (about `0.03`).
#' @param kb Clumping distance cutoff. The default, `NULL`, is the LD Store's own window (`10000`).
#' @param opengwas_jwt Deprecated and ignored. Credentials come from `src`: [opengwasr::og_api()] reads `OPENGWAS_JWT`.
#' @param force_server Deprecated and ignored. Clumping always runs on the server.
#' @param src A Source from opengwasr. The default, [opengwasr::og_api()], is made only when it is needed.
#'
#' @export
#' @return data frame
extract_instruments <- function(
  outcomes,
  p1 = 5e-8,
  clump = 1,
  r2 = NULL,
  kb = NULL,
  opengwas_jwt = NULL,
  force_server = FALSE,
  src = opengwasr::og_api()
) {
  if (!(clump %in% c(0, 1, FALSE, TRUE))) {
    stop("The clump argument should be 0 or 1.")
  }
  if (!missing(opengwas_jwt)) {
    deprecate_opengwas_jwt("extract_instruments")
  }
  if (!missing(force_server)) {
    deprecated_argument_message(paste0(
      "`force_server` is no longer used by `extract_instruments()` and is ignored. ",
      "Clumping always runs on the OpenGWAS server."
    ))
  }

  d <- opengwasr::og_tophits(
    src,
    as.character(unique(outcomes)),
    threshold = p1,
    limit = 10000,
    clump = as.logical(clump),
    r2 = r2,
    kb = kb
  )

  if (nrow(d) == 0) {
    return(NULL)
  }
  d$SNP <- variant_name(d)
  d$phenotype <- paste0(d$analysis_label, " || id:", d$id)
  key <- paste(d$id, tolower(d$SNP))
  res <- format_data(
    d,
    type = "exposure",
    snps = NULL,
    phenotype_col = "phenotype",
    snp_col = "SNP",
    chr_col = "chromosome",
    pos_col = "position",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col = "n",
    ncase_col = "ncase",
    ncontrol_col = "ncontrol",
    min_pval = 1e-200,
    id_col = "id"
  )
  index <- match(paste(res$id.exposure, res$SNP), key)
  # format_data() lower-cases SNP, which would break an ALID.
  alid <- is.na(d$rsid[index])
  res$SNP[alid] <- d$alid[index][alid]
  res <- carry_columns(
    res,
    d,
    index,
    c("effect_scale", "association_status", "clump_ld", "clump_ld_choice", "clump_r2", "clump_kb"),
    "exposure"
  )
  res$data_source.exposure <- "igd"
  return(res)
}
