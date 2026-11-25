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
#' * On a Local Source ([opengwasr::og_local()]), the data are read from Store Releases you have
#'   downloaded, and nothing calls the API. Clumping uses the Source's own LD resource: PLINK, with
#'   the PLINK Reference given to `og_local()` as `bfile`, matching variants by rsID. A PLINK
#'   Reference has no thresholds of its own, so `r2` and `kb` must be set (`r2 = 0.001, kb = 10000`
#'   are TwoSampleMR's old ones); without them opengwasr refuses, with error class
#'   `opengwasr_plink_settings_required`. A Local Source with no LD resource can't clump (error class
#'   `opengwasr_ld_resource_required`): give it one, or use `clump = 0`. Store Releases hold canonical
#'   ids only, so a legacy id is refused (error class `opengwasr_legacy_id_local`): translate it once
#'   with `opengwasr::og_analyses(opengwasr::og_api(), ids)`.
#'
#' @param outcomes Array of outcome IDs (see [available_outcomes()]).
#' @param p1 Significance threshold. The default is `5e-8`. OpenGWAS indexes `5e-8`, `5e-6` and `5e-4`.
#' @param clump Whether to clump results (`1`) or not (`0`). Default is `1`. (`TRUE` and `FALSE` are also allowed for backwards compatibility.) Unclumped, or clumped, you get at most 10,000 variants per Analysis (opengwasr warns if there were more).
#' @param r2 Clumping r2 cut off. The default, `NULL`, is the LD Store's own threshold (about `0.03`). Required with a PLINK Reference (a Local Source's `bfile`).
#' @param kb Clumping distance cutoff. The default, `NULL`, is the LD Store's own window (`10000`). Required with a PLINK Reference (a Local Source's `bfile`).
#' @param opengwas_jwt Deprecated and ignored. Credentials come from `src`: [opengwasr::og_api()] reads `OPENGWAS_JWT`.
#' @param force_server Deprecated and ignored. Clumping always runs on the server.
#' @param src A Source from opengwasr: [opengwasr::og_api()], or a Local Source, [opengwasr::og_local()], which reads Store Releases you have downloaded. The default, [opengwasr::og_api()] (named as TwoSampleMR's), is made only when it is needed.
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
  src = twosamplemr_source()
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

  # Resolved first, so that an id opengwasr refuses (a legacy id on a Local Source) is the error
  # the user sees, before anything about clumping is checked.
  analyses <- resolve_analyses(src, outcomes)
  opengwasr_instruments(
    src,
    analyses,
    p1 = p1,
    clump = clump,
    r2 = r2,
    kb = kb
  )
}
