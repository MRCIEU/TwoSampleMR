#' Get list of studies with available GWAS summary statistics through API
#'
#' Lists every Analysis (study) that the OpenGWAS Source `src` holds, as
#' [opengwasr::og_analyses()] returns them: `id` is the Analysis's name in the form
#' `<store_id>:<analysis_id>`, which the other functions accept, `analysis_label` is its trait, and
#' `assigned_ancestry`, `sample_size`, `n_cases`, `n_controls`, `stored_effect_scale`,
#' `original_sd` and the rest of the Analysis's metadata follow. For scripts that read the old
#' column, `trait` repeats `analysis_label`.
#'
#' @param opengwas_jwt Deprecated and ignored. Credentials come from `src`: [opengwasr::og_api()] reads `OPENGWAS_JWT`.
#' @param src A Source from opengwasr: [opengwasr::og_api()], or a Local Source, [opengwasr::og_local()], which reads Store Releases you have downloaded. The default, [opengwasr::og_api()] (named as TwoSampleMR's), is made only when it is needed.
#'
#' @export
#' @return Dataframe of details for all available studies
available_outcomes <- function(opengwas_jwt = NULL, src = twosamplemr_source()) {
  if (!missing(opengwas_jwt)) {
    deprecate_opengwas_jwt("available_outcomes")
  }
  catalogue(src)
}


# Extract summary statistics from MySQL db through API given a list of SNPs and outcomes
#'
#' Supply the output from [read_exposure_data()] and all the SNPs therein will be queried against the requested outcomes in remote database using API.
#'
#' The data come from opengwasr (see [opengwasr::og_associations()]).
#'
#' * `snps` are rsIDs, or ALIDs (`chromosome:position:effect allele:other allele`, on GRCh38; see
#'   [opengwasr::og_alid()]). `SNP` in the result is the rsID, or the ALID for a variant you asked
#'   for by ALID or that has no rsID.
#' * LD proxies are no longer looked up. The Store Releases opengwasr reads are Reference-Completed:
#'   they hold an Imputed Association for most variants a study did not report. Those come back
#'   as the variant itself, and `association_status.outcome` says `"imputed"` for them (and
#'   `"observed"` for the others). [harmonise_data()] reports how many were Imputed. Passing
#'   `proxies`, `rsq`, `align_alleles`, `palindromes`, `maf_threshold` or `proxy_splitsize` prints a message and
#'   changes nothing.
#' * Outcomes are Analyses, named `<store_id>:<analysis_id>`, such as `OGS-00011:GCST003116`. Legacy
#'   OpenGWAS ids such as `ebi-a-GCST003116` are accepted where they name one Analysis. An id that
#'   names no Analysis is an error.
#' * Effects of continuous traits are per SD of the trait, and `effect_scale.outcome` says so.
#' * On a Local Source ([opengwasr::og_local()]) the Associations are read from Store Releases you
#'   have downloaded. A legacy id is refused there (error class `opengwasr_legacy_id_local`).
#'
#' @param snps Array of SNP rs IDs.
#' @param outcomes Array of IDs (see \code{id} column in output from [available_outcomes()]).
#' @param proxies Deprecated and ignored: Imputed Associations replace proxies.
#' @param rsq Deprecated and ignored.
#' @param align_alleles Deprecated and ignored.
#' @param palindromes Deprecated and ignored.
#' @param maf_threshold Deprecated and ignored.
#' @param opengwas_jwt Deprecated and ignored. Credentials come from `src`: [opengwasr::og_api()] reads `OPENGWAS_JWT`.
#' @param splitsize Deprecated and ignored: opengwasr splits large queries itself.
#' @param proxy_splitsize Deprecated and ignored.
#' @param src A Source from opengwasr: [opengwasr::og_api()], or a Local Source, [opengwasr::og_local()], which reads Store Releases you have downloaded. The default, [opengwasr::og_api()] (named as TwoSampleMR's), is made only when it is needed.
#'
#' @export
#' @return Dataframe of summary statistics for all available outcomes
extract_outcome_data <- function(
  snps,
  outcomes,
  proxies = TRUE,
  rsq = 0.8,
  align_alleles = 1,
  palindromes = 1,
  maf_threshold = 0.3,
  opengwas_jwt = NULL,
  splitsize = 10000,
  proxy_splitsize = 500,
  src = twosamplemr_source()
) {
  passed <- c(
    proxies = !missing(proxies),
    rsq = !missing(rsq),
    align_alleles = !missing(align_alleles),
    palindromes = !missing(palindromes),
    maf_threshold = !missing(maf_threshold),
    proxy_splitsize = !missing(proxy_splitsize)
  )
  if (any(passed)) {
    deprecate_proxy_arguments(names(passed)[passed])
  }
  if (!missing(splitsize)) {
    deprecate_splitsize()
  }
  if (!missing(opengwas_jwt)) {
    deprecate_opengwas_jwt("extract_outcome_data")
  }

  snps <- unique(snps)
  outcomes <- resolve_analyses(src, outcomes)
  message(
    "Extracting data for ",
    length(snps),
    " SNP(s) from ",
    length(outcomes),
    " GWAS(s)"
  )
  opengwasr_outcomes(src, snps, outcomes)
}


#' Avoid issues in MR by finding impossible vals and setting to NA
#'
#' @param d Data frame
#' @return Cleaned data frame
#' @keywords internal
cleanup_outcome_data <- function(d) {
  d$se.outcome[d$se.outcome <= 0] <- NA
  d$eaf.outcome[d$eaf.outcome <= 0 | d$eaf.outcome >= 1] <- NA
  d$beta.outcome[d$beta.outcome == -9] <- NA
  return(d)
}


#' Get SE from effect size and p-value
#'
#' @param eff effect size
#' @param pval p-values
#' @return array
#' @export
get_se <- function(eff, pval) {
  abs(eff) / abs(stats::qnorm(pval / 2))
}
