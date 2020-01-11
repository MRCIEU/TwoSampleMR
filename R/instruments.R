#' Find instruments for use in MR from the MR Base database
#'
#' This function searches for GWAS significant SNPs (for a given p-value) for a specified set of outcomes.
#' It then performs LD based clumping to return only independent significant associations.
#'
#' @param outcomes Array of outcome IDs (see \code{available_outcomes})
#' @param p1 = 5e-8 Significance threshold
#' @param clump = TRUE Whether to clump results
#' @param p2 = 5e-8 Secondary clumping threshold
#' @param r2 = 0.001 Clumping r2 cut off
#' @param kb = 10000 Clumping distance cutoff
#' @param access_token = ieugwasr::check_access_token() Google OAuth2 access token. Used to authenticate level of access to data
#' @param force_server Force the analysis to extract results from the server rather than the MRInstruments package
#' @param force_server_if_empty Some of the newly added MR-Base datasets don't have pre-calculated clumped results yet. This option is soon to be deprecated but temporarily we are forcing the search for instruments when an outcome doesn't have any precalculated results. Default = TRUE.
#'
#' @export
#' @return data frame
extract_instruments <- function(outcomes, p1 = 5e-8, clump = TRUE, p2 = 5e-8, r2 = 0.001, kb = 10000, access_token = ieugwasr::check_access_token(), force_server=FALSE)
{
	# .Deprecated("ieugwasr::tophits()")
	outcomes <- ieugwasr::legacy_ids(unique(outcomes))

	d <- ieugwasr::tophits(outcomes, pval=p1, clump=clump, r2=r2, kb=kb, force_server=FALSE, access_token=access_token)

	# d$phenotype.deprecated <- paste0(d$trait, " || ", d$consortium, " || ", d$year, " || ", d$unit)
	if(nrow(d) == 0) return(NULL)
	d$phenotype <- paste0(d$trait, " || id:", d$id)
	d <- format_data(
		d,
		type="exposure",
		snps=NULL,
		phenotype_col="phenotype",
		snp_col="rsid",
		chr_col="chr",
		pos_col="position",
		beta_col="beta",
		se_col="se",
		eaf_col="eaf",
		effect_allele_col="ea",
		other_allele_col="nea",
		pval_col="p",
		samplesize_col="n",
		min_pval=1e-200,
		id_col="id"
	)
	d$data_source.exposure <- "igd"
	return(d)
}
