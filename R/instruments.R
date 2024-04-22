#' Find instruments for use in MR from the MR Base database
#'
#' This function searches for GWAS significant SNPs (for a given p-value) for a specified set of outcomes.
#' It then performs LD based clumping to return only independent significant associations.
#'
#' @param outcomes Array of outcome IDs (see [available_outcomes()]).
#' @param p1 Significance threshold. The default is `5e-8`.
#' @param clump Logical; whether to clump results. The default is `TRUE`.
#' @param p2 Secondary clumping threshold. The default is `5e-8`.
#' @param r2 Clumping r2 cut off. The default is `0.001`.
#' @param kb Clumping distance cutoff. The default is `10000`.
#' @param opengwas_jwt Used to authenticate protected endpoints. Login to <https://api.opengwas.io> to obtain a jwt. Provide the jwt string here, or store in .Renviron under the keyname OPENGWAS_JWT.
#' @param force_server Force the analysis to extract results from the server rather than the MRInstruments package.
#'
#' @export
#' @return data frame
extract_instruments <- function(outcomes, p1 = 5e-8, clump = TRUE, p2 = 5e-8, r2 = 0.001, kb = 10000, opengwas_jwt=ieugwasr::get_opengwas_jwt(), force_server=FALSE)
{
	# .Deprecated("ieugwasr::tophits()")
	outcomes <- ieugwasr::legacy_ids(unique(outcomes))

	d <- ieugwasr::tophits(outcomes, pval=p1, clump=clump, r2=r2, kb=kb, force_server=FALSE, opengwas_jwt=opengwas_jwt)

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
