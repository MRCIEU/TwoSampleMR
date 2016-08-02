#' Find instruments for use in MR from the MR Base database
#'
#' This function searches for GWAS significant SNPs (for a given p-value) for a specified set of outcomes.
#' It then performs LD based clumping to return only independent significant associations.
#'
#' @param outcomes Array of outcome IDs (see \code{available_outcomes})
#' @param p1 = 5e-8 Significance threshold
#' @param clump = TRUE Whether to clump results
#' @param p2 = 5e-8 Secondary clumping threshold
#' @param r2 = 0.1 Clumping r2 cut off
#' @param kb = 5000 Clumping distance cutoff
#' @param access_token = get_mrbase_access_token() Google OAuth2 access token. Used to authenticate level of access to data
#'
#' @export
#' @return data frame
extract_instruments <- function(outcomes, p1 = 5e-8, clump = TRUE, p2 = 5e-8, r2 = 0.001, kb = 10000, access_token = get_mrbase_access_token())
{
	outcomes <- unique(outcomes)
	message("Extracting data from ", length(unique(outcomes)), " GWAS(s)")
	if(clump) message("and performing LD clumping")

	# outcomes <- paste(outcomes, collapse=",")

	d <- list()
	for(i in 1:length(outcomes))
	{
		message(" [>] ", i, " of ", length(outcomes))
		url <- paste0(options()$mrbaseapi, "extract_instruments?access_token=", access_token,
			"&outcomes=", outcomes[i], 
			"&pval=", p1,
			"&clump=", ifelse(clump, "yes", "no"),
			"&p2=", p2,
			"&r2=", r2,
			"&kb=", kb
		)
		out <- fromJSON_safe(url)
		if(!is.data.frame(out)) out <- data.frame()
		d[[i]] <- out
	}
	d <- plyr::rbind.fill(d)

	if(length(d) == 0)
	{
		message("None of the requested outcomes had GWAS hits at the specified threshold.")
		return(NULL)
	}



	d$phenotype <- paste0(d$trait, " || ", d$consortium, " || ", d$year, " || ", d$unit)
	d <- format_data(
		d,
		type="exposure",
		snps=NULL,
		phenotype_col="phenotype",
		snp_col="name",
		beta_col="beta",
		se_col="se",
		eaf_col="effect_allelel_freq",
		effect_allele_col="effect_allele",
		other_allele_col="other_allele",
		pval_col="p",
		units_col="unit",
		ncase_col="ncase",
		ncontrol_col="ncontrol",
		samplesize_col="n",
		min_pval=1e-200
	)
	d$data_source.exposure <- "mrbase"
	return(d)
}
