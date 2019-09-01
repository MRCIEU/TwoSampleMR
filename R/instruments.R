#' Find instruments for use in MR from the MR Base database
#'
#' This function searches for GWAS significant SNPs (for a given p-value) for a specified set of outcomes.
#' It then performs LD based clumping to return only independent significant associations.
#'
#' @param outcomes Array of outcome IDs (see \code{available_outcomes})
#' @param p1 = 5e-8 Significance threshold
#' @param clump = 1 Whether to clump results
#' @param p2 = 5e-8 Secondary clumping threshold
#' @param r2 = 0.001 Clumping r2 cut off
#' @param kb = 10000 Clumping distance cutoff
#' @param access_token = get_mrbase_access_token() Google OAuth2 access token. Used to authenticate level of access to data
#' @param force_server Force the analysis to extract results from the server rather than the MRInstruments package
#'
#' @export
#' @return data frame
extract_instruments <- function(outcomes, p1 = 5e-8, clump = 1, p2 = 5e-8, r2 = 0.001, kb = 10000, access_token = get_mrbase_access_token(), force_server=FALSE)
{
	outcomes <- unique(outcomes)

	if(clump & p1 == 5e-8 & r2 == 0.001 & kb == 10000 & !force_server)
	{
		message("Requesting default values. Extracting from pre-clumped data")
		a <- require(MRInstruments)
		if(!a)
		{
			message("MRInstruments package not available")
			message("To install: devtools::install_github('MRCIEU/MRInstruments')")
			message("and then try again")
			return(NULL)
		}

		data(mrbase_instruments, envir=environment())
		a <- exists("mrbase_instruments")
		if(!a)
		{
			message("Pre-clumped dataset is not available. You might have an old version of the MRInstruments package")
			message("To update: devtools::install_github('MRCIEU/MRInstruments')")
			message("and then try again")
			return(NULL)
		}

		a <- subset(mrbase_instruments, id.exposure %in% outcomes)

		if(nrow(a) == 0)
		{
			message("None of the requested outcomes had GWAS hits at the specified threshold.")
			return(NULL)
		}
		a$exposure <- paste0(a$exposure, " || id:", a$id)
		return(a)
	}

	message("Extracting data from ", length(unique(outcomes)), " GWAS(s)")
	if(clump) message("and performing LD clumping")

	# outcomes <- paste(outcomes, collapse=",")

	d <- list()
	for(i in 1:length(outcomes))
	{
		message(" [>] ", i, " of ", length(outcomes))		
		out <- api_query('tophits', query = list(id = outcomes[i], pval = p1, clump = as.numeric(clump), r2 = r2, kb = kb), access_token=access_token)
		if(!is.data.frame(out)) out <- data.frame()
		d[[i]] <- out
	}
	d <- plyr::rbind.fill(d)

	if(length(d) == 0)
	{
		message("None of the requested outcomes had GWAS hits at the specified threshold.")
		return(NULL)
	}

	# d$phenotype.deprecated <- paste0(d$trait, " || ", d$consortium, " || ", d$year, " || ", d$unit)
	d$phenotype <- paste0(d$trait, " || id:", d$id)
	if("ncase" %in% names(d)) d$ncase <- as.numeric(d$ncase)
	if("ncontrol" %in% names(d)) d$ncontrol <- as.numeric(d$ncontrol)
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
		min_pval=1e-200,
		id_col="id"
	)
	d$data_source.exposure <- "mrbase"

	return(d)
}
