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
#' @param access_token = get_mrbase_access_token() Google OAuth2 access token. Used to authenticate level of access to data
#' @param force_server Force the analysis to extract results from the server rather than the MRInstruments package
#' @param force_server_if_empty Some of the newly added MR-Base datasets don't have pre-calculated clumped results yet. This option is soon to be deprecated but temporarily we are forcing the search for instruments when an outcome doesn't have any precalculated results. Default = TRUE.
#'
#' @export
#' @return data frame
extract_instruments <- function(outcomes, p1 = 5e-8, clump = TRUE, p2 = 5e-8, r2 = 0.001, kb = 10000, access_token = get_mrbase_access_token(), force_server=FALSE, force_server_if_empty=TRUE)
{
	outcomes <- unique(outcomes)

	# ieu_index <- ! grepl("[A-Z]+-[a-z]+:", outcomes)
	# outcomes[ieu_index] <- paste0("IEU-a:", outcomes[ieu_index])

	# TODO: UKB-b traits don't have MRInstruments entry

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
			if(!force_server_if_empty)
			{
				message("None of the requested outcomes had GWAS hits at the specified threshold.")
				return(NULL)
			}
		}
		if(nrow(a) > 0)
		{
			a$id <- gsub("IEU-a:", "", a$id)
			a$exposure <- paste0(a$exposure, " || id:", a$id)
		}
		missing_outcomes <- outcomes[!outcomes %in% a$id]
		message(length(outcomes) - length(missing_outcomes), " out of ", length(outcomes), " requested outcomes have pre-calculated instruments. ")
		if(length(missing_outcomes) == 0)
		{
			return(a)
		} else {
			if(!force_server_if_empty)
			{
				message("force_server_if_empty=FALSE, so not checking server for outcomes with no instruments")
				return(a)
			} else {
				found_outcomes <- outcomes[outcomes %in% a$id]
				outcomes <- missing_outcomes
			}
		}
	}

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

	# d$phenotype.deprecated <- paste0(d$trait, " || ", d$consortium, " || ", d$year, " || ", d$unit)
	d$phenotype <- paste0(d$trait, " || id:", d$id)
	d$ncase <- as.numeric(d$ncase)
	d$ncontrol <- as.numeric(d$ncontrol)
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
	d$id.exposure <- gsub("IEU-a:", "", d$id.exposure)

	if(force_server_if_empty)
	{
		d <- plyr::rbind.fill(d, a)
	}

	return(d)
}
