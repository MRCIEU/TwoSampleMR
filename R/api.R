

#' Toggle API address between development and release
#'
#' @param where Which API to use. Choice between "local", "release", "test". Default = "local"
#'
#' @export
#' @return NULL
toggle_api <- function(where="test")
{
	.Deprecated("ieugwasr::toggle_api")
	ieugwasr::toggle_api(where=where)
}


#' Get access token for OAuth2 access to MR Base
#'
#'
#' @export
#' @return access token string
get_mrbase_access_token <- function()
{
	.Deprecated("ieugwasr::get_access_token")
	ieugwasr::get_access_token()
}


#' Revoke access token for MR Base
#'
#' @export
#' @return NULL
revoke_mrbase_access_token <- function()
{
	.Deprecated("ieugwasr::revoke_access_token")
	ieugwasr::revoke_access_token()
}


#' Get list of studies with available GWAS summary statistics through API
#'
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data
#'
#' @export
#' @return Dataframe of details for all available studies
available_outcomes <- function(access_token = get_mrbase_access_token())
{
	.Deprecated("ieugwasr::gwasinfo")
	ieugwasr::gwasinfo(access_token=access_token)
}


# Extract summary statistics from MySQL db through API given a list of SNPs and outcomes
#'
#' Supply the output from \code{read_exposure_data} and all the SNPs therein will be queried against the requested outcomes in remote database using API.
#'
#' @param snps Array of SNP rs IDs
#' @param outcomes Array of IDs (see \code{id} column in output from \code{available_outcomes})
#' @param proxies Look for LD tags? Default is TRUE.
#' @param rsq Minimum LD rsq value (if proxies = 1). Default = 0.8.
#' @param align_alleles = 1 Try to align tag alleles to target alleles (if proxies = 1). 1 = yes, 0 = no
#' @param palindromes = 1 Allow palindromic SNPs (if proxies = 1). 1 = yes, 0 = no
#' @param maf_threshold = 0.3 MAF threshold to try to infer palindromic SNPs.
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data
#' @param ... Other options to pass to internal extraction function
#'
#' @export
#' @return Dataframe of summary statistics for all available outcomes
extract_outcome_data <- function(snps, outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token = get_mrbase_access_token(), splitsize=10000, proxy_splitsize=500)
{
	outcomes <- unique(outcomes)
	snps <- unique(snps)
	firstpass <- extract_outcome_data_internal(snps, outcomes, proxies = FALSE, access_token=access_token, splitsize = splitsize)


	if(proxies)
	{
		for(i in 1:length(outcomes))
		{
			if(class(firstpass) == "NULL")
			{
				missedsnps <- snps
			} else {
				missedsnps <- snps[!snps %in% subset(firstpass, id.outcome == outcomes[i])$SNP]
			}
			if(length(missedsnps)>0)
			{
				message("Finding proxies for ", length(missedsnps), " SNPs in outcome ", outcomes[i])
				temp <- extract_outcome_data_internal(missedsnps, outcomes[i], proxies = TRUE, rsq, align_alleles, palindromes, maf_threshold, access_token = access_token, splitsize = proxy_splitsize)
				if(!is.null(temp))
				{
					firstpass <- plyr::rbind.fill(firstpass, temp)
				}
			}
		}

	}

	return(firstpass)
}



extract_outcome_data_internal <- function(snps, outcomes, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token = get_mrbase_access_token(), splitsize=10000)
{
	snps <- unique(snps)
	message("Extracting data for ", length(snps), " SNP(s) from ", length(unique(outcomes)), " GWAS(s)")
	outcomes <- unique(outcomes)

	stopifnot(proxies %in% 0:1)

	if((length(snps) < splitsize & length(outcomes) < splitsize) | (length(outcomes) < splitsize & length(snps) < splitsize))
	{
		d <- ieugwasr::api_query('associations',
			query = list(
				id = outcomes,
				rsid = snps,
				proxies = proxies,
				r2 = rsq,
				align_alleles = align_alleles,
				palindromes = palindromes,
				maf_threshold = maf_threshold
			),
			access_token=access_token
		)
		if(!is.data.frame(d)) d <- data.frame()

	} else if(length(snps) > length(outcomes)) {

		# Split snps 
		n <- length(snps)		
		splits <- data.frame(snps=snps, chunk_id=rep(1:(ceiling(n/splitsize)), each=splitsize)[1:n])
		d <- list()
		for(i in 1:length(outcomes))
		{
			message(i, " of ", length(outcomes), " outcomes")
			
			d[[i]] <- plyr::ddply(splits, c("chunk_id"), function(x)
			{
				x <- plyr::mutate(x)
				message(" [>] ", x$chunk_id[1], " of ", max(splits$chunk_id), " chunks")
				out <- ieugwasr::api_query('associations',
					query = list(
						id = outcomes,
						rsid = snps,
						proxies = proxies,
						r2 = rsq,
						align_alleles = align_alleles,
						palindromes = palindromes,
						maf_threshold = maf_threshold
					),
					access_token=access_token
				)
				if(!is.data.frame(out)) out <- data.frame()
				return(out)
			})
		}

		d <- plyr::rbind.fill(d)

	} else {
		# Split outcomes
		n <- length(outcomes)
		splits <- data.frame(outcomes=outcomes, chunk_id=rep(1:(ceiling(n/splitsize)), each=splitsize)[1:n])
		d <- list()
		for(i in 1:length(snps))
		{
			message(i, " of ", length(snps), " snps")
			
			d[[i]] <- plyr::ddply(splits, c("chunk_id"), function(x)
			{
				x <- plyr::mutate(x)
				message(" [>] ", x$chunk_id[1], " of ", max(splits$chunk_id), " chunks")
				out <- ieugwasr::api_query('associations',
					query = list(
						id = outcomes,
						rsid = snps,
						proxies = proxies,
						r2 = rsq,
						align_alleles = align_alleles,
						palindromes = palindromes,
						maf_threshold = maf_threshold
					),
					access_token=access_token
				)
				if(!is.data.frame(out)) out <- data.frame()
				return(out)
			})
		}

		d <- plyr::rbind.fill(d)

	}
	if(is.null(nrow(d)) | nrow(d) == 0)
	{
		# message("None of the requested SNPs were available in the specified GWASs.")
		return(NULL)
	}
	d <- format_d(d)
	if (nrow(d)>0){
		d$data_source.outcome <- "mrbase"
		return(d)
	} else {
		return(NULL)
	}
}


#' Avoid issues in MR by finding impossible vals and setting to NA
#'
#' @param d Data frame
#' @return Cleaned data frame
cleanup_outcome_data <- function(d)
{
	d$se.outcome[d$se.outcome <= 0] <- NA
	d$eaf.outcome[d$eaf.outcome <= 0 | d$eaf.outcome >= 1] <- NA
	d$beta.outcome[d$beta.outcome == -9] <- NA
	return(d)
}


#' Get SE from effect size and pval
#'
#' @param eff effect size
#' @param pval pvals
#' @return array
get_se <- function(eff, pval)
{
	abs(eff) / abs(qnorm(pval / 2))
}


#' Format the returned table from the MySQL database
#'
#' @param d Data frame
#' @return Data frame
format_d <- function(d)
{
	nom <- c(
		SNP = "name",
		beta.outcome = "beta",
		se.outcome = "se",
		samplesize.outcome = "n",
		ncase.outcome = "ncase",
		ncontrol.outcome = "ncontrol",
		pval.outcome = "p",
		eaf.outcome = "effect_allele_freq",
		effect_allele.outcome = "effect_allele",
		other_allele.outcome = "other_allele",
		units.outcome = "unit",
		outcome = "trait",
		consortium.outcome = "consortium",
		year.outcome = "year",
		pmid.outcome = "pmid",
		id.outcome = "id",
		category.outcome = "category",
		subcategory.outcome = "subcategory",
		proxy.outcome = "proxy",
		target_snp.outcome = "target_snp",
		proxy_snp.outcome = "proxy_snp",
		target_a1.outcome = "target_a1",
		target_a2.outcome = "target_a2",
		proxy_a1.outcome = "proxy_a1",
		proxy_a2.outcome = "proxy_a2"
	)

	d <- subset(d, select=names(d)[names(d) %in% nom])

	index <- match(names(d), nom)
	names(d) <- names(nom)[index]
	d$originalname.outcome <- 0

	if("proxy.outcome" %in% names(d))
	{
		# If two SNPs have the same proxy SNP then one has to be removed
		d <- plyr::ddply(d, c("id.outcome"), function(x)
		{
			x <- plyr::mutate(x)
			subset(x, !duplicated(proxy_snp.outcome))
		})
	}

	if(nrow(d) == 0)
	{
		message("No matches")
		return(d)
	}

	d$originalname.outcome <- d$outcome
	d$outcome.deprecated <- paste0(d$outcome, " || ", d$consortium.outcome, " || ", d$year.outcome)
	d$outcome <- paste0(d$outcome, " || id:", d$id.outcome)

	rem <- is.na(d$beta.outcome) & is.na(d$pval.outcome)
	d <- subset(d, !rem)

	# For any that have missing SE but available beta and pval, infer SE
	index <- (is.na(d$se.outcome) | d$se.outcome == 0) & (!is.na(d$beta.outcome) & !is.na(d$pval.outcome))
	if(any(index))
	{
		d$se.outcome[index] <- get_se(d$beta.outcome[index], d$pval.outcome[index])
	}

	d <- cleanup_outcome_data(d)

	mrcols <- c("beta.outcome", "se.outcome", "effect_allele.outcome")
	d$mr_keep.outcome <- apply(d[, mrcols], 1, function(x) !any(is.na(x)))
	if(any(!d$mr_keep.outcome))
	{
		warning("The following SNP(s) are missing required information for the MR tests and will be excluded\n", paste(subset(d, !mr_keep.outcome)$SNP, collapse="\n"))
	}
	if(all(!d$mr_keep.outcome))
	{
		warning("None of the provided SNPs can be used for MR analysis, they are missing required information.")
	}

	return(d)	
}

