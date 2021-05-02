
#' Get list of studies with available GWAS summary statistics through API
#'
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data
#'
#' @export
#' @return Dataframe of details for all available studies
available_outcomes <- function(access_token = ieugwasr::check_access_token())
{
	# .Deprecated("ieugwasr::gwasinfo()")
	a <- ieugwasr::gwasinfo(access_token=access_token)	
	return(a)
}



# Extract summary statistics from MySQL db through API given a list of SNPs and outcomes
#'
#' Supply the output from \code{\link{read_exposure_data}} and all the SNPs therein will be queried against the requested outcomes in remote database using API.
#'
#' @md
#' @param snps Array of SNP rs IDs.
#' @param outcomes Array of IDs (see \code{id} column in output from \code{\link{available_outcomes}}).
#' @param proxies Look for LD tags? Default is `TRUE`.
#' @param rsq Minimum LD rsq value (if proxies = 1). Default = `0.8`.
#' @param align_alleles Try to align tag alleles to target alleles (if proxies = 1). `1` = yes, `0` = no. The default is `1`.
#' @param palindromes Allow palindromic SNPs (if proxies = 1). `1` = yes, `0` = no. The default is `1`.
#' @param maf_threshold MAF threshold to try to infer palindromic SNPs. The default is `0.3`.
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data.
#' @param splitsize The default is `10000`.
#' @param proxy_splitsize The default is `500`.
#'
#' @export
#' @return Dataframe of summary statistics for all available outcomes
extract_outcome_data <- function(snps, outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token = ieugwasr::check_access_token(), splitsize=10000, proxy_splitsize=500)
{
	# .Deprecated("ieugwasr::associations()")
	outcomes <- ieugwasr::legacy_ids(unique(outcomes))

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



extract_outcome_data_internal <- function(snps, outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token = ieugwasr::check_access_token(), splitsize=10000)
{
	snps <- unique(snps)
	message("Extracting data for ", length(snps), " SNP(s) from ", length(unique(outcomes)), " GWAS(s)")
	outcomes <- unique(outcomes)

	if(proxies == FALSE)
	{
		proxies <- 0
	} else if(proxies == TRUE)
	{
		proxies <- 1
	} else {
		stop("'proxies' argument should be TRUE or FALSE")
	}

	if((length(snps) < splitsize & length(outcomes) < splitsize) | (length(outcomes) < splitsize & length(snps) < splitsize))
	{

		d <- ieugwasr::associations(
			variants = snps, 
			id = outcomes,
			proxies = proxies,
			r2 = rsq,
			align_alleles = align_alleles,
			palindromes = palindromes,
			maf_threshold = maf_threshold,
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
				out <- ieugwasr::associations(
					variants = x$snps, 
					id = outcomes[i],
					proxies = proxies,
					r2 = rsq,
					align_alleles = align_alleles,
					palindromes = palindromes,
					maf_threshold = maf_threshold,
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

				out <- ieugwasr::associations(
					variants = snps[i], 
					id = x$outcomes,
					proxies = proxies,
					r2 = rsq,
					align_alleles = align_alleles,
					palindromes = palindromes,
					maf_threshold = maf_threshold,
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
		d$data_source.outcome <- "igd"
		return(d)
	} else {
		return(NULL)
	}
}



#' Avoid issues in MR by finding impossible vals and setting to NA
#'
#' @param d Data frame
#' @return Cleaned data frame
#' @keywords internal
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
#' @importFrom stats qnorm
#' @export
get_se <- function(eff, pval)
{
	abs(eff) / abs(qnorm(pval / 2))
}


#' Format the returned table from the MySQL database
#'
#' @param d Data frame
#' @return Data frame
#' @keywords internal
format_d <- function(d)
{

	d1 <- data.frame(
		SNP = as.character(d$rsid),
		chr = as.character(d$chr),
		pos = as.character(d$position),
		beta.outcome = as.numeric(d$beta),
		se.outcome = as.numeric(d$se),
		samplesize.outcome = as.numeric(d$n),
		pval.outcome = as.numeric(d$p),
		eaf.outcome = as.numeric(d$eaf),
		effect_allele.outcome = as.character(d$ea),
		other_allele.outcome = as.character(d$nea),
		outcome = as.character(d$trait),
		id.outcome = as.character(d$id),
		stringsAsFactors=FALSE
	)

	if("proxy" %in% names(d))
	{
		p <- data.frame(
			proxy.outcome = d$proxy,
			target_snp.outcome = d$target_snp,
			proxy_snp.outcome = d$proxy_snp,
			target_a1.outcome = d$target_a1,
			target_a2.outcome = d$target_a2,
			proxy_a1.outcome = d$proxy_a1,
			proxy_a2.outcome = d$proxy_a2,
			stringsAsFactors = FALSE
		)
		d <- cbind(d1, p)

		# If two SNPs have the same proxy SNP then one has to be removed
		d <- plyr::ddply(d, c("outcome"), function(x)
		{
			x <- plyr::mutate(x)
			subset(x, !duplicated(proxy_snp.outcome))
		})

	} else {
		d <- d1
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




