
#' Toggle API address between development and release
#'
#' @param where Which API to use. Choice between "local", "release", "test". Default = "local"
#'
#' @export
#' @return NULL
toggle_dev <- function(where="local")
{
	stopifnot(where %in% c("local", "release", "test"))
	release <- "http://api.mrbase.org/"

	url <- switch(where,
		local = "http://localhost/",
		test = "http://apitest.mrbase.org/",
		release = "http://api.mrbase.org/"
	)

	options(mrbaseapi=url)
	message("Changing API to ", where, ": ", url)
}





#' Get access token for OAuth2 access to MR Base
#'
#'
#' @export
#' @return access token string
get_mrbase_access_token <- function()
{
	tf <- basename(tempfile())
	check <- suppressWarnings(file.create(tf))
	if(!check)
	{
		stop("You are currently in a directory which doesn't have write access.\n",
			"  In order to authenticate we need to store the credentials in a file called '.httr-oauth'.\n",
			"  Please setwd() to a different directory where you have write access.")
	} else {
		unlink(tf)
	}
	a <- googleAuthR::gar_auth()
	if(! a$validate())
	{
		a$refresh()
	}
	return(a$credentials$access_token)
}


#' Revoke access token for MR Base
#'
#' @export
#' @return NULL
revoke_mrbase_access_token <- function()
{
	a <- googleAuthR::gar_auth()
	a$revoke()
}


fromJSON_safe <- function(url)
{
	# r <- readLines(url, warn=FALSE)
	r <- RCurl::getURL(url, timeout=300)
	return(jsonlite::fromJSON(r))
}

#' Check MR Base access level
#'
#' In order to be granted access to a particular dataset that is not generally available please contact the developers.
#'
#'
#' @export
#' @return access level string
check_mrbase_access <- function()
{
	url <- paste0(options()$mrbaseapi, "check_token?access_token=", get_mrbase_access_token())
	d <- fromJSON_safe(url)
	return(d)
}

#' Get number of studies in database
#'
#'
#' @export
#' @return integer
get_mrbase_status <- function()
{
	url <- paste0(options()$mrbaseapi, "get_status")
	d <- fromJSON_safe(url)[1,1]
	return(d)	
}

#' Get list of studies with available GWAS summary statistics through API
#'
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data
#'
#' @export
#' @return Dataframe of details for all available studies
available_outcomes <- function(access_token = get_mrbase_access_token())
{
	url <- paste0(options()$mrbaseapi, "get_studies?access_token=", access_token)
	return(fromJSON_safe(url))
}



#' Upload a file using POST request through API
#'
#' @param x Vector that will be written to a file to be posted
#' @param max_file_size Maximum file size permitted to be uploaded in bytes (16Mb default)
#' @param header Header in the file?
#' @return basename of file
upload_file_to_api <- function(x, max_file_size=16*1024*1024, header=FALSE)
{
	uri <- paste0(options()$mrbaseapi, "upload")
	filename <- paste0(tempfile(), ".txt")
	f <- file(filename, open="wb")
	write.table(x, file=f, row=F, col=header, qu=F, eol="\n")
	close(f)
	if(file.size(filename) > max_file_size) stop("File size is too large, your request has too many SNPs")
	suppressWarnings(RCurl::postForm(uri, file=RCurl::fileUpload(filename=filename)))
	unlink(filename)
	return(basename(filename))
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
#'
#' @export
#' @return Dataframe of summary statistics for all available outcomes
extract_outcome_data <- function(snps, outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token = get_mrbase_access_token())
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

	splitsize <- 50

	if((length(snps) < 5 & length(outcomes) < splitsize) | (length(outcomes) < 5 & length(snps) < splitsize))
	{
		snpfile <- upload_file_to_api(snps)
		outcomefile <- upload_file_to_api(outcomes)

		url <- paste0(options()$mrbaseapi, "get_effects_from_file?",
			"access_token=", access_token,
			"&outcomefile=", outcomefile,
			"&snpfile=", snpfile,
			"&proxies=", proxies,
			"&rsq=", rsq,
			"&align_alleles=", align_alleles,
			"&palindromes=", palindromes,
			"&maf_threshold=", maf_threshold
		)
		d <- fromJSON_safe(url)

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
				snpfile <- upload_file_to_api(x$snps)
				outcomefile <- upload_file_to_api(outcomes[i])

				url <- paste0(options()$mrbaseapi, "get_effects_from_file?",
					"access_token=", access_token,
					"&outcomefile=", outcomefile,
					"&snpfile=", snpfile,
					"&proxies=", proxies,
					"&rsq=", rsq,
					"&align_alleles=", align_alleles,
					"&palindromes=", palindromes,
					"&maf_threshold=", maf_threshold
				)
				out <- fromJSON_safe(url)
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
				snpfile <- upload_file_to_api(snps[i])
				outcomefile <- upload_file_to_api(x$outcomes)

				url <- paste0(options()$mrbaseapi, "get_effects_from_file?",
					"access_token=", access_token,
					"&outcomefile=", outcomefile,
					"&snpfile=", snpfile,
					"&proxies=", proxies,
					"&rsq=", rsq,
					"&align_alleles=", align_alleles,
					"&palindromes=", palindromes,
					"&maf_threshold=", maf_threshold
				)
				out <- fromJSON_safe(url)
				if(!is.data.frame(out)) out <- data.frame()
				return(out)
			})
		}

		d <- plyr::rbind.fill(d)

	}
	if(is.null(nrow(d)) | nrow(d) == 0)
	{
		message("None of the requested SNPs were available in the specified GWASs.")
		return(NULL)
	}
	d <- format_d(d)
	d$data_source.outcome <- "mrbase"
	return(d)
}



extract_outcome_data_simple <- function(snps, outcomes, proxies = 0, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token = get_mrbase_access_token())
{
	snps <- unique(snps)
	message("Extracting data for ", length(snps), " SNP(s) from ", length(unique(outcomes)), " GWAS(s)")
	outcomes <- unique(outcomes)

	snpfile <- upload_file_to_api(snps)
	outcomefile <- upload_file_to_api(outcomes)

	url <- paste0(options()$mrbaseapi, "get_effects_from_file?",
		"access_token=", access_token, 
		"&outcomefile=", outcomefile, 
		"&snpfile=", snpfile,
		"&proxies=", proxies,
		"&rsq=", rsq,
		"&align_alleles=", align_alleles,
		"&palindromes=", palindromes,
		"&maf_threshold=", maf_threshold
	)
	d <- fromJSON_safe(url)

	if(length(d) == 0)
	{
		message("None of the requested SNPs were available in the specified GWASs.")
		return(NULL)
	}
	d <- format_d(d)
	d$data_source.outcome <- "mrbase"
	return(d)
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
	d1 <- data.frame(
		SNP = as.character(d$name),
		beta.outcome = as.numeric(d$beta),
		se.outcome = as.numeric(d$se),
		samplesize.outcome = as.numeric(d$n),
		ncase.outcome = as.numeric(d$ncase),
		ncontrol.outcome = as.numeric(d$ncontrol),
		pval.outcome = as.numeric(d$p),
		eaf.outcome = as.numeric(d$effect_allelel_freq),
		effect_allele.outcome = as.character(d$effect_allele),
		other_allele.outcome = as.character(d$other_allele),
		units.outcome = as.character(d$unit),
		outcome = as.character(d$trait),
		consortium.outcome = as.character(d$consortium),
		year.outcome = as.numeric(d$year),
		pmid.outcome = as.numeric(d$pmid),
		id.outcome = as.character(d$id),
		category.outcome = as.character(d$category),
		subcategory.outcome = as.character(d$subcategory),
		originalname.outcome = 0,
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
	d$outcome <- paste0(d$outcome, " || ", d$consortium.outcome, " || ", d$year.outcome)

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




# Extract outcome data given a set of SNPs
#'
#' Supply the output from \code{read_exposure_data} and all the SNPs therein will be queried against the requested outcomes in remote database using API.
#' WARNING: This is unlikely to work correctly if there are a large number of SNPs
#'
#' @param exposure_dat Output from \code{read_exposure_data}
#' @param outcomes Array of IDs (see \code{id} column in output from \code{available_outcomes})
#' @export
#' @return Dataframe of summary statistics for all available outcomes
extract_outcome_data_using_get <- function(snps, outcomes)
{
	access_token <- get_mrbase_access_token()
	snps <- unique(snps)
	message("Extracting data for ", length(snps), " SNP(s) from ", length(unique(outcomes)), " GWAS(s)")
	outcomes <- unique(outcomes)
	snps <- paste(snps, collapse=",")
	outcomes <- paste(outcomes, collapse=",")

	url <- paste0(options()$mrbaseapi, "get_effects?access_token=", access_token, "&outcomes=", outcomes, "&snps=", snps)
	d <- fromJSON_safe(url)
	if(length(d) == 0)
	{
	  message("None of the requested SNPs were available in the specified GWASs.")
	  return(NULL)
	}
	return(format_d(d))
}
