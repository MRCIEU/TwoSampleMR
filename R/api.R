

#' Toggle API address between development and release
#'
#' @param where Which API to use. Choice between "local", "release", "test". Default = "local"
#'
#' @export
#' @return NULL
toggle_api <- function(where="test")
{
	url <- switch(where,
		test = "http://ieu-db-interface.epi.bris.ac.uk:8084/",
		dev = "http://localhost:8019/"
	)
	if(is.null(url))
	{
		url <- options()$mrbaseapi
		warning("A valid API was not selected. No change")
	}

	options(mrbaseapi=url)
	message("API: ", where, ": ", url)
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
	a <- googleAuthR::gar_auth("mrbase.oauth")
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
	a <- googleAuthR::gar_auth("mrbase.oauth")
	a$revoke()
}


#' Wrapper for sending queries and payloads to API
#'
#' There are a number of different GET and POST endpoints in the GWAS database API. This is a generic way to access them
#'
#' @param path Either a full query path (e.g. for get) or an endpoint (e.g. for post) queries
#' @param query If post query, provide a list of arguments as the payload. NULL by default
#' @param access_token=get_mrbase_access_token()
#'
#' @export
#' @return Parsed json output from query, often in form of data frame
api_query <- function(path, query=NULL, access_token=get_mrbase_access_token())
{
	ntry <- 0
	ntries <- 3
	headers <- httr::add_headers(
		# 'Content-Type'='application/json; charset=UTF-8',
		'X-Api-Token'=access_token,
		'X-Api-Source'=ifelse(is.null(options()$mrbase.environment), 'R/TwoSampleMR', 'mr-base-shiny')
	)

	while(ntry <= ntries)
	{
		if(!is.null(query))
		{
			r <- try(
				httr::POST(
					paste0(options()$mrbaseapi, path),
					body = query, 
					headers,
					encode="json",
					httr::timeout(300)
				),
				silent=TRUE
			)
		} else {
			r <- try(
				httr::GET(
					paste0(options()$mrbaseapi, path),
					headers,
					httr::timeout(300)
				),
				silent=TRUE
			)			
		}
		if(class(r) == 'try-error')
		{
			if(grepl("Timeout", as.character(attributes(r)$condition)))
			{
				stop("The query to MR-Base exceeded 300 seconds and timed out. Please simplify the query")
			}
		}
		if(class(r) != 'try-error')
		{
			break
		}
		ntry <- ntry + 1
	}
	if(class(r) == 'try-error')
	{
		if(grepl("Could not resolve host", as.character(attributes(r)$condition)))
		{
			stop("The MR-Base server appears to be down, the following error was received:\n", as.character(attributes(r)$condition))
		} else {
			stop("The following error was encountered in trying to query the MR-Base server:\n",
				as.character(attributes(r)$condition)
			)
		}
	}

	if(httr::status_code(r) >= 200 & httr::status_code(r) < 300)
	# if(httr::status_code(r) >= 200)
	{
		return(jsonlite::fromJSON(httr::content(r, "text", encoding='UTF-8')))
	} else {
		return(r)
		stop("error code: ", httr::status_code(r), "\n  message: ", jsonlite::fromJSON(httr::content(r, "text", encoding='UTF-8')))
	}
}

# error_codes <- function(code)
# {
# 	codes <- list(
# 		data_frame(code=400, message = "Incorrect"),
# 		data_frame(code=400, message = ""),
# 	)
# }


#' MR-Base server status
#'
#' @export
#' @return list of values regarding status
api_status <- function()
{
	o <- api_query('status')
	class(o) <- "ApiStatus"
	return(o)
}

print.ApiStatus <- function(x)
{
	lapply(names(x), function(y) cat(format(paste0(y, ":"), width=30, justify="right"), x[[y]], "\n"))
}

#' Get list of studies with available GWAS summary statistics through API
#'
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data
#'
#' @export
#' @return Dataframe of details for all available studies
available_outcomes <- function(access_token = get_mrbase_access_token())
{
	message("DEPRECATED. Use gwasinfo() instead.")
	return(api_query('gwasinfo', access_token=access_token))
}


#' Get list of studies with available GWAS summary statistics through API
#'
#' @param id List of MR-Base IDs to retrieve. If NULL (default) retrieves all available datasets
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data
#'
#' @export
#' @return Dataframe of details for all available studies
gwasinfo <- function(id=NULL, access_token = get_mrbase_access_token())
{
	if(!is.null(id))
	{
		stopifnot(is.vector(id))
		out <- api_query('gwasinfo', query = list(id=id), access_token=access_token)
	} else {
		out <- api_query('gwasinfo', access_token=access_token)
	}
	out <- dplyr::bind_rows(out) %>%
		dplyr::select(id, trait, sample_size, nsnp, year, consortium, author, dplyr::everything())
	class(out) <- c("GwasInfo", class(out))
	return(out)
}

# print.GwasInfo <- function(x)
# {
# 	dplyr::glimpse(x)
# }


#' Upload a file using POST request through API
#'
#' @param x Vector that will be written to a file to be posted
#' @param max_file_size Maximum file size permitted to be uploaded in bytes (16Mb default)
#' @param header Header in the file?
#' @return basename of file
upload_file_to_api <- function(x, max_file_size=16*1024*1024, header=FALSE)
{
	message("DEPRECATED. This may no longer work. Please check the manual")
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
		d <- api_query('associations/',
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
				out <- api_query('associations/',
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
				out <- api_query('associations/',
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

