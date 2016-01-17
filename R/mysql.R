
#' Get list of studies with available GWAS summary statistics through API
#'
#' @param password If correct password is supplied then access to restricted studies will be supplied
#'
#' @export
#' @return Dataframe of details for all available studies
available_outcomes <- function(password=NULL)
{
	password <- get_password(password)
	url <- paste0("http://scmv-webapps.epi.bris.ac.uk:5000/get_studies?password=", password)
	d <- fromJSON(url)
	return(d)
}


#' Upload a file using POST request through API
#'
#' @param x Vector that will be written to a file to be posted
#' @param max_file_size Maximum file size permitted to be uploaded in bytes (16Mb default)
#' @return
upload_file_to_api <- function(x, max_file_size=16*1024*1024, header=FALSE)
{
	require(RCurl)
	uri <- "http://scmv-webapps.epi.bris.ac.uk:5000/upload"
	filename <- paste0(tempfile(), ".txt")
	write.table(x, filename, row=F, col=header, qu=F)
	if(file.size(filename) > max_file_size) stop("File size is too large, your request has too many SNPs")
	suppressWarnings(postForm(uri, file=fileUpload(filename=filename)))
	unlink(filename)
	return(basename(filename))
}



# Extract summary statistics from MySQL db through API given a list of SNPs and outcomes
#'
#' Supply the output from \code{read_exposure_data} and all the SNPs therein will be queried against the requested outcomes in remote database using API.
#'
#' @param exposure_dat Output from \code{read_exposure_data}
#' @param outcomes Array of IDs (see \code{id} column in output from \code{available_outcomes})
#' @param password If correct password is supplied then access to restricted studies will be supplied
#' @export
#' @return Dataframe of summary statistics for all available outcomes
extract_outcome_data <- function(exposure_dat, outcomes, password=NULL)
{
	password <- get_password(password)
	message("Extracting data for ", nrow(exposure_dat), " SNP(s) from ", length(unique(outcomes)), " GWAS(s)")
	snps <- unique(exposure_dat$SNP)
	outcomes <- unique(outcomes)

	snpfile <- upload_file_to_api(snps)
	outcomefile <- upload_file_to_api(outcomes)

	url <- paste0("http://scmv-webapps.epi.bris.ac.uk:5000/get_effects_from_file?password=", password, "&outcomefile=", outcomefile, "&snpfile=", snpfile)
	d <- fromJSON(url)
	return(format_d(d))
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
	d <- data.frame(
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
		outcome = as.character(d$trait),
		consortium.outcome = as.character(d$consortium),
		year.outcome = as.numeric(d$year),
		pmid.outcome = as.numeric(d$pmid),
		id.outcome = as.numeric(d$id),
		displayname.outcome = 0
	)

	if(nrow(d) == 0)
	{
		message("No matches")
		return(d)
	}

	d$displayname.outcome <- paste0(d$outcome, " (", d$consortium.outcome, " ", d$year.outcome, ")")

	rem <- is.na(d$beta.outcome) & is.na(d$pval.outcome)
	d <- subset(d, !rem)

	# For any that have missing SE but available beta and pval, infer SE
	index <- (is.na(d$se.outcome) | d$se.outcome == 0) & (!is.na(d$beta.outcome) & !is.na(d$pval.outcome))
	if(any(index))
	{
		d$se.outcome[index] <- get_se(d$beta.outcome[index], d$pval.outcome[index])
	}

	d <- cleanup_outcome_data(d)
	return(d)	
}


##
## Deprecated
##


#' Get list of studies with available GWAS summary statistics
#'
#' @param user User name
#' @param password Password
#' @param dbname Database
#' @param host Host
#' @return Dataframe of study details
available_outcomes_mysql <- function(user, password, dbname, host, port)
{
	mydb <- dbConnect(MySQL(), user=user, password=password, dbname=dbname, host=host)
	query <- "SELECT * from study;"
	out <- dbSendQuery(mydb, query)
	d <- fetch(out, n=-1)
	
	dbClearResult(dbListResults(mydb)[[1]])
	dbDisconnect(mydb)
	return(d)
}

#' Extract SNP effects from GWAS summary statistics
#'
#' Any SNPs that have no beta or p value will be removed from the outcome results
#' 
#' @param exposure_dat Output from \code{read_exposure_data}
#' @param outcomes List of study names to search for, obtained from \code{available_outcomes}
#' @param user User name
#' @param password Password
#' @param dbname Database
#' @param host Host
#' @return Dataframe of summary statistics for all available outcomes
extract_outcome_data_mysql <- function(exposure_dat, outcomes, user, password, dbname, host, port)
{
	message("Extracting data for ", nrow(exposure_dat), " SNPs")
	snps <- paste(exposure_dat$SNP, collapse="', '")
	outcomes <- paste(outcomes, collapse="', '")

	message("Connecting to database")
	mydb <- dbConnect(MySQL(), user=user, password=password, dbname=dbname, host=host)

    message("Extracting data")
	query <- paste(
		"SELECT a.effect_allele, a.other_allele, a.effect_allelel_freq, a.beta, a.se, a.p, a.n, b.name, c.* ",
		"FROM assoc a, snp b, study c ",
		"WHERE a.snp=b.id AND a.study=c.id ",
		"AND a.study IN ('", outcomes, "') ",
		"AND b.name IN ('", snps, "') ",
		"ORDER BY a.study;", sep="")

	out <- dbSendQuery(mydb, query)
	d <- fetch(out, n=-1)
	dbClearResult(dbListResults(mydb)[[1]])
	dbDisconnect(mydb)

	return(format_d(d))
}



# Extract outcome data given a set of SNPs
#'
#' Supply the output from \code{read_exposure_data} and all the SNPs therein will be queried against the requested outcomes in remote database using API.
#' WARNING: This is unlikely to work correctly if there are a large number of SNPs
#'
#' @param exposure_dat Output from \code{read_exposure_data}
#' @param outcomes Array of IDs (see \code{id} column in output from \code{available_outcomes})
#' @param password If correct password is supplied then access to restricted studies will be supplied
#' @return Dataframe of summary statistics for all available outcomes
extract_outcome_data_using_get <- function(exposure_dat, outcomes, password=NULL)
{
	password <- get_password(password)
	message("Extracting data for ", nrow(exposure_dat), " SNPs")
	snps <- paste(exposure_dat$SNP, collapse=",")
	outcomes <- paste(outcomes, collapse=",")

	url <- paste0("http://scmv-webapps.epi.bris.ac.uk:5000/get_effects?password=", password, "&outcomes=", outcomes, "&snps=", snps)
	d <- fromJSON(url)

	return(format_d(d))
}
