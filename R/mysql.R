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
	require(RMySQL)
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
	require(RMySQL)
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


