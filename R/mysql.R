#' Get list of studies with available GWAS summary statistics
#'
#' @param user="mruser" User name
#' @param password="TMG_F1WnTL" Password
#' @param dbname="mrbase" Database
#' @param host="epi-franklin.epi.bris.ac.uk" Host
#' @export
#' @return Dataframe of study details
available_outcomes <- function(user="mruser", password="TMG_F1WnTL", dbname="mrbase", host="127.0.0.1", port="3306")
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
#' 
#' @param exposure_dat Output from \code{read_exposure_data}
#' @param outcomes List of study names to search for, obtained from \code{available_outcomes}
#' @param user="mruser" User name
#' @param password="TMG_F1WnTL" Password
#' @param dbname="mrbase" Database
#' @param host="epi-franklin.epi.bris.ac.uk" Host
#' @export
#' @return Dataframe of summary statistics for 
extract_outcome_data <- function(exposure_dat, outcomes, user="mruser", password="TMG_F1WnTL", dbname="mrbase", host="127.0.0.1", port="3306")
{
	message("Extracting data for ", nrow(exposure_dat), " SNPs")
	snps <- paste(exposure_dat$SNP, collapse="', '")
	outcomes <- paste(outcomes, collapse="', '")

	message("Connecting to database")
	mydb <- dbConnect(MySQL(), user=user, password=password, dbname=dbname, host=host)

    message("Extracting data")
	query <- paste(
		"SELECT a.*, b.*, c.* ",
		"FROM assoc a, snp b, study c ",
		"WHERE a.snp=b.id AND a.study=c.id ",
		"AND a.study IN ('", outcomes, "') ",
		"AND b.name IN ('", snps, "') ",
		"ORDER BY a.study;", sep="")

	out <- dbSendQuery(mydb, query)
	d <- fetch(out, n=-1)
	d <- subset(d, select=c(
		name,
		beta,
		se,
		n,
		p,
		effect_allelel_freq,
		effect_allele,
		other_allele,
		trait
	))
	names(d) <- c("SNP", "beta.outcome", "se.outcome", "samplesize.outcome", "pval.outcome", "eaf.outcome", "effect_allele.outcome", "other_allele.outcome", "outcome")

	rem <- is.na(d$beta.outcome) & is.na(d$pval.outcome)
	d <- subset(d, !rem)

	# For any that have missing SE but available beta and pval, infer SE
	index <- (is.na(d$se.outcome) | d$se.outcome == 0) & (!is.na(d$beta.outcome) & !is.na(d$pval.outcome))
	if(any(index))
	{
		message("Inferring SE")
		d$se.outcome[index] <- get_se(d$beta.outcome[index], d$pval.outcome[index])
	}


	message("Formatting")
	d$SNP <- as.character(d$SNP)
	d$beta.outcome <- as.numeric(d$beta.outcome)
	d$se.outcome <- as.numeric(d$se.outcome)
	d$samplesize.outcome <- as.numeric(d$samplesize.outcome)
	d$pval.outcome <- as.numeric(d$pval.outcome)
	d$eaf.outcome <- as.numeric(d$eaf.outcome)
	d$effect_allele.outcome <- as.character(d$effect_allele.outcome)
	d$other_allele.outcome <- as.character(d$other_allele.outcome)
	d$outcome <- as.character(d$outcome)

	d <- cleanup_outcome_data(d)

	dbClearResult(dbListResults(mydb)[[1]])
	dbDisconnect(mydb)
	return(d)
}



cleanup_outcome_data <- function(d)
{
	d$se.outcome[d$se.outcome <= 0] <- NA
	d$eaf.outcome[d$eaf.outcome <= 0 | d$eaf.outcome >= 1] <- NA
	d$beta.outcome[d$beta.outcome == -9] <- NA
	return(d)
}


get_se <- function(eff, pval)
{
	abs(eff) / abs(qnorm(pval / 2))
}
