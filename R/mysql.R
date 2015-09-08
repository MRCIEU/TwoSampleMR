#' Get list of studies with available GWAS summary statistics
#'
#' @param user="mruser" User name
#' @param password="TMG_F1WnTL" Password
#' @param dbname="mrbase" Database
#' @param host="epi-franklin.epi.bris.ac.uk" Host
#' @export
#' @return Dataframe of study details
available_outcomes <- function(user="mruser", password="TMG_F1WnTL", dbname="mrbase", host="epi-franklin.epi.bris.ac.uk")
{
	mydb <- dbConnect(MySQL(), user=user, password=password, dbname=dbname, host=host)
	query <- "SELECT * from study;"
	out <- dbSendQuery(mydb, query)
	d <- fetch(out, n=-1)
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
extract_outcome_data <- function(exposure_dat, outcomes, user="mruser", password="TMG_F1WnTL", dbname="mrbase", host="epi-franklin.epi.bris.ac.uk")
{
	snps <- paste(exposure_dat$SNP, collapse="', '")
	outcomes <- paste(outcomes, collapse="', '")

	mydb <- dbConnect(MySQL(), user=user, password=password, dbname=dbname, host=host)

	query <- paste(
		"SELECT a.*, b.*, c.*",
		"FROM assoc a, snps b, study c",
		"WHERE a.snp=b.id AND a.study=c.id",
		"AND b.name IN ('", snps, "')",
		"AND c.filename IN ('", outcomes, "')",
		"ORDER BY filename;")

	out <- dbSendQuery(mydb, query)
	d <- fetch(out, n=-1)

	d <- subset(d, select=c(name, beta, se, n, p, freq, effect, other, filename,))
	names(d) <- c("SNP", "beta.outcome", "se.outcome", "samplesize.outcome", "pval.outcome", "eaf.outcome", "effect_allele.outcome", "other_allele.outcome", "outcome")

	index <- d$se.outcome==0
	d$se[index] <- d$beta.outcome[index] / pt(d$pval.outcome[index] / 2, df = d$n[index], low=FALSE)
	return(outcome_dat)
}

