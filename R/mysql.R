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
#' Any SNPs that have no beta or p value will be removed from the outcome results
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
		"SELECT a.effect_allele, a.other_allele, a.effect_allelel_freq, a.beta, a.se, a.p, a.n, b.name, c.* ",
		"FROM assoc a, snp b, study c ",
		"WHERE a.snp=b.id AND a.study=c.id ",
		"AND a.study IN ('", outcomes, "') ",
		"AND b.name IN ('", snps, "') ",
		"ORDER BY a.study;", sep="")

	out <- dbSendQuery(mydb, query)
	d <- fetch(out, n=-1)

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
		id.outcome = as.numeric(d$id)
	)

	d$displayname.outcome <- paste0(d$outcome, "; ", d$consortium.outcome, " (", d$year.outcome, ")")

	d$eaf.outcome[is.na(d$eaf.outcome)] <- 0.5

	rem <- is.na(d$beta.outcome) & is.na(d$pval.outcome)
	d <- subset(d, !rem)

	# For any that have missing SE but available beta and pval, infer SE
	# index <- (is.na(d$se.outcome) | d$se.outcome == 0) & (!is.na(d$beta.outcome) & !is.na(d$pval.outcome))
	# if(any(index))
	# {
	# 	message("Inferring SE")
	# 	d$se.outcome[index] <- get_se(d$beta.outcome[index], d$pval.outcome[index])
	# }

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
