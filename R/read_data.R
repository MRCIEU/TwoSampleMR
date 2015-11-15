#' Read in summary statistics for exposure trait
#'
#' @param filename Filename
#' @param exposure Name of exposure trait
#' @param quote Character used for quotes
#' @param sep Character used to delimit columns
#' @export
#' @return Data frame of exposure summary stats
read_exposure_data <- function(filename, exposure, quote='"', sep=" ")
{
	exposure_dat <- read.csv(filename, header=T, stringsAsFactors=FALSE, quote=quote, sep=sep)
	format_exposure_dat(exposure_dat, exposure)

}


#' Format exposure_dat into the right shape
#'
#' @param exposure_dat Data frame
#' @param exposure Name of exposure trait
#' @export
format_exposure_dat <- function(exposure_dat, exposure, pval_only=FALSE)
{
	# Check all the columns are there as expected
	stopifnot(
		all(c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele") %in% names(exposure_dat)) |
		all(c("SNP", "P_value") %in% names(exposure_dat))
	)

	mrcols <- c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele")
	pcols <- c("SNP", "P_value")
	allcols <- c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "P_value")

	if(all(allcols %in% names(exposure_dat)))
	{
		type <- "all"
	} else if(all(mrcols) %in% names(exposure_dat)) {
		type <- "mr"
	} else if(all(pcols) %in% names(exposure_dat)) {
		type <- "pval"
	} else {
		stop(
			"Error: Unexpected columns.\n\n"
			"For MR the data must have: \n",
			paste(mrcols, collapse=" "), "\n\n",
			"For p-value based tests the data must have: \n",
			paste(pcols, collapse=" "), "\n\n",
			"For all tests the data must have: \n",
			paste(allcols, collapse=" "), "\n"
		)
	}

	exposure_dat$SNP <- tolower(exposure_dat$SNP)
	exposure_dat$SNP <- gsub("[[:space:]]", "", exposure_dat$SNP)


	if(type %in% c("all", "mr"))
	{
		# Do some checks
		stopifnot(all(is.numeric(exposure_dat$beta), na.rm=TRUE))
		stopifnot(all(is.numeric(exposure_dat$se), na.rm=TRUE))
		stopifnot(all(is.numeric(exposure_dat$eaf), na.rm=TRUE))
		stopifnot(all(exposure_dat$eaf > 0, na.rm=TRUE))
		stopifnot(all(exposure_dat$eaf < 1, na.rm=TRUE))
		exposure_dat$exposure <- exposure
		exposure_dat$effect_allele <- toupper(exposure_dat$effect_allele)
		exposure_dat$other_allele <- toupper(exposure_dat$other_allele)

		exposure_dat$keep <- apply(exposure_dat, 1, function(x) any(is.na(x)))
		if(any(!exposure_dat$keep))
		{
			message("Warning: The following SNP(s) are missing required information for the MR tests. They will be excluded. Sorry. We did all we could.")
			message("Atenção: O SNP (s) seguinte estão faltando informações necessárias. Eles serão excluídos. Desculpe. Fizemos tudo o que podíamos.")			
			print(subset(exposure_dat, !keep))
		}
	}

	if(type %in% c("all", "pval"))
	{
		stopifnot(all(is.numeric(exposure_dat$P_value), na.rm=TRUE))
		stopifnot(all(exposure_dat$P_value >= 0, na.rm=TRUE))
		stopifnot(all(exposure_dat$P_value < 1, na.rm=TRUE))
	}


	# Check for missing values 

	column_set1 <- c(
		"eaf",
		"effect_allele",
		"other_allele",
		"beta",
		"se"
	)
	column_set2 <- "P_value"

	apply(exposure_dat[,columns_to_check], 1, function(x) )

	i1 <- !is.na(exposure_dat$eaf)
	i2 <- !is.na(exposure_dat$effect_allele)
	i3 <- !is.na(exposure_dat$other_allele)
	i4 <- !is.na(exposure_dat$beta)
	i5 <- !is.na(exposure_dat$se)




	i6 <- !is.na(exposure_dat$P_value)
	
	if(pval_only)
	{
		exposure_dat$keep <- i6
		if(any(!exposure_dat$keep))
		{
			message("Warning: The following SNP(s) are missing required p_value information. They will be excluded. Sorry. We did all we could.")
			print(subset(exposure_dat, !keep))
		}
	}else{
		exposure_dat$keep <- i1 & i2 & i3 & i4 & i5
	}
	exposure_dat <- subset(exposure_dat, keep)

	stopifnot(nrow(exposure_dat) > 0)

	# Get SNP positions
	bm <- ensembl_get_position(exposure_dat$SNP)
	missing <- exposure_dat$SNP[! exposure_dat$SNP %in% bm$refsnp_id]
	if(length(missing) > 0)
	{
		message("Warning: The following SNP(s) were not present in ensembl GRCh37. They will be excluded. Sorry. This is Matt's fault.")
		message("Atenção: O SNP (s) seguinte não estavam presentes no GRCh37 Ensembl. Eles serão excluídos. Desculpe. Isso é culpa do Matt.")
		print(missing)
	}
	stopifnot(nrow(exposure_dat) > 0)

	i6 <- exposure_dat$effect_allele %in% c("A", "C", "T", "G")
	i7 <- exposure_dat$other_allele %in% c("A", "C", "T", "G")


	exposure_dat <- exposure_dat[,names(exposure_dat)!="chr_name"]
	exposure_dat <- exposure_dat[,names(exposure_dat)!="chrom_start"]
	exposure_dat <- merge(exposure_dat, bm, by.x="SNP", by.y="refsnp_id", all.x=TRUE)

	names(exposure_dat)[names(exposure_dat) == "effect_allele"] <- "effect_allele.exposure"
	names(exposure_dat)[names(exposure_dat) == "other_allele"] <- "other_allele.exposure"
	names(exposure_dat)[names(exposure_dat) == "beta"] <- "beta.exposure"
	names(exposure_dat)[names(exposure_dat) == "se"] <- "se.exposure"
	names(exposure_dat)[names(exposure_dat) == "eaf"] <- "eaf.exposure"
	names(exposure_dat)[names(exposure_dat) == "P_value"] <- "pval.exposure"

	return(exposure_dat)
}

ensembl_get_position <- function(snp)
{
	library(biomaRt)
	Mart <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
	Attr<-listAttributes(Mart)
	ensembl<-getBM(attributes=c("refsnp_id","chr_name","chrom_start"),filters="snp_filter",values=snp,mart=Mart)
	ensembl<-ensembl[nchar(ensembl$chr_name)<=2,]
	ensembl
}


#' Read in summary statistics for outcome trait
#'
#' @param filename Filename
#' @param outcome Name of exposure trait
#' @param quote Character used for quotes
#' @param sep Character used to delimit columns
#' @export
#' @return Data frame of exposure summary stats
read_outcome_data <- function(filename, outcome, quote='"', sep=" ")
{
	outcome_dat <- read.csv(filename, header=T, stringsAsFactors=FALSE, quote=quote, sep=sep)

	# Check all the columns are there as expected
	stopifnot(all(c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele") %in% names(outcome_dat)))
	names(outcome_dat)[names(outcome_dat) == "effect_allele"] <- "effect_allele.outcome"
	names(outcome_dat)[names(outcome_dat) == "other_allele"] <- "other_allele.outcome"
	names(outcome_dat)[names(outcome_dat) == "beta"] <- "beta.outcome"
	names(outcome_dat)[names(outcome_dat) == "se"] <- "se.outcome"
	names(outcome_dat)[names(outcome_dat) == "eaf"] <- "eaf.outcome"
	outcome_dat$outcome <- outcome

	return(outcome_dat)
}



newfunction <- function()
{
	print("hello / bom dia!")
}