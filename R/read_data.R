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

	# Check all the columns are there as expected
	stopifnot(all(c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele") %in% names(exposure_dat)))

	# Do some checks
	stopifnot(all(is.numeric(exposure_dat$beta)))
	stopifnot(all(is.numeric(exposure_dat$se)))
	stopifnot(all(is.numeric(exposure_dat$eaf)))
	stopifnot(all(exposure_dat$eaf > 0))
	stopifnot(all(exposure_dat$eaf < 1))
	exposure_dat$exposure <- exposure

	exposure_dat$keep <- TRUE

	exposure_dat$effect_allele <- toupper(exposure_dat$effect_allele)
	exposure_dat$other_allele <- toupper(exposure_dat$other_allele)
	exposure_dat$SNP <- tolower(exposure_dat$SNP)
	exposure_dat$SNP <- gsub("[[:space:]]", "", exposure_dat$SNP)

	# Check for missing values 
	i1 <- !is.na(exposure_dat$eaf)
	i2 <- !is.na(exposure_dat$effect_allele)
	i3 <- !is.na(exposure_dat$other_allele)
	i4 <- !is.na(exposure_dat$beta)
	i5 <- !is.na(exposure_dat$se)
	exposure_dat$keep <- i1 & i2 & i3 & i4 & i5
	if(any(!exposure_dat$keep))
	{
		message("Warning: The following SNP(s) are missing required information. They will be excluded. Sorry. We did all we could.")
		message("Atenção: O SNP (s) seguinte estão faltando informações necessárias. Eles serão excluídos. Desculpe. Fizemos tudo o que podíamos.")
		
		print(subset(exposure_dat, !keep))
	}
	exposure_dat <- subset(exposure_dat, keep)

	# Get SNP positions
	bm <- ensembl_get_position(exposure_dat$SNP)
	missing <- exposure_dat$SNP[! exposure_dat$SNP %in% bm$refsnp_id]
	if(length(missing) > 0)
	{
		message("Warning: The following SNP(s) were not present in ensembl GRCh37. They will be excluded. Sorry. This it's Matt's fault.")
		message("Atenção: O SNP (s) seguinte não estavam presentes no GRCh37 Ensembl. Eles serão excluídos. Desculpe. Isso é culpa do Matt.")
		print(missing)
	}

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
