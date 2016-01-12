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
#' Check that SNP is there
#' If SNP isn't there, error
#' If SNP is there, check that beta, se, eaf, a1, a2 are there
#' For all SNPs with those vals, say MR possible
#' For all without, say that MR is not possible but fisher test is
#' If some p values provided, fill in p-values where possible
#' If any SNPs not found in ensembl then remove
#'
#' @param exposure_dat Data frame
#' @param exposure Name of exposure trait
#' @export
format_exposure_dat <- function(exposure_dat, exposure)
{
	# Check all the columns are there as expected

	mrcols <- c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele")

	if(! "SNP" %in% names(exposure_dat))
	{
		stop(
			"Must have at least 'SNP' column to perform any analysis\n",
			"For MR analysis must have the following columns:\n",
			paste(mrcols, collapse=" "), "\n\n"	
		)
	}

	exposure_dat$SNP <- tolower(exposure_dat$SNP)
	exposure_dat$SNP <- gsub("[[:space:]]", "", exposure_dat$SNP)

	domr <- FALSE
	if(all(mrcols %in% names(exposure_dat)))
	{
		message("All necessary columns present for MR analysis")
		domr <- TRUE
	} else {
		message("MR analysis can't be performed. Must have the following columns:\n", paste(mrcols, collapse=" "), "\n")
	}

	if(domr)
	{
		# Do some checks
		stopifnot(all(is.numeric(exposure_dat$beta), na.rm=TRUE))
		stopifnot(all(is.numeric(exposure_dat$se), na.rm=TRUE))
		stopifnot(all(is.numeric(exposure_dat$eaf), na.rm=TRUE))
		stopifnot(all(exposure_dat$eaf > 0, na.rm=TRUE))
		stopifnot(all(exposure_dat$eaf < 1, na.rm=TRUE))
		exposure_dat$effect_allele <- toupper(exposure_dat$effect_allele)
		exposure_dat$other_allele <- toupper(exposure_dat$other_allele)
		stopifnot(all(exposure_dat$effect_allele %in% c("A", "C", "T", "G")))
		stopifnot(all(exposure_dat$other_allele %in% c("A", "C", "T", "G")))
		exposure_dat$exposure <- exposure

		exposure_dat$mr_keep <- apply(exposure_dat[,mrcols], 1, function(x) !any(is.na(x)))
		if(any(!exposure_dat$mr_keep))
		{
			message("Warning: The following SNP(s) are missing required information for the MR tests. They will be excluded. Sorry. We did all we could.")
			message("Atenção: O SNP (s) seguinte estão faltando informações necessárias. Eles serão excluídos. Desculpe. Fizemos tudo o que podíamos.")			
			print(subset(exposure_dat, !mr_keep))
		}
	}

	if("P_value" %in% names(exposure_dat))
	{
		badp <- !is.numeric(exposure_dat$P_value) | exposure_dat$P_value <= 0 | exposure_dat$P_value >= 1 | !is.finite(exposure_dat$P_value)
		if(any(badp))
		{
			message("P-values have been provided but the following SNPs have issues with the p-values:\n", paste(exposure_dat$SNP[badp], collapse="\n"))
		}
		if(all(c("beta", "se") %in% names(exposure_dat)))
		{
			goodp <- badp & is.finite(exposure_dat$beta) & is.finite(exposure_dat$se)
			exposure_dat$P_value[goodp] <- pnorm(abs(exposure_dat$beta[goodp]) / exposure_dat$se[goodp], lower.tail=FALSE)
		}
	} else if(all(c("beta", "se") %in% names(exposure_dat))) {
		exposure_dat$P_value <- pnorm(abs(exposure_dat$beta) / exposure_dat$se, lower.tail=FALSE)
	} else {
		exposure_dat$P_value <- 0.99
	}

	exposure_dat$P_value[!is.finite(exposure_dat$P_value)] <- 0.99


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



	exposure_dat <- exposure_dat[,names(exposure_dat)!="chr_name", drop=FALSE]
	exposure_dat <- exposure_dat[,names(exposure_dat)!="chrom_start", drop=FALSE]
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




#' Get data selected from GWAS catalog into correct format
#'
#' Subset the GWAS catalogue to have the rows you require for instrumenting a particular exposure and then run this command.
#' Be careful to avoid using different phenotypes, phenotype types, or units together.
#' @param gwas_catalog_subset Subset of rows from \code{data(gwas_catalog)}
#' @param type Are these data used as "exposure" or "outcome"? Default is "exposure"
#' @param  traitname If specified, will name the exposure/outcome this variable. Otherwise (default) will name it based on the Phenotype columnin \code{gwas_catalog_subset}
#' @export
#' @return Data frame
#' @examples \dontrun{
#' data(gwas_catalog)
#' bmi <- subset(gwas_catalog, Phenotype=="Body mass index" & Year==2010 & grepl("kg", Units)
#' bmi <- format_gwas_catalog(bmi)
#'}
format_gwas_catalog <- function(gwas_catalog_subset, type="exposure", traitname=NULL)
{
	stopifnot(type %in% c("exposure", "outcome"))
	if(is.null(traitname))
	{
		ph <- paste(gwas_catalog_subset$Phenotype, gwas_catalog_subset[["Phenotype info"]])
		if(length(unique(ph)) > 1)
		{
			warning("There is more than one Phenotype / Phenotype info combination. This can cause problems in MR. Using the first entry.")
		}
		traitname <- gwas_catalog_subset$Phenotype[1]
	}
	if(length(unique(gwas_catalog_subset$Units)) > 1)
	{
		warning("The effect sizes selected are in different units. This will cause problems in MR")
	}
	gwas_catalog_subset <- subset(gwas_catalog_subset, select=c("SNP", "Effect", "eaf", "Allele", "other_allele", "SE", "P-value"))
	names(gwas_catalog_subset) <- c("SNP", "beta", "eaf", "effect_allele", "other_allele", "se", "P_value")
	format_exposure_dat(gwas_catalog_subset, traitname)
}
