

#' Read outcome data
#'
#' Reads in outcome data. Checks and organises columns for use with MR or enrichment tests. Infers p-values when possible from beta and se.
#'
#' @param filename Filename. Must have header with at least SNP column present.
#' @param snps=NULL SNPs to extract. If NULL then doesn't extract any and keeps all.
#' @param sep=" " Specify delimeter in file
#' @param phenotype_col="Phenotype" Optional column name for the column with phenotype name corresponding the the SNP. If not present then will be created with the value "Outcome"
#' @param snp_col="SNP" Required name of column with SNP rs IDs
#' @param beta_col="beta" Required for MR. Name of column with effect sizes
#' @param se_col="se" Required for MR. Name of column with standard errors
#' @param eaf_col="eaf" Required for MR. Name of column with effect allele frequency
#' @param effect_allele_col="effect_allele" Required for MR. Name of column with effect allele. Must be "A", "C", "T" or "G"
#' @param other_allele_col="other_allele" Required for MR. Name of column with non effect allele. Must be "A", "C", "T" or "G"
#' @param pval_col="pval" Required for enrichment tests. Name of column with p-value.
#' @param ncase_col="ncase" Optional column name for number of cases.
#' @param ncontrol_col="ncontrol" Optional column name for number of controls.
#' @param samplesize_col="samplesize" Optional column name for sample size.
#' @param gene_col="gene" Optional column name for gene name.
#' @param min_pval=1e-100 Minimum allowed pval
#'
#' @export
#' @return data frame
read_outcome_data <- function(filename, snps=NULL, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", min_pval=1e-100)
{
	outcome_dat <- fread(filename, header=TRUE, sep=sep)
	outcome_dat <- format_data(
		as.data.frame(outcome_dat),
		type="outcome",
		snps=snps,
		sep=sep,
		phenotype_col=phenotype_col,
		snp_col=snp_col,
		beta_col=beta_col,
		se_col=se_col,
		eaf_col=eaf_col,
		effect_allele_col=effect_allele_col,
		other_allele_col=other_allele_col,
		pval_col=pval_col,
		ncase_col=ncase_col,
		ncontrol_col=ncontrol_col,
		samplesize_col=samplesize_col,
		gene_col=gene_col,
		min_pval=min_pval
	)
	return(outcome_dat)
}

#' Read exposure data
#'
#' Reads in exposure data. Checks and organises columns for use with MR or enrichment tests. Infers p-values when possible from beta and se. Looks up SNPs in biomaRt to get basic info.
#'
#' @param filename Filename. Must have header with at least SNP column present.
#' @param sep=" " Specify delimeter in file
#' @param phenotype_col="Phenotype" Optional column name for the column with phenotype name corresponding the the SNP. If not present then will be created with the value "Outcome"
#' @param snp_col="SNP" Required name of column with SNP rs IDs
#' @param beta_col="beta" Required for MR. Name of column with effect sizes
#' @param se_col="se" Required for MR. Name of column with standard errors
#' @param eaf_col="eaf" Required for MR. Name of column with effect allele frequency
#' @param effect_allele_col="effect_allele" Required for MR. Name of column with effect allele. Must be "A", "C", "T" or "G"
#' @param other_allele_col="other_allele" Required for MR. Name of column with non effect allele. Must be "A", "C", "T" or "G"
#' @param pval_col="pval" Required for enrichment tests. Name of column with p-value.
#' @param ncase_col="ncase" Optional column name for number of cases.
#' @param ncontrol_col="ncontrol" Optional column name for number of controls.
#' @param samplesize_col="samplesize" Optional column name for sample size.
#' @param gene_col="gene" Optional column name for gene name.
#' @param min_pval=1e-100 Minimum allowed pval
#'
#' @export
#' @return data frame
read_exposure_data <- function(filename, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", min_pval=1e-100)
{
	exposure_dat <- fread(filename, header=TRUE, sep=sep)
	exposure_dat <- format_data(
		as.data.frame(exposure_dat),
		type="exposure",
		snps=NULL,
		sep=sep,
		phenotype_col=phenotype_col,
		snp_col=snp_col,
		beta_col=beta_col,
		se_col=se_col,
		eaf_col=eaf_col,
		effect_allele_col=effect_allele_col,
		other_allele_col=other_allele_col,
		pval_col=pval_col,
		ncase_col=ncase_col,
		ncontrol_col=ncontrol_col,
		samplesize_col=samplesize_col,
		gene_col=gene_col,
		min_pval=min_pval
	)
	return(exposure_dat)
}

#' Read exposure data
#'
#' Reads in exposure data. Checks and organises columns for use with MR or enrichment tests. Infers p-values when possible from beta and se. If it is the exposure then looks up SNPs in biomaRt to get basic info.
#'
#' @param dat Data frame. Must have header with at least SNP column present.
#' @param type="exposure". Is this the exposure or the outcome data that is being read in?
#' @param snps=NULL SNPs to extract. If NULL then doesn't extract any and keeps all.
#' @param sep=" " Specify delimeter in file
#' @param phenotype_col="Phenotype" Optional column name for the column with phenotype name corresponding the the SNP. If not present then will be created with the value "Outcome"
#' @param snp_col="SNP" Required name of column with SNP rs IDs
#' @param beta_col="beta" Required for MR. Name of column with effect sizes
#' @param se_col="se" Required for MR. Name of column with standard errors
#' @param eaf_col="eaf" Required for MR. Name of column with effect allele frequency
#' @param effect_allele_col="effect_allele" Required for MR. Name of column with effect allele. Must be "A", "C", "T" or "G"
#' @param other_allele_col="other_allele" Required for MR. Name of column with non effect allele. Must be "A", "C", "T" or "G"
#' @param pval_col="pval" Required for enrichment tests. Name of column with p-value.
#' @param ncase_col="ncase" Optional column name for number of cases.
#' @param ncontrol_col="ncontrol" Optional column name for number of controls.
#' @param samplesize_col="samplesize" Optional column name for sample size.
#' @param gene_col="gene" Optional column name for gene name.
#' @param min_pval=1e-100 Minimum allowed pval
#'
#' @export
#' @return data frame
format_data <- function(dat, type="exposure", snps=NULL, sep=" ", header=TRUE, phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", min_pval=1e-100)
{
	all_cols <- c(phenotype_col, snp_col, beta_col, se_col, eaf_col, effect_allele_col, other_allele_col, pval_col, ncase_col, ncontrol_col, samplesize_col)

	i <- names(dat) %in% all_cols
	if(sum(i) == 0)
	{
		stop("None of the specified columns present")
	}
	dat <- dat[,i]

	if(! snp_col %in% names(dat))
	{
		stop("SNP column not found")
	}
	names(dat)[names(dat) == snp_col] <- "SNP"
	dat$SNP <- tolower(dat$SNP)
	dat$SNP <- gsub("[[:space:]]", "", dat$SNP)
	dat <- subset(dat, !is.na(SNP))

	if(!is.null(snps))
	{
		dat <- subset(dat, SNP %in% snps)
	}
	
	if(! phenotype_col %in% names(dat))
	{
		message("No phenotype name specified, defaulting to '", type, "'.")
		dat[[type]] <- type
	} else {
		dat[[type]] <- dat[[phenotype_col]]
		if(phenotype_col != type)
		{
			dat <- subset(dat, select=-c(phenotype_col))
		}
	}

	# Remove duplicated SNPs
	dat <- ddply(dat, type, function(x){
		x <- mutate(x)
		dup <- duplicated(x$SNP)
		if(any(dup))
		{
			warning("Duplicated SNPs present in exposure data for phenotype '", x[[type]][1], ". Just keeping the first instance of the following:\n", paste(x$SNP[dup], collapse="\n"))
			x <- x[!dup,]
		}
		return(x)		
	})

	# Check if columns required for MR are present
	mr_cols <- c(snp_col, beta_col, se_col, eaf_col, effect_allele_col, other_allele_col)
	if(! all(mr_cols %in% names(dat)))
	{
		warning("The following columns are not present and are required for MR analysis\n", paste(mr_cols[!mr_cols %in% names(dat)]), collapse="\n")
		dat$mr_keep.outcome <- FALSE
	} else {
		dat$mr_keep.outcome <- TRUE
	}

	# Check beta
	i <- which(names(dat) == beta_col)[1]
	if(!is.na(i))
	{
		names(dat)[i] <- "beta.outcome"
		if(!is.numeric(dat$beta.outcome))
		{
			warning("beta column is not numeric. Coercing...")
			dat$beta.outcome <- as.numeric(dat$beta.outcome)
		}
		index <- !is.finite(dat$beta.outcome)
		index[is.na(index)] <- TRUE
		dat$beta.outcome[index] <- NA
	}

	# Check se
	i <- which(names(dat) == se_col)[1]
	if(!is.na(i))
	{
		names(dat)[i] <- "se.outcome"
		if(!is.numeric(dat$se.outcome))
		{
			warning("se column is not numeric. Coercing...")
			dat$se.outcome <- as.numeric(dat$se.outcome)
		}
		index <- !is.finite(dat$se.outcome) | dat$se.outcome <= 0
		index[is.na(index)] <- TRUE
		dat$se.outcome[index] <- NA
	}

	# Check eaf
	i <- which(names(dat) == eaf_col)[1]
	if(!is.na(i))
	{
		names(dat)[i] <- "eaf.outcome"
		if(!is.numeric(dat$eaf.outcome))
		{
			warning("eaf column is not numeric. Coercing...")
			dat$eaf.outcome <- as.numeric(dat$eaf.outcome)
		}
		index <- !is.finite(dat$eaf.outcome) | dat$eaf.outcome <= 0 | dat$eaf.outcome >= 1
		index[is.na(index)] <- TRUE
		dat$eaf.outcome[index] <- NA
	}

	# Check effect_allele
	i <- which(names(dat) == effect_allele_col)[1]
	if(!is.na(i))
	{
		names(dat)[i] <- "effect_allele.outcome"
		if(!is.character(dat$effect_allele.outcome))
		{
			warning("effect_allele column is not character data. Coercing...")
			dat$effect_allele.outcome <- as.character(dat$effect_allele.outcome)
		}

		dat$effect_allele.outcome <- toupper(dat$effect_allele.outcome)
		index <- ! dat$effect_allele.outcome %in% c("A", "C", "T", "G")
		index[is.na(index)] <- TRUE
		dat$effect_allele.outcome[index] <- NA
	}


	# Check other_allele
	i <- which(names(dat) == other_allele_col)[1]
	if(!is.na(i))
	{
		names(dat)[i] <- "other_allele.outcome"
		if(!is.character(dat$other_allele.outcome))
		{
			warning("other_allele column is not character data. Coercing...")
			dat$other_allele.outcome <- as.character(dat$other_allele.outcome)
		}

		dat$other_allele.outcome <- toupper(dat$other_allele.outcome)
		index <- ! dat$other_allele.outcome %in% c("A", "C", "T", "G")
		index[is.na(index)] <- TRUE
		dat$other_allele.outcome[index] <- NA
	}


	# Check pval
	i <- which(names(dat) == pval_col)[1]
	if(!is.na(i))
	{
		names(dat)[i] <- "pval.outcome"
		if(!is.numeric(dat$pval.outcome))
		{
			warning("pval column is not numeric. Coercing...")
			dat$pval.outcome <- as.numeric(dat$pval.outcome)
		}
		index <- !is.finite(dat$pval.outcome) | dat$pval.outcome < 0 | dat$pval.outcome > 1
		index[is.na(index)] <- TRUE
		dat$pval.outcome[index] <- NA
		index <- dat$pval.outcome < min_pval
		index[is.na(index)] <- FALSE
		dat$pval.outcome[index] <- min_pval

		dat$pval_origin.outcome <- "reported"
		if(any(is.na(dat$pval.outcome)))
		{
			if("beta.outcome" %in% names(dat) & "se.outcome" %in% names(dat))
			{
				index <- is.na(dat$pval.outcome)
				dat$pval.outcome[index] <- pnorm(abs(dat$beta.outcome[index])/dat$se.outcome[index], lower=FALSE)
				dat$pval_origin.outcome[index] <- "inferred"
			}
		}
	}

	# If no pval column then create it from beta and se if available
	if("beta.outcome" %in% names(dat) & "se.outcome" %in% names(dat) & ! "pval.outcome" %in% names(dat))
	{
		message("Inferring p-values")
		dat$pval.outcome <- pnorm(abs(dat$beta.outcome)/dat$se.outcome, lower=FALSE)
		dat$pval_origin.outcome <- "inferred"
	}

	if(ncase_col %in% names(dat))
	{
		names(dat)[which(names(dat) == ncase_col)[1]] <- "ncase.outcome"
	}
	if(ncontrol_col %in% names(dat))
	{
		names(dat)[which(names(dat) == ncontrol_col)[1]] <- "ncontrol.outcome"
	}
	if(samplesize_col %in% names(dat))
	{
		names(dat)[which(names(dat) == samplesize_col)[1]] <- "samplesize.outcome"
	}

	if(any(dat$mr_keep.outcome))
	{
		mrcols <- c("beta.outcome", "se.outcome", "eaf.outcome", "effect_allele.outcome", "other_allele.outcome")
		dat$mr_keep.outcome <- apply(dat[, mrcols], 1, function(x) !any(is.na(x)))
		if(any(!dat$mr_keep.outcome))
		{
			warning("The following SNP(s) are missing required information for the MR tests and will be excluded\n", paste(subset(dat, !mr_keep.outcome)$SNP, collapse="\n"))
		}
		if(all(!dat$mr_keep.outcome))
		{
			warning("None of the provided SNPs can be used for MR analysis, they are missing required information.")
		}

	}

	dat$id.outcome <- as.numeric(as.factor(dat[[type]]))

	# Get SNP positions if exposure
	if(type == "exposure")
	{
		dat$row_index <- 1:nrow(dat)
		snp <- unique(dat$SNP)

		if(length(snp) > 500)
		{
			message("Looking up SNP info for ", length(snp), " SNPs, this could take some time.")
		}

		bm <- ensembl_get_position(snp)
		missing <- dat$SNP[! dat$SNP %in% bm$refsnp_id]
		if(length(missing) > 0)
		{
			warning("The following SNP(s) were not present in ensembl GRCh37. They will be excluded.", paste(missing, collapse="\n"))
			dat <- subset(dat, SNP %in% bm$refsnp_id)
		}
		stopifnot(nrow(dat) > 0)

		dat <- dat[,names(dat)!="chr_name", drop=FALSE]
		dat <- dat[,names(dat)!="chrom_start", drop=FALSE]
		dat <- merge(dat, bm, by.x="SNP", by.y="refsnp_id", all.x=TRUE)
		dat <- dat[order(dat$row_index), ]
		dat <- subset(dat, select=-c(row_index))
	}

	names(dat) <- gsub("outcome", type, names(dat))

	return(dat)
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
format_gwas_catalog <- function(gwas_catalog_subset, type="exposure")
{
	stopifnot(type %in% c("exposure", "outcome"))


	gwas_catalog_subset[[type]] <- paste(gwas_catalog_subset$Phenotype, gwas_catalog_subset[["Phenotype info"]], gwas_catalog_subset$Units, sep=" || ")
	if(length(unique(gwas_catalog_subset[[type]])) > 1)
	{
		message("Separating the entries into the following phenotypes:\n", paste(unique(gwas_catalog_subset[[type]]), collapse="\n"))
	}

	gwas_catalog_subset <- subset(gwas_catalog_subset, select=c("SNP", "Effect", "eaf", "Allele", "other_allele", "SE", "P-value", type))
	names(gwas_catalog_subset) <- c("SNP", "beta", "eaf", "effect_allele", "other_allele", "se", "pval", type)
	if(type == "exposure")
	{
		format_data(gwas_catalog_subset, type=type, phenotype_col=type)
	} else {

	}
}



ucsc_get_position <- function(snp)
{
	snp <- paste(snp, collapse="', '")
	require(RMySQL)
	message("Connecting to UCSC MySQL database")
	mydb <- dbConnect(MySQL(), user="genome", dbname="hg19", host="genome-mysql.cse.ucsc.edu")

	query <- paste0(
		"SELECT * from snp144 where name in ('", snp, "');"
	)
	message(query)
	out <- dbSendQuery(mydb, query)
	d <- fetch(out, n=-1)
	dbClearResult(dbListResults(mydb)[[1]])
	dbDisconnect(mydb)
	return(d)

}


ensembl_get_position <- function(snp)
{
	library(biomaRt)
	library(stringr)
	Mart <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
	Attr <- listAttributes(Mart)
	ensembl <- getBM(attributes=c("refsnp_id","chr_name","chrom_start","allele","minor_allele","minor_allele_freq"),filters="snp_filter",values=snp,mart=Mart)

	# Sort out chromosome name
	ensembl$chr_name <- str_match(ensembl$chr_name, "(HSCHR)?([0-9X]*)")[,3]
	ensembl$chr_name[ensembl$chr_name == "X"] <- 23
	ensembl$chr_name[! ensembl$chr_name %in% 1:23] <- NA
	ensembl$chr_name <- as.numeric(ensembl$chr_name)

	# Remove SNPs that are problematic:
	# - No chromosome name
	# - No position
	# - Not normal alleles
	# - Not biallelic
	# - No MAF
	# - Duplicate SNPs
	# - minor allele doesn't match the alleles
	# Add major allele
	remove <- is.na(ensembl$chr_name) |
		is.na(ensembl$chrom_start) |
		nchar(ensembl$allele) != 3 |
		! ensembl$minor_allele %in% c("A", "C", "T", "G") |
		is.na(ensembl$minor_allele_freq)
	ensembl <- ensembl[!remove, ]
	ensembl <- subset(ensembl, !duplicated(refsnp_id))
	al <- do.call(rbind, strsplit(ensembl$allele, split="/"))
	i1 <- al[,1] == ensembl$minor_allele
	i2 <- al[,2] == ensembl$minor_allele
	i <- (i1 | i2)
	ensembl <- ensembl[i, ]
	al <- al[i, ]
	i1 <- i1[i]
	i2 <- i2[i]
	ensembl <- subset(ensembl, select=-c(allele))
	ensembl$major_allele[!i1] <- al[!i1, 1]
	ensembl$major_allele[!i2] <- al[!i2, 2]

	return(ensembl)
}
