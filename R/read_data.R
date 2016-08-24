

#' Read outcome data
#'
#' Reads in outcome data. Checks and organises columns for use with MR or enrichment tests. Infers p-values when possible from beta and se.
#'
#' @param filename Filename. Must have header with at least SNP column present.
#' @param snps Default=NULL. SNPs to extract. If NULL then doesn't extract any and keeps all.
#' @param sep Default=" ". Specify delimeter in file
#' @param phenotype_col Default="Phenotype". Optional column name for the column with phenotype name corresponding the the SNP. If not present then will be created with the value "Outcome"
#' @param snp_col Default="SNP". Required name of column with SNP rs IDs
#' @param beta_col Default="beta". Required for MR. Name of column with effect sizes
#' @param se_col Default="se". Required for MR. Name of column with standard errors
#' @param eaf_col Default="eaf". Required for MR. Name of column with effect allele frequency
#' @param effect_allele_col Default="effect_allele". Required for MR. Name of column with effect allele. Must be "A", "C", "T" or "G"
#' @param other_allele_col Default="other_allele". Required for MR. Name of column with non effect allele. Must be "A", "C", "T" or "G"
#' @param pval_col Default="pval". Required for enrichment tests. Name of column with p-value.
#' @param units_col Default="units". Optional column name for units.
#' @param ncase_col Default="ncase". Optional column name for number of cases.
#' @param ncontrol_col Default="ncontrol". Optional column name for number of controls.
#' @param samplesize_col Default="samplesize". Optional column name for sample size.
#' @param gene_col Default="gene". Optional column name for gene name.
#' @param min_pval Default=1e-200. Minimum allowed pval
#'
#' @export
#' @return data frame
read_outcome_data <- function(filename, snps=NULL, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", units_col="units", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", min_pval=1e-200)
{
	outcome_dat <- data.table::fread(filename, header=TRUE, sep=sep)
	outcome_dat <- format_data(
		as.data.frame(outcome_dat),
		type="outcome",
		snps=snps,
		phenotype_col=phenotype_col,
		snp_col=snp_col,
		beta_col=beta_col,
		se_col=se_col,
		eaf_col=eaf_col,
		effect_allele_col=effect_allele_col,
		other_allele_col=other_allele_col,
		pval_col=pval_col,
		units_col=units_col,
		ncase_col=ncase_col,
		ncontrol_col=ncontrol_col,
		samplesize_col=samplesize_col,
		gene_col=gene_col,
		min_pval=min_pval
	)
	outcome_dat$data_source.outcome <- "textfile"
	return(outcome_dat)
}

#' Read exposure data
#'
#' Reads in exposure data. Checks and organises columns for use with MR or enrichment tests. Infers p-values when possible from beta and se. Looks up SNPs in biomaRt to get basic info.
#'
#' @param filename Filename. Must have header with at least SNP column present.
#' @param sep Default=" ". Specify delimeter in file
#' @param phenotype_col Default="Phenotype". Optional column name for the column with phenotype name corresponding the the SNP. If not present then will be created with the value "Outcome"
#' @param snp_col Default="SNP". Required name of column with SNP rs IDs
#' @param beta_col Default="beta". Required for MR. Name of column with effect sizes
#' @param se_col Default="se". Required for MR. Name of column with standard errors
#' @param eaf_col Default="eaf". Required for MR. Name of column with effect allele frequency
#' @param effect_allele_col Default="effect_allele". Required for MR. Name of column with effect allele. Must be "A", "C", "T" or "G"
#' @param other_allele_col Default="other_allele". Required for MR. Name of column with non effect allele. Must be "A", "C", "T" or "G"
#' @param pval_col Default="pval". Required for enrichment tests. Name of column with p-value.
#' @param units_col Default="units". Optional column name for units.
#' @param ncase_col Default="ncase". Optional column name for number of cases.
#' @param ncontrol_col Default="ncontrol". Optional column name for number of controls.
#' @param samplesize_col Default="samplesize". Optional column name for sample size.
#' @param gene_col Default="gene". Optional column name for gene name.
#' @param min_pval Default=1e-200. Minimum allowed pval
#'
#' @export
#' @return data frame
read_exposure_data <- function(filename, clump=FALSE, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", units_col="units", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", min_pval=1e-200)
{
	exposure_dat <- data.table::fread(filename, header=TRUE, sep=sep)
	exposure_dat <- format_data(
		as.data.frame(exposure_dat),
		type="exposure",
		snps=NULL,
		phenotype_col=phenotype_col,
		snp_col=snp_col,
		beta_col=beta_col,
		se_col=se_col,
		eaf_col=eaf_col,
		effect_allele_col=effect_allele_col,
		other_allele_col=other_allele_col,
		pval_col=pval_col,
		units_col=units_col,
		ncase_col=ncase_col,
		ncontrol_col=ncontrol_col,
		samplesize_col=samplesize_col,
		gene_col=gene_col,
		min_pval=min_pval
	)
	exposure_dat$data_source.exposure <- "textfile"
	if(clump)
	{
		exposure_dat <- clump_data(exposure_dat)
	}
	return(exposure_dat)
}

#' Read exposure or outcome data
#'
#' Reads in exposure data. Checks and organises columns for use with MR or enrichment tests. Infers p-values when possible from beta and se. If it is the exposure then looks up SNPs in biomaRt to get basic info.
#'
#' @param dat Data frame. Must have header with at least SNP column present.
#' @param type="exposure". Is this the exposure or the outcome data that is being read in?
#' @param snps=NULL SNPs to extract. If NULL then doesn't extract any and keeps all.
#' @param phenotype_col="Phenotype" Optional column name for the column with phenotype name corresponding the the SNP. If not present then will be created with the value "Outcome"
#' @param snp_col="SNP" Required name of column with SNP rs IDs
#' @param beta_col="beta" Required for MR. Name of column with effect sizes
#' @param se_col="se" Required for MR. Name of column with standard errors
#' @param eaf_col="eaf" Required for MR. Name of column with effect allele frequency
#' @param effect_allele_col="effect_allele" Required for MR. Name of column with effect allele. Must be "A", "C", "T" or "G"
#' @param other_allele_col="other_allele" Required for MR. Name of column with non effect allele. Must be "A", "C", "T" or "G"
#' @param pval_col="pval" Required for enrichment tests. Name of column with p-value.
#' @param units_col="units" Optional column name for units.
#' @param ncase_col="ncase" Optional column name for number of cases.
#' @param ncontrol_col="ncontrol" Optional column name for number of controls.
#' @param samplesize_col="samplesize" Optional column name for sample size.
#' @param gene_col="gene" Optional column name for gene name.
#' @param min_pval=1e-200 Minimum allowed pval
#'
#' @export
#' @return data frame
format_data <- function(dat, type="exposure", snps=NULL, header=TRUE, phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", units_col="units", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", min_pval=1e-200)
{
	all_cols <- c(phenotype_col, snp_col, beta_col, se_col, eaf_col, effect_allele_col, other_allele_col, pval_col, units_col, ncase_col, ncontrol_col, samplesize_col, gene_col)

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
	snp_col <- "SNP"
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
			dat <- dat[,-which(names(dat)==phenotype_col)]
		}
	}

	# Remove duplicated SNPs
	dat <- plyr::ddply(dat, type, function(x){
		x <- plyr::mutate(x)
		dup <- duplicated(x$SNP)
		if(any(dup))
		{
			warning("Duplicated SNPs present in exposure data for phenotype '", x[[type]][1], ". Just keeping the first instance:\n", paste(x$SNP[dup], collapse="\n"))
			x <- x[!dup,]
		}
		return(x)		
	})

	# Check if columns required for MR are present
	mr_cols_required <- c(snp_col, beta_col, se_col, effect_allele_col) 
	mr_cols_desired <- c(other_allele_col, eaf_col)
	if(! all(mr_cols_required %in% names(dat)))
	{
		warning("The following columns are not present and are required for MR analysis\n", paste(mr_cols_required[!mr_cols_required %in% names(dat)]), collapse="\n")
		dat$mr_keep.outcome <- FALSE
	} else {
		dat$mr_keep.outcome <- TRUE
	}

	if(! all(mr_cols_desired %in% names(dat)))
	{
		warning("The following columns are not present but are helpful for harmonisation\n", paste(mr_cols_desired[!mr_cols_desired %in% names(dat)]), collapse="\n")
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
		if(is.logical(dat$effect_allele.outcome))
		{
			dat$effect_allele.outcome <- substr(as.character(dat$effect_allele.outcome), 1, 1)
		}
		if(!is.character(dat$effect_allele.outcome))
		{
			warning("effect_allele column is not character data. Coercing...")
			dat$effect_allele.outcome <- as.character(dat$effect_allele.outcome)
		}

		dat$effect_allele.outcome <- toupper(dat$effect_allele.outcome)
		index <- ! dat$effect_allele.outcome %in% c("A", "C", "T", "G")
		index[is.na(index)] <- TRUE
		if(any(index))
		{
			warning("effect_allele column has some values that are not A/C/T/G. These SNPs will be excluded")
			dat$effect_allele.outcome[index] <- NA
			dat$mr_keep.outcome[index] <- FALSE
		}
	}


	# Check other_allele
	i <- which(names(dat) == other_allele_col)[1]
	if(!is.na(i))
	{
		names(dat)[i] <- "other_allele.outcome"
		if(is.logical(dat$other_allele.outcome))
		{
			dat$other_allele.outcome <- substr(as.character(dat$other_allele.outcome), 1, 1)
		}
		if(!is.character(dat$other_allele.outcome))
		{
			warning("other_allele column is not character data. Coercing...")
			dat$other_allele.outcome <- as.character(dat$other_allele.outcome)
		}

		dat$other_allele.outcome <- toupper(dat$other_allele.outcome)
		index <- ! dat$other_allele.outcome %in% c("A", "C", "T", "G")
		index[is.na(index)] <- TRUE
		if(any(index))
		{
			warning("other_allele column has some values that are not A/C/T/G. These SNPs will be excluded")
			dat$other_allele.outcome[index] <- NA
			dat$mr_keep.outcome[index] <- FALSE
		}
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
		if(!is.numeric(dat$ncase.outcome))
		{
			warning(ncase_col, " column is not numeric")
			dat$ncase.outcome <- as.numeric(dat$ncase.outcome)
		}
	}
	if(ncontrol_col %in% names(dat))
	{
		names(dat)[which(names(dat) == ncontrol_col)[1]] <- "ncontrol.outcome"
		if(!is.numeric(dat$ncontrol.outcome))
		{
			warning(ncontrol_col, " column is not numeric")
			dat$ncontrol.outcome <- as.numeric(dat$ncontrol.outcome)
		}
	}

	

	if(samplesize_col %in% names(dat))
	{
		names(dat)[which(names(dat) == samplesize_col)[1]] <- "samplesize.outcome"
		if(!is.numeric(dat$samplesize.outcome))
		{
			warning(samplesize_col, " column is not numeric")
			dat$samplesize.outcome <- as.numeric(dat$samplesize.outcome)
		}

		if("ncontrol.outcome" %in% names(dat) & "ncase.outcome" %in% names(dat))
		{
			index <- is.na(dat$samplesize.outcome) & !is.na(dat$ncase.outcome) & !is.na(dat$ncontrol.outcome)
			if(any(index))
			{
				message("Generating sample size from ncase and ncontrol")
				dat$samplesize.outcome[index] <- dat$ncase.outcome[index] + dat$ncontrol.outcome[index]			
			}
		}
	} else if("ncontrol.outcome" %in% names(dat) & "ncase.outcome" %in% names(dat))
	{
		message("Generating sample size from ncase and ncontrol")
		dat$samplesize.outcome <- dat$ncase.outcome + dat$ncontrol.outcome
	}

	if(gene_col %in% names(dat))
	{
		names(dat)[which(names(dat) == gene_col)[1]] <- "gene.outcome"
	}
	

	if(units_col %in% names(dat))
	{
		names(dat)[which(names(dat) == units_col)[1]] <- "units.outcome"
		dat$units.outcome_dat <- as.character(dat$units.outcome)
		temp <- check_units(dat, type, "units.outcome")
		if(any(temp$ph))
		{
			dat[[type]] <- paste0(dat[[type]], " (", dat$units.outcome, ")")
		}
	}

	# Create id column
	dat$id.outcome <- create_ids(dat[[type]])

	if(any(dat$mr_keep.outcome))
	{
		mrcols <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome")
		mrcols_present <- mrcols[mrcols %in% names(dat)]
		dat$mr_keep.outcome <- apply(dat[, mrcols_present], 1, function(x) !any(is.na(x)))
		if(any(!dat$mr_keep.outcome))
		{
			warning("The following SNP(s) are missing required information for the MR tests and will be excluded\n", paste(subset(dat, !mr_keep.outcome)$SNP, collapse="\n"))
		}
	}
	if(all(!dat$mr_keep.outcome))
	{
		warning("None of the provided SNPs can be used for MR analysis, they are missing required information.")
	}


	# Add in missing MR cols
	for(col in c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome"))
	{
		if(! col %in% names(dat))
		{
			dat[[col]] <- NA
		}
	}

	names(dat) <- gsub("outcome", type, names(dat))
	rownames(dat) <- NULL
	return(dat)
}

check_units <- function(x, id, col)
{
	temp <- plyr::ddply(x, id, function(x1)
	{
		ph <- FALSE
		if(length(unique(x1[[col]])) > 1)
		{
			warning("More than one type of unit specified for ", x1[[id]][1])
			x1 <- plyr::mutate(x1)
			ph <- TRUE
		}
		return(data.frame(ph=ph[1], stringsAsFactors=FALSE))
	})
	return(temp)
}


#' Get data selected from GWAS catalog into correct format
#'
#' Subset the GWAS catalogue to have the rows you require for instrumenting a particular exposure and then run this command.
#' Be careful to avoid using different phenotypes, phenotype types, or units together.
#'
#' @param gwas_catalog_subset Subset of rows from \code{data(gwas_catalog)}
#' @param type Are these data used as "exposure" or "outcome"? Default is "exposure"
#'
#' @export
#' @return Data frame
#' @examples \dontrun{
#' data(gwas_catalog)
#' require(MRInstruments)
#' bmi <- subset(gwas_catalog, Phenotype=="Body mass index" & Year==2010 & grepl("kg", Units)
#' bmi <- format_gwas_catalog(bmi)
#'}
format_gwas_catalog <- function(gwas_catalog_subset, type="exposure")
{
	stopifnot(type %in% c("exposure", "outcome"))

	if(length(unique(gwas_catalog_subset$Phenotype)) > 1)
	{
		message("Separating the entries into the following phenotypes:\n", paste(unique(gwas_catalog_subset$Phenotype), collapse="\n"))
	}
	
	dat <- format_data(gwas_catalog_subset, type=type)
	dat[[paste0("data_source.", type)]] <- "gwas_catalog"

	return(dat)
}


#' Get data from eQTL catalog into correct format
#'
#' See \code{format_data}
#'
#' @param gtex_eqtl_subset Selected rows from \code{gtex_eqtl} data loaded from \code{MRInstruments} package
#' @param type Are these data used as "exposure" or "outcome"? Default is "exposure"
#'
#' @export
#' @return Data frame
format_gtex_eqtl <- function(gtex_eqtl_subset, type="exposure")
{
	stopifnot(type %in% c("exposure", "outcome"))
	gtex_eqtl_subset[[type]] <- paste0(gtex_eqtl_subset$gene_name, " (", gtex_eqtl_subset$tissue, ")")

	if(length(unique(gtex_eqtl_subset[[type]])) > 1)
	{
		message("Separating the entries into the following phenotypes:\n", paste(unique(gtex_eqtl_subset[[type]]), collapse="\n"))
	}

	dat <- format_data(gtex_eqtl_subset, type=type, phenotype_col=type, pval_col="P_value", samplesize_col="n")
	dat[[paste0("data_source.", type)]] <- "gtex_eqtl"

	return(dat)

}


#' Get data from metabolomic QTL results
#'
#' See \code{format_data}
#'
#' @param metab_qtls_subset Selected rows from \code{metab_qtls} data loaded from \code{MRInstruments} package
#' @param type Are these data used as "exposure" or "outcome"? Default is "exposure"
#'
#' @export
#' @return Data frame
format_metab_qtls <- function(metab_qtls_subset, type="exposure")
{
	stopifnot(type %in% c("exposure", "outcome"))
	

	if(length(unique(metab_qtls_subset$phenotype)) > 1)
	{
		message("Separating the entries into the following phenotypes:\n", paste(unique(metab_qtls_subset$phenotype), collapse="\n"))
	}

	dat <- format_data(metab_qtls_subset, type=type, phenotype_col="phenotype")
	dat[[paste0("data_source.", type)]] <- "metab_qtls"

	return(dat)
}



#' Get data from proteomic QTL results
#'
#' See \code{format_data}
#'
#' @param proteomic_qtls_subset Selected rows from \code{proteomic_qtls} data loaded from \code{MRInstruments} package
#' @param type Are these data used as "exposure" or "outcome"? Default is "exposure"
#'
#' @export
#' @return Data frame
format_proteomic_qtls <- function(proteomic_qtls_subset, type="exposure")
{
	stopifnot(type %in% c("exposure", "outcome"))
	

	if(length(unique(proteomic_qtls_subset$analyte)) > 1)
	{
		message("Separating the entries into the following phenotypes:\n", paste(unique(proteomic_qtls_subset$analyte), collapse="\n"))
	}

	dat <- format_data(proteomic_qtls_subset, type=type, phenotype_col="analyte")
	dat[[paste0("data_source.", type)]] <- "proteomic_qtls"

	return(dat)
}



#' Get data from methylation QTL results
#'
#' See \code{format_data}
#'
#' @param aries_mqtl_subset Selected rows from \code{aries_mqtl} data loaded from \code{MRInstruments} package
#' @param type Are these data used as "exposure" or "outcome"? Default is "exposure"
#'
#' @export
#' @return Data frame
format_aries_mqtl <- function(aries_mqtl_subset, type="exposure")
{
	stopifnot(type %in% c("exposure", "outcome"))
	
	aries_mqtl_subset$Phenotype <- paste0(aries_mqtl_subset$cpg, " (", aries_mqtl_subset$age, ")")

	if(length(unique(aries_mqtl_subset$Phenotype)) > 1)
	{
		message("Separating the entries into the following phenotypes:\n", paste(unique(aries_mqtl_subset$Phenotype), collapse="\n"))
	}

	dat <- format_data(aries_mqtl_subset, type=type)
	dat[[paste0("data_source.", type)]] <- "aries_mqtl"

	return(dat)
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
	require(biomaRt)
	require(stringr)
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
	ensembl <- ensembl[!remove, , drop=FALSE]
	ensembl <- subset(ensembl, !duplicated(refsnp_id))
	al <- do.call(rbind, strsplit(ensembl$allele, split="/"))
	i1 <- al[,1] == ensembl$minor_allele
	i2 <- al[,2] == ensembl$minor_allele
	i <- (i1 | i2)
	ensembl <- ensembl[i, ]
	al <- al[i, , drop=FALSE]
	i1 <- i1[i]
	i2 <- i2[i]
	ensembl <- subset(ensembl, select=-c(allele))
	ensembl$major_allele[!i1] <- al[!i1, 1]
	ensembl$major_allele[!i2] <- al[!i2, 2]

	return(ensembl)
}


random_string <- function(n=1, len=6)
{
	randomString <- c(1:n)
	for (i in 1:n)
	{
		randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
		len, replace=TRUE),
		collapse="")
	}
	return(randomString)
}

create_ids <- function(x)
{
	a <- as.factor(x)
	levels(a) <- random_string(length(levels(a)))
	a <- as.character(a)
	return(a)
}


#' Combine data
#'
#' Taking exposure or outcome data (returned from \code{format_data})
#' combine multiple datasets together so they can be analysed in one
#' batch. Removes duplicate SNPs, preferentially keeping those usable
#' in MR analysis
#'
#' @param x List of data frames returned from \code{format_data}
#'
#' @export
#' @return data frame
combine_data <- function(x)
{
	stopifnot(is.list(x))
	if("exposure" %in% names(x[[1]])) type <- "exposure"
	else if("outcome" %in% names(x[[1]])) type <- "outcome"
	else stop("Datasets must be generated from format_data")

	check <- all(sapply(x, function(i) {
		type %in% names(i)}))

	if(!check)
	{
		stop("Not all datasets or of type '", type, "'")
	}

	id_col <- paste0("id.", type)
	mr_keep_col <- paste0("mr_keep.", type)
	x <- plyr::rbind.fill(x)

	x <- plyr::ddply(x, id_col, function(x)
	{
		x <- plyr::mutate(x)
		x <- x[order(x[[mr_keep_col]], decreasing=TRUE), ]
		x <- subset(x, !duplicated(SNP))
		return(x)
	})

	rownames(x) <- NULL
	return(x)
}



