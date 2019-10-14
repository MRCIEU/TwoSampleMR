#' Perform LD clumping on SNP data
#'
#' Uses PLINK clumping method, where SNPs in LD within a particular window will be pruned. The SNP with the lowest p-value is retained.
#'
#' @param dat Output from \code{format_data}. Must have a SNP name column (SNP), SNP chromosome column (chr_name), SNP position column (chrom_start). If id.exposure or pval.exposure not present they will be generated
#' @param clump_kb=10000 Clumping window 
#' @param clump_r2=0.001 Clumping r2 cutoff. Note that this default value has recently changed from 0.01.
#' @param clump_p1=1 Clumping sig level for index SNPs
#' @param clump_p2=1 Clumping sig level for secondary SNPs
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data
#'
#' @export
#' @return Data frame
clump_data <- function(dat, clump_kb=10000, clump_r2=0.001, clump_p1=1, clump_p2=1, access_token=get_mrbase_access_token())
{
	.Deprecated("ieugwasr::ld_clump")

	if(missing(clump_r2))
	{
		message("Warning: since v0.4.2 the default r2 value for clumping has changed from 0.01 to 0.001")
	}
	if(!is.data.frame(dat))
	{
		stop("Expecting data frame returned from format_data")
	}

	if(! "pval.exposure" %in% names(dat))
	{
		dat$pval.exposure <- 0.99
	}

	if(! "id.exposure" %in% names(dat))
	{
		dat$id.exposure <- random_string(1)
	}

	res <- plyr::ddply(dat, c("id.exposure"), function(x)
	{
		x <- plyr::mutate(x)
		if(nrow(x) == 1)
		{
			message("Only one SNP for ", x$id.exposure[1])
			return(x)
		} else {
			message("Clumping ", x$id.exposure[1], ", ", nrow(x), " SNPs")
			return(ld_pruning_api(x, clump_kb=clump_kb, clump_r2=clump_r2, clump_p1=clump_p1, clump_p2=clump_p2, access_token=access_token))
		}
	})
	return(res)
}


#' Perform clumping on the chosen SNPs using through API
#'
#' @param dat Output from \code{read_exposure_data}. Must have a SNP name column (SNP), SNP chromosome column (chr_name), SNP position column (chrom_start) and p-value column (pval.exposure)
#' @param clump_kb=10000 Clumping window 
#' @param clump_r2=0.1 Clumping r2 cutoff
#' @param clump_p1=1 Clumping sig level for index SNPs
#' @param clump_p2=1 Clumping sig level for secondary SNPs
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data#' @return Data frame of only independent SNPs
ld_pruning_api <- function(dat, clump_kb=10000, clump_r2=0.1, clump_p1=1, clump_p2=1, access_token=get_mrbase_access_token())
{
	res <- api_query('ld/clump',
			query = list(
				rsid = dat$SNP,
				pval = dat$pval.exposure,
				pthresh = clump_p1,
				r2 = clump_r2,
				kb = clump_kb
			),
			access_token=access_token
		)
	y <- subset(dat, !SNP %in% res)
	if(nrow(y) > 0)
	{
		message("Removing the following SNPs due to LD with other SNPs:\n", paste(y$SNP, collapse="\n"), sep="\n")
	}
	return(subset(dat, SNP %in% res))
}


#' Get LD matrix for list of SNPs
#'
#' This function takes a list of SNPs and searches for them in 502 European samples from 1000 Genomes phase 3 data
#' It then creates an LD matrix of r values (signed, and not squared)
#' All LD values are with respect to the major alleles in the 1000G dataset. You can specify whether the allele names are displayed
#'
#' @param snps List of SNPs
#' @param with_alleles Whether to append the allele names to the SNP names. Default: TRUE
#'
#' @export
#' @return Matrix of LD r values
ld_matrix <- function(snps, with_alleles=TRUE)
{
	.Deprecated("ieugwasr::ld_matrix")
	if(length(snps) > 500)
	{
		stop("SNP list must be smaller than 500")
	}

	res <- api_query('ld/matrix', query = list(rsid=snps), access_token=NULL)

	if(all(is.na(res))) stop("None of the requested SNPs were found")
	snps2 <- res$snplist
	res <- res$matrix
	res <- matrix(as.numeric(res), nrow(res), ncol(res))
	snps3 <- do.call(rbind, strsplit(snps2, split="_"))
	if(with_alleles)
	{
		rownames(res) <- snps2
		colnames(res) <- snps2
	} else {
		rownames(res) <- snps3[,1]
		colnames(res) <- snps3[,1]
	}
	missing <- snps[!snps %in% snps3[,1]]
	if(length(missing) > 0)
	{
		warning("The following SNPs are not present in the LD reference panel\n", paste(missing, collapse="\n"))
	}
	ord <- match(snps3[,1], snps)
	res <- res[order(ord), order(ord)]
	return(res)
}


