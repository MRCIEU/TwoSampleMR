#' Perform LD clumping on SNP data
#'
#' Uses PLINK clumping method, where SNPs in LD within a particular window will be pruned. The SNP with the lowest p-value is retained.
#' 
#' @details
#' This function interacts with the OpenGWAS API, which houses LD reference panels for the 5 super-populations in the 1000 genomes reference panel. It includes only bi-allelic SNPs with MAF > 0.01, so it's quite possible that a variant you want to include in the clumping process will be absent. If it is absent, it will be automatically excluded from the results.
#' 
#' You can check if your variants are present in the LD reference panel using ieugwasr::ld_reflookup()
#'
#' This function does put load on the OpenGWAS servers, which makes life more difficult for other users. We have implemented a method and made available the LD reference panels to perform clumping locally, see ieugwasr::ld_clump() and related vignettes for details.
#'
#' @md
#' @param dat Output from [`format_data`]. Must have a SNP name column (SNP), SNP chromosome column (chr_name), SNP position column (chrom_start). If id.exposure or pval.exposure not present they will be generated.
#' @param clump_kb Clumping window, default is `10000`.
#' @param clump_r2 Clumping r2 cutoff. Note that this default value has recently changed from `0.01` to `0.001`.
#' @param clump_p1 Clumping sig level for index SNPs, default is `1`.
#' @param clump_p2 Clumping sig level for secondary SNPs, default is `1`.
#' @param pop Super-population to use as reference panel. Default = "EUR". Options are EUR, SAS, EAS, AFR, AMR. 'legacy' also available - which is a previously used verison of the EUR panel with a slightly different set of markers
#'
#' @export
#' @return Data frame
clump_data <- function(dat, clump_kb=10000, clump_r2=0.001, clump_p1=1, clump_p2=1, pop="EUR")
{
	# .Deprecated("ieugwasr::ld_clump()")

	pval_column <- "pval.exposure"

	if(!is.data.frame(dat))
	{
		stop("Expecting data frame returned from format_data")
	}

	if("pval.exposure" %in% names(dat) & "pval.outcome" %in% names(dat))
	{
		message("pval.exposure and pval.outcome columns present. Using pval.exposure for clumping.")
	} else if(!"pval.exposure" %in% names(dat) & "pval.outcome" %in% names(dat))
	{
		message("pval.exposure column not present, using pval.outcome column for clumping.")
		pval_column <- "pval.outcome"
	} else if(! "pval.exposure" %in% names(dat))
	{
		message("pval.exposure not present, setting clumping p-value to 0.99 for all variants")
		dat$pval.exposure <- 0.99
	} else {
		pval_column <- "pval.exposure"
	}
	
	if(! "id.exposure" %in% names(dat))
	{
		dat$id.exposure <- random_string(1)
	}

	d <- data.frame(rsid=dat$SNP, pval=dat[[pval_column]], id=dat$id.exposure)
	out <- ieugwasr::ld_clump(d, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p1, pop=pop)
	keep <- paste(dat$SNP, dat$id.exposure) %in% paste(out$rsid, out$id)
	return(dat[keep, ])
}


#' Get LD matrix for list of SNPs
#'
#' This function takes a list of SNPs and searches for them in a specified super-population in the 1000 Genomes phase 3 reference panel.
#' It then creates an LD matrix of r values (signed, and not squared).
#' All LD values are with respect to the major alleles in the 1000G dataset. You can specify whether the allele names are displayed.
#'
#' @details
#' The data used for generating the LD matrix includes only bi-allelic SNPs with MAF > 0.01, so it's quite possible that a variant you want to include will be absent. If it is absent, it will be automatically excluded from the results.
#' 
#' You can check if your variants are present in the LD reference panel using ieugwasr::ld_reflookup()
#'
#' This function does put load on the OpenGWAS servers, which makes life more difficult for other users, and has been limited to analyse only up to 500 variants at a time. We have implemented a method and made available the LD reference panels to perform the operation locally, see ieugwasr::ld_matrix() and related vignettes for details.
#'
#' @param snps List of SNPs.
#' @param with_alleles Whether to append the allele names to the SNP names. The default is `TRUE`.
#' @param pop Super-population to use as reference panel. Default = "EUR". Options are EUR, SAS, EAS, AFR, AMR. 'legacy' also available - which is a previously used verison of the EUR panel with a slightly different set of markers
#'
#' @export
#' @return Matrix of LD r values
ld_matrix <- function(snps, with_alleles=TRUE, pop="EUR")
{
	# .Deprecated("ieugwasr::ld_matrix()")
	ieugwasr::ld_matrix(variants=snps, with_alleles=with_alleles, pop=pop)
}

