plink_clump <- function(snps, pvals, refdat, clump_kb, clump_r2, clump_p1, clump_p2, plink_bin, tempdir)
{
	# Make textfile
	fn <- tempfile(tmpdir = tempdir)
	write.table(data.frame(SNP=snps, P=pvals), file=fn, row=F, col=T, qu=F)

	fun2 <- paste(plink_bin, " --bfile ", refdat, " --clump ", fn, " --clump-p1 ", clump_p1, " --clump-p2 ", clump_p2, " --clump-r2 ", clump_r2, " --clump-kb ", clump_kb, " --out ", fn, sep="")
	system(fun2)
	a <- read.table(paste(fn, ".clumped", sep=""), he=T)
	unlink(paste(fn, "*", sep=""))
	return(a)
}

#' Perform clumping on the chosen SNPs
#'
#' @param dat Output from \code{read_exposure_data}
#' @param refdat="data_maf0.01" Location of reference genotype data
#' @param clump_kb=10000 Clumping window 
#' @param clump_r2=0.1 Clumping r2 cutoff
#' @param clump_p1=1 Clumping sig level for index SNPs
#' @param clump_p2=1 Clumping sig level for secondary SNPs
#' @param plink_bin="plink1.90" plink2 exe
#' @param tempdir = "." Where to print results
#'
#' @export
#' @return Data frame of only independent SNPs
ld_pruning_all <- function(dat, refdat="data_maf0.01", clump_kb=10000, clump_r2=0.1, clump_p1=1, clump_p2=1, plink_bin="plink1.90", tempdir = ".")
{
	dat$p.exposure <- pnorm(abs(dat$beta.exposure / dat$se.exposure), lower.tail=FALSE)
	snpcode <- paste("chr", dat$chr_name, ":", dat$chrom_start, sep="")
	res <- plink_clump(
		snpcode,
		dat$p.exposure,
		refdat,
		clump_kb, clump_r2, clump_p1, clump_p2, plink_bin, tempdir)

	index <- match(res$SNP, snpcode)
	res$SNP <- dat$SNP[index]
	y <- subset(dat, !SNP %in% res$SNP)
	if(nrow(y) > 0)
	{
		message("Removing the following SNPs due to LD with other SNPs:", paste(y$SNP, collapse="\n"), sep="\n")
	}
	return(subset(dat, SNP %in% res$SNP))
}


