plink_clump <- function(snps, pvals, refdat, clump-kb, clump-r2, clump-p1, clump-p2, plink-bin, tempdir)
{
	# Make textfile
	fn <- tempfile(tmpdir = tempdir)
	write.table(data.frame(SNP=snps, P=pvals), file=fn, row=F, col=T, qu=F)

	fun1 <- paste(plink-bin, " --bfile ", refdat, " --extract ", fn, " --make-bed --out ", fn, sep="")
	system(fun1)

	fun2 <- paste(plink-bin, " --bfile ", fn, " --clump ", fn, " --clump-p1 ", clump-p1, " --clump-p2 ", clump-p2, " --clump-r2 ", clump-r2, " --clump-kb ", clump-kb, " --out ", fn, sep="")
	system(fun2)
	a <- read.table(paste(fn, ".clumped", sep=""), he=T)
	unlink(paste(fn, "*", sep=""))
	return(a)
}

# Use clumping in plink2
ld_pruning <- function(dat, refdat=paste("/path/chr", 1:22, sep=""), clump-kb=10000, clump-r2=0.1, clump-p1=1, clump-p2=1, plink-bin, tempdir = ".")
{

	dat$p.exposure <- pnorm(abs(dat$beta.exposure / dat$se.exposure), lower.tail=FALSE)
	ddply(dat, .(chr_name), function(x)
	{
		res <- plink_clump(
			paste(x$chr_name, ":", x$chrom_start, sep=""),
			dat$p.exposure,
			refdat[as.numeric(x$chr_name)],
			clump-kb, clump-r2, clump-p1, clump-p2, plink-bin, tempdir)
		y <- subset(x, SNP %in% res$SNP)
		if(nrow(y) > 0)
		{
			message("Removing", y$SNP, "using clumping due to LD with other SNPs")
		}
		return(subset(x, SNP %in% res$SNP))
	})
}
