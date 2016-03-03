#' Perform LD clumping on SNP data
#'
#' Uses PLINK clumping method, where SNPs in LD within a particular window will be pruned. The SNP with the lowest p-value is retained.
#' @param dat Output from \code{read_exposure_data}. Must have a SNP name column (SNP), SNP chromosome column (chr_name), SNP position column (chrom_start) and p-value column (pval.exposure)
#' @param  where Where to perform clumping - if using MR-Base servers then specify "remote" (default), if providing your own reference data and plink binary then specify "local"
#' @param refdat Path to binary plink files to use for LD referencing. If NULL (default) and where="local" then looks for data provided by the Eur1000Genomes R package.
#' @param plink_bin Path to plink binary. If NULL (default) and where="local" then looks for executable in Eur1000Genomes R package.
#' @param clump_kb=10000 Clumping window 
#' @param clump_r2=0.1 Clumping r2 cutoff
#' @param clump_p1=1 Clumping sig level for index SNPs
#' @param clump_p2=1 Clumping sig level for secondary SNPs
#' @param tempdir = "." Where to print results
#' @export
#' @return Data frame
clump_data <- function(dat, where="remote", refdat=NULL, plink_bin=NULL, clump_kb=10000, clump_r2=0.1, clump_p1=1, clump_p2=1, tempdir = getwd())
{
	if(! where %in% c("remote", "local"))
	{
		stop("where variable must be remote or local")
	}
	if(where == "remote")
	{
		return(ld_pruning_api(dat, clump_kb=10000, clump_r2=0.1, clump_p1=1, clump_p2=1))
	} else {
		return(ld_pruning_local(dat, refdat, plink_bin, clump_kb, clump_r2, clump_p1, clump_p2, tempdir))
	}
}


get_plink_exe <- function()
{
    os <- Sys.info()['sysname']
    a <- paste0("exe/plink_", os)
    if(os == "Windows") a <- paste0(a, ".exe")
    plink_bin <- system.file(a, package="Eur1000Genomes")
	if(!file.exists(plink_bin))
	{
		stop("No plink2 executable available for OS '", os, "'. Please provide your own plink2 executable file using the plink_bin argument.")
	}
	return(plink_bin)
}


plink_clump <- function(snps, pvals, refdat, clump_kb, clump_r2, clump_p1, clump_p2, plink_bin, tempdir)
{
	# Make textfile
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	fn <- tempfile(tmpdir = tempdir)
	write.table(data.frame(SNP=snps, P=pvals), file=fn, row=F, col=T, qu=F)

	fun2 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bfile ", shQuote(refdat, type=shell),
		" --clump ", shQuote(fn, type=shell), 
		" --clump-p1 ", clump_p1, 
		" --clump-p2 ", clump_p2, 
		" --clump-r2 ", clump_r2, 
		" --clump-kb ", clump_kb, 
		" --out ", shQuote(fn, type=shell)
	)
	system(fun2)
	a <- read.table(paste(fn, ".clumped", sep=""), he=T)
	unlink(paste(fn, "*", sep=""))
	return(a)
}


#' Perform clumping on the chosen SNPs using through API
#'
#' @param dat Output from \code{read_exposure_data}. Must have a SNP name column (SNP), SNP chromosome column (chr_name), SNP position column (chrom_start) and p-value column (pval.exposure)
#' @param clump_kb=10000 Clumping window 
#' @param clump_r2=0.1 Clumping r2 cutoff
#' @param clump_p1=1 Clumping sig level for index SNPs
#' @param clump_p2=1 Clumping sig level for secondary SNPs
#' @return Data frame of only independent SNPs
ld_pruning_api <- function(dat, clump_kb=10000, clump_r2=0.1, clump_p1=1, clump_p2=1)
{
	snpcode <- paste0("chr", dat$chr_name, ":", dat$chrom_start)
	snpfile <- upload_file_to_api(data.frame(SNP=snpcode, P=dat$pval.exposure), header=TRUE)
	url <- paste0("http://scmv-webapps.epi.bris.ac.uk:5000/clump?snpfile=", snpfile, 
		"&p1=", clump_p1,
		"&p2=", clump_p2,
		"&r2=", clump_r2,
		"&kb=", clump_kb)
	res <- read.table(url, header=TRUE)
	index <- match(res$SNP, snpcode)
	res$SNP <- dat$SNP[index]
	y <- subset(dat, !SNP %in% res$SNP)
	if(nrow(y) > 0)
	{
		message("Removing the following SNPs due to LD with other SNPs:\n", paste(y$SNP, collapse="\n"), sep="\n")
	}
	return(subset(dat, SNP %in% res$SNP))
}


#' Perform clumping on the chosen SNPs using local data
#'
#' @param dat Output from \code{read_exposure_data}
#' @param refdat Path to binary plink files to use for LD referencing. If NULL (default) then looks for data provided by the Eur1000Genomes R package.
#' @param plink_bin Path to plink binary. If NULL (default) then looks for executable in Eur1000Genomes R package.
#' @param clump_kb=10000 Clumping window 
#' @param clump_r2=0.1 Clumping r2 cutoff
#' @param clump_p1=1 Clumping sig level for index SNPs
#' @param clump_p2=1 Clumping sig level for secondary SNPs
#' @param tempdir = "." Where to print results
#'
#' @return Data frame of only independent SNPs
ld_pruning_local <- function(dat, refdat=NULL, clump_kb=10000, clump_r2=0.1, clump_p1=1, clump_p2=1, plink_bin=NULL, tempdir = getwd())
{
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	if(is.null(refdat))
	{
		require(Eur1000Genomes)
		refdat <- gsub(".bed$", "", system.file("data/data_maf0.01.bed", package="Eur1000Genomes"))
	} else {
		stopifnot(file.exists(paste(refdat, ".bed", sep="")))
		stopifnot(file.exists(paste(refdat, ".bim", sep="")))
		stopifnot(file.exists(paste(refdat, ".fam", sep="")))
	}

	if(is.null(plink_bin))
	{
		require(Eur1000Genomes)
		plink_bin <- get_plink_exe()
	} else {
		stopifnot(file.exists(plink_bin))
	}

	snpcode <- paste("chr", dat$chr_name, ":", dat$chrom_start, sep="")

	res <- plink_clump(
		snpcode,
		dat$pval.exposure,
		refdat,
		clump_kb, clump_r2, clump_p1, clump_p2, plink_bin, tempdir)

	index <- match(res$SNP, snpcode)
	res$SNP <- dat$SNP[index]
	y <- subset(dat, !SNP %in% res$SNP)
	if(nrow(y) > 0)
	{
		message("Removing the following SNPs due to LD with other SNPs:", paste(y$SNP, collapse="\n"), sep="\n")
	}
	return(subset(dat, SNP %in% res$SNP, select=-c(pval.exposure)))
}
