#' Perform LD clumping on SNP data
#'
#' Uses PLINK clumping method, where SNPs in LD within a particular window will be pruned. The SNP with the lowest p-value is retained.
#'
#' @param dat Output from \code{format_data}. Must have a SNP name column (SNP), SNP chromosome column (chr_name), SNP position column (chrom_start). If id.exposure or pval.exposure not present they will be generated
#' @param clump_kb=10000 Clumping window 
#' @param clump_r2=0.1 Clumping r2 cutoff
#' @param clump_p1=1 Clumping sig level for index SNPs
#' @param clump_p2=1 Clumping sig level for secondary SNPs
#'
#' @export
#' @return Data frame
clump_data <- function(dat, clump_kb=10000, clump_r2=0.1, clump_p1=1, clump_p2=1)
{
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

	res <- ddply(dat, .(id.exposure), function(x)
	{
		x <- mutate(x)
		message("Clumping ", x$id.exposure[1], ", ", nrow(x), " SNPs")
		return(ld_pruning_api(x, clump_kb=clump_kb, clump_r2=clump_r2, clump_p1=clump_p1, clump_p2=clump_p2))
	})
	return(res)
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
	snpfile <- upload_file_to_api(data.frame(SNP=dat$SNP, P=dat$pval.exposure), header=TRUE)
	url <- paste0(options()$mrbaseapi, "clump?snpfile=", snpfile,
		"&p1=", clump_p1,
		"&p2=", clump_p2,
		"&r2=", clump_r2,
		"&kb=", clump_kb)
	res <- fromJSON_safe(url)
	y <- subset(dat, !SNP %in% res$SNP)
	if(nrow(y) > 0)
	{
		message("Removing the following SNPs due to LD with other SNPs:\n", paste(y$SNP, collapse="\n"), sep="\n")
	}
	return(subset(dat, SNP %in% res$SNP))
}



get_snp_positions_biomart <- function(dat)
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
	return(dat)
}

