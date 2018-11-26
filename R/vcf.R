#' Check a GWAS dataset against a reference known to be on the forward strand
#'
#' 
#' Assuming reference data is all on forward strand, check if 
#' the GWAS is also.
#' Use some threshold e.g. if more than 90% of alleles don't 
#' need to be flipped then it's likely that the dataset is on
#' the forward strand
#'
#' This function can be used to evaluate how strict harmonisation should be
#' The trade off if you assume we are not on the forward strand then palindromic SNPs are dropped within a particular frequency range
#' But you could instead have some small probability of error for whether palindromic SNPs are on the forward strand, and avoid dropping too many variants.
#'
#' @param gwas_snp Vector of SNP names for the dataset being checked
#' @param gwas_a1 Vector of alleles
#' @param gwas_a2 Vector of alleles
#' @param ref_snp Vector of SNP names for the reference dataset
#' @param ref_a1 Vector of alleles
#' @param ref_a2 Vector of alleles
#' @param threshold=0.9 If the proportion of allele strands match is above this threshold, then declare the dataset to be on the forward strand
#'
#' @export
#' @return 1 = Forward strand; 2 = Not on forward strand
is_forward_strand <- function(gwas_snp, gwas_a1, gwas_a2, ref_snp, ref_a1, ref_a2, threshold=0.9)
{
	requireNamespace("dplyr", quietly=TRUE)
	if(is.null(gwas_a1) | is.null(gwas_a2))
	{
		message("No info for both GWAS alleles")
		return(2)
	}

	if(1-(sum(is.na(gwas_a1)) / length(gwas_a1)) < threshold)
	{
		message("Too many missing values for gwas A1")
		return(2)
	}
	if(1-(sum(is.na(gwas_a2)) / length(gwas_a2)) < threshold)
	{
		message("Too many missing values for gwas A2")
		return(2)
	}

	g <- dplyr::data_frame(SNP=gwas_snp, A1=toupper(gwas_a1), A2=toupper(gwas_a2))
	r <- dplyr::data_frame(SNP=ref_snp, A1=toupper(ref_a1), A2=toupper(ref_a2))

	gr <- dplyr::inner_join(g,r,by="SNP")
	index <- (gr$A1.x == gr$A1.y & gr$A2.x == gr$A2.y) | (gr$A1.x == gr$A2.y & gr$A2.x == gr$A1.y)
	diff <- gr[!index,]
	diff$A1.x[diff$A1.x == "C"] <- "g"
	diff$A1.x[diff$A1.x == "G"] <- "c"
	diff$A1.x[diff$A1.x == "T"] <- "a"
	diff$A1.x[diff$A1.x == "A"] <- "t"
	diff$A2.x[diff$A2.x == "C"] <- "g"
	diff$A2.x[diff$A2.x == "G"] <- "c"
	diff$A2.x[diff$A2.x == "T"] <- "a"
	diff$A2.x[diff$A2.x == "A"] <- "t"
	diff$A1.x <- toupper(diff$A1.x)
	diff$A2.x <- toupper(diff$A2.x)

	index2 <- (diff$A1.x == diff$A1.y & diff$A2.x == diff$A2.y) | (diff$A1.x == diff$A2.y & diff$A2.x == diff$A1.y)

	# Number that match initially
	message("SNPs that match: ", sum(index))
	message("SNPs that match after flipping: ", sum(index2))
	message("SNPs that never match: ", sum(!index2))

	prop <- 1 - sum(index2) / sum(index)
	message("Proportion on forward strand: ", prop)

	return(ifelse(prop > threshold, 1, 2))
}



#' Generate VCF header based on genome build
#'
#' Only possible to use b37 or b38 at the moment. Obtained chromosome lengths obtained frmo fai files here http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/
#'
#' @param build='b37' Must be 'b37' or 'b38'
#'
#' @export
#' @return Vector of strings to be used in vcfR package
vcf_header <- function(build='b37')
{
	stopifnot(build %in% c("b37", "b38"))
	info <- c('##fileformat=VCFv4.2',
		'##INFO=<ID=B,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele(s)">',
		'##INFO=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">',
		'##INFO=<ID=PVAL,Number=A,Type=Float,Description="P-value for effect estimate">',
		'##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">',
		'##INFO=<ID=N1,Number=A,Type=Float,Description="Number of cases. 0 if continuous trait">',
		'##INFO=<ID=N0,Number=A,Type=Float,Description="Number of controls. Total sample size if continuous trait">')

	# From fai files here http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/
	chr <- c(1:22, "X", "Y", "MT")
	b37 <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16569)
	b38 <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569)

	contig <- paste0("##contig=<ID=", chr,",length=", get(build),",assembly=",build,">")
	return(c(info, contig))
}


#' Generate VCF format from GWAS summary dataset
#'
#' Takes relevant vectors and produces a vcfR object. This can then be written with \code{write_vcf}. 
#'
#' @param CHROM Vector of chromosome names
#' @param POS Vector of chromosome positions
#' @param ID Vector of SNP IDs
#' @param REF Vector of non-effect alleles - known also as reference allele
#' @param ALT Vector of effect alleles - known also as alternative allele
#' @param QUAL Vector of quality values. Can be NAs.
#' @param FILTER Vector of filter values. Can be NAs.
#' @param B Vector of effects for each SNP, relative to the ALT allele
#' @param SE Vector of standard errors
#' @param PVAL Vector
#' @param AF Vector of allele frequencies for the alternative allele
#' @param N1 Vector of number of cases for each SNP. Vector of 0s if continuous trait
#' @param N0 Vector of number of controls for each SNP. Total sample size if continuous trait
#' @param build="b37" Used in CHROM and POS. This argument used to generate the header. Must be 'b37' or 'b38'.
#'
#' @export
#' @return vcfR object (with empty gt slot)
make_vcf <- function(CHROM, POS, ID, REF, ALT, QUAL, FILTER, B, SE, PVAL, AF, N1, N0, build="b37")
{
	requireNamespace("dplyr", quietly=TRUE)
	stopifnot(all(!is.na(CHROM)))
	stopifnot(all(!is.na(POS)))
	stopifnot(all(!is.na(REF)))
	stopifnot(all(!is.na(ALT)))
	stopifnot(is.numeric(CHROM))
	stopifnot(is.numeric(POS))
	stopifnot(length(POS) == length(CHROM))
	stopifnot(length(ID) == length(CHROM))
	stopifnot(length(REF) == length(CHROM))
	stopifnot(length(ALT) == length(CHROM))
	stopifnot(length(QUAL) == length(CHROM))
	stopifnot(length(FILTER) == length(CHROM))
	stopifnot(length(B) == length(CHROM))
	stopifnot(length(SE) == length(CHROM))
	stopifnot(length(PVAL) == length(CHROM))
	stopifnot(length(AF) == length(CHROM))
	stopifnot(length(N1) == length(CHROM))
	stopifnot(length(N0) == length(CHROM))

	fixed <- dplyr::data_frame(CHROM, POS, ID, REF, ALT, QUAL, FILTER)

	info <- list(B=B, SE=SE, AF=AF, PVAL=PVAL, N1=N1, N0=N0)
	for(i in names(info))
	{
		x <- as.character(info[[i]])
		x[is.na(x)] <- "."
		info[[i]] <- paste0(i, "=", x)
	}
	fixed$INFO <- paste(
		info$B, info$SE, info$PVAL, info$AF, info$N1, info$N0, sep=";"
	)
	fixed <- dplyr::arrange(fixed, CHROM, POS)
	fixed$CHROM <- as.character(fixed$CHROM)
	fixed$POS <- as.character(fixed$POS)
	fixed$ID <- as.character(fixed$ID)
	fixed$REF <- as.character(fixed$REF)
	fixed$ALT <- as.character(fixed$ALT)
	fixed$QUAL <- as.character(fixed$QUAL)
	fixed$FILTER <- as.character(fixed$FILTER)
	fixed[is.na(fixed)] <- "."
	
	vcf <- new(Class="vcfR")
	vcf@meta <- vcf_header(build)
	vcf@fix <- as.matrix(fixed)
	vcf@gt <- matrix("a", nrow=0, ncol=0)
	return(vcf)
}


#' Check if bcftools is available on the path
#'
#' @export
#' @return Status code from running `bcftools --help`. Expect 0 if bcftools available
check_bcftools <- function()
{
	o <- system("bcftools --help", ignore.stdout=TRUE)
	if(o == 0)
	{
		message("bcftools is available")
	} else {
		message("bcftools not found. Please install https://samtools.github.io/bcftools/")
	}
	return(o)
}


#' Write vcf file of summary GWAS data
#'
#' Writes vcf file and indexes. Note, this is only for datasets that have no genotype data. For full implementation of vcf file writing see \code{vcfR::write.vcf}. But for GWAS data it is better to use write_vcf instead of vcfR::write.vcf because the latter creates an ordinary gzip file, and we want to make a bgzip file or bcf file.
#'
#' @param vcf Output from \code{make_vcf}
#' @param filename Must have suffix of .vcf, .vcf.gz or .bcf. If the latter two then index is also created.
#'
#' @export
#' @return NULL
write_vcf <- function(vcf, filename)
{
	if(grepl("vcf$", filename))
	{
		root <- filename
		fext <- "vcf"
	} else if(grepl("vcf\\.gz$", filename)) {
		stopifnot(check_bcftools() == 0)
		root <- gsub("vcf.gz$", "vcf", filename)
		fext <- "vcf.gz"
	} else if(grepl("bcf$", filename)) {
		stopifnot(check_bcftools() == 0)
		root <- gsub("bcf$", "vcf", filename)
		fext <- "bcf"
	} else {
		stop("filename suffix must be vcf, vcf.gz or bcf")
	}

	con <- file(root)
	writeLines(con=con,
		c(vcf@meta, paste0("#", paste(colnames(vcf@fix), collapse="\t"))))
	close(con)
	write.table(vcf@fix, file=root, append=TRUE, row=FALSE, col=FALSE, qu=FALSE, sep="\t")
	if(fext == "vcf.gz")
	{
		cmd <- paste0("bcftools view ", root, " -Oz -o ", filename)
		system(cmd)
		cmd <- paste0("bcftools index ", filename)
		system(cmd)
	}
	if(fext == "bcf")
	{
		cmd <- paste0("bcftools view ", root, " -Ob -o ", filename)
		system(cmd)
		unlink(root)
		cmd <- paste0("bcftools index ", filename)
		system(cmd)
	}
}


