#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(TwoSampleMR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vcfR))
library(methods)
library(utils)

# create parser object
parser <- ArgumentParser()

parser$add_argument('--ref-file', required=TRUE)
parser$add_argument('--ref-build', required=TRUE, default="b37")
parser$add_argument('--gwas-file', required=TRUE)
parser$add_argument('--gwas-header', required=TRUE, type="logical", default=FALSE)
parser$add_argument('--gwas-snp', type="integer", required=TRUE)
parser$add_argument('--gwas-ref', type="integer", required=FALSE)
parser$add_argument('--gwas-alt', type="integer", required=TRUE)
parser$add_argument('--gwas-af', type="integer", required=FALSE)
parser$add_argument('--gwas-beta', type="integer", required=FALSE)
parser$add_argument('--gwas-se', type="integer", required=FALSE)
parser$add_argument('--gwas-pval', type="integer", required=FALSE)
parser$add_argument('--gwas-n0', type="integer", required=FALSE)
parser$add_argument('--gwas-n1', type="integer", required=FALSE)
parser$add_argument('--out', required=TRUE)

args <- parser$parse_args()

print(args)

read_dat <- function(filename, type, header, snp, ref, alt, af, beta, se, pval, n0, n1)
{
	if(grepl("gz$", filename))
	{
		dat <- data.table::fread(paste0("gunzip -c ", filename), header=header)
	} else {
		dat <- data.table::fread(filename, header=header)
	}
	nc <- ncol(dat)
	if(snp == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		snp <- ncol(dat)
	}
	if(ref == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		ref <- ncol(dat)
	}
	if(alt == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		alt <- ncol(dat)
	}
	if(af == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		af <- ncol(dat)
	}
	if(beta == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		beta <- ncol(dat)
	}
	if(se == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		se <- ncol(dat)
	}
	if(pval == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		pval <- ncol(dat)
	}
	if(n0 == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		n0 <- ncol(dat)
	}
	if(n1 == 0)
	{
		dat[[paste0("V", ncol(dat)+1)]] <- rep(NA, nrow(dat))
		n1 <- ncol(dat)
	}

	o <- format_data(
		dat, 
		type=type, 
		phenotype_col=type,
		snp_col=names(dat)[snp],
		beta_col=names(dat)[beta],
		se_col=names(dat)[se],
		effect_allele_col=names(dat)[alt],
		other_allele_col=names(dat)[ref],
		eaf_col=names(dat)[af],
		pval_col=names(dat)[pval],
		ncase_col=names(dat)[n1],
		ncontrol_col=names(dat)[n0]
	)
	return(o)
}



# Read in gwas data
gwas <- read_dat(
	args[["gwas_file"]],
	type="outcome",
	header=args[["gwas_header"]],
	snp=args[["gwas_snp"]],
	ref=args[["gwas_ref"]],
	alt=args[["gwas_alt"]],
	af=args[["gwas_af"]],
	beta=args[["gwas_beta"]],
	se=args[["gwas_se"]],
	pval=args[["gwas_pval"]],
	n0=args[["gwas_n0"]],
	n1=args[["gwas_n1"]]
)


# Read in ref

ref <- data.table::fread(paste0("gunzip -c ", args[["ref_file"]]))
stopifnot(all(c("CHROM", "ID", "REF", "ALT", "AF", "POS") %in% names(ref)))

# For simplicity just keeping SNP Ids that are in common
ref <- subset(ref, ID %in% gwas$SNP)

# Put in some dummy variables for the reference for harmonising
ref$beta <- 1
ref$se <- 0.1
ref$pval <- 0.1
a <- TwoSampleMR::format_data(
        ref,
        type="exposure",
        snp_col="ID",
        effect_allele_col="ALT",
        other_allele_col="REF",
        eaf_col="AF"
)

# Check strand
action <- TwoSampleMR::is_forward_strand(gwas$SNP, gwas$effect_allele.outcome, gwas$other_allele.outcome, ref$ID, ref$ALT, ref$REF, threshold=0.9)

# Harmonise
dat <- TwoSampleMR::harmonise_data(a, gwas, action)


gwas_h <- dat %$%
	dplyr::data_frame(
		ID=SNP,
		ALT=effect_allele.exposure,
		REF=other_allele.exposure,
		BETA=beta.outcome,
		SE=se.outcome,
		PVALUE=pval.outcome,
		AF=eaf.outcome,
		N=samplesize.outcome,
		NCASE=ncase.outcome,
		NCONTROL=ncontrol.outcome) %>%
	dplyr::inner_join(subset(ref, select=c(ID,REF,ALT,CHROM,POS)), by=c("ID", "REF", "ALT"))


# Create vcf format
vcf <- TwoSampleMR::make_vcf(
                ID = gwas_h$ID,
                ALT = gwas_h$ALT,
                REF = gwas_h$REF,
                B = gwas_h$BETA,
                SE = gwas_h$SE,
                PVAL = gwas_h$PVALUE,
                N0 = gwas_h$NCONTROL,
                N1 = gwas_h$NCASE,
                CHROM = gwas_h$CHROM,
                POS = gwas_h$POS,
                AF = gwas_h$AF,
                QUAL = rep(NA, nrow(gwas_h)),
                FILTER = rep("PASS", nrow(gwas_h)),
                build = args[["ref_build"]]
        )

# Write vcf
TwoSampleMR::write_vcf(vcf, args[["out"]])


