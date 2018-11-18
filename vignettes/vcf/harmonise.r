#!/usr/bin/env Rscript

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(argparse)
library(magrittr)

parser <- ArgumentParser()
parser$add_argument('--ref', required=TRUE)
parser$add_argument('--gwas', required=FALSE)
parser$add_argument('--out', required=TRUE)

args <- parser$parse_args()



# 1. Read in the GWAS
# 2. Read in the reference data
# 3. Harmonise GWAS against the reference
# 4. Write to VCF


# Read the GWAS
# Just assuming the format used for uploading to elastic
gwas <- data.table::fread(paste0("gunzip -c ", args[["gwas"]]))
names(gwas) <- c("snp_col", "ea_col", "oa_col", "eaf_col", "beta_col", "se_col", "pval_col", "ncontrol_col")
# This is a continuous GWAS so no ncase column
gwas$ncase_col <- NA


# Read the reference
ref <- data.table::fread(paste0("gunzip -c ", args[["ref"]]))
stopifnot(c("CHROM", "ID", "REF", "ALT", "MAF", "POS") %in% names(ref))

# For simplicity just keeping SNP Ids that are in common
ref <- subset(ref, ID %in% gwas$snp_col)

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
	eaf_col="MAF"
)

b <- TwoSampleMR::format_data(gwas, type="outcome", 
	snp_col="snp_col",
	beta_col="beta_col",
	se_col="se_col",
	effect_allele_col="ea_col",
	other_allele_col="oa_col",
	eaf_col="eaf_col",
	ncase_col="ncase_col",
	ncontrol_col="ncontrol_col",
	pval_col="pval_col"
)

# Is the gwas on the forward strand?
action <- is_forward_strand(gwas$snp_col, gwas$ea_col, gwas$oa_col, ref$ID, ref$ALT, ref$REF)

# Harmonise the gwas according to the reference panel
ab <- TwoSampleMR::harmonise_data(a, b, action=action)

gwas_h <- ab %$% 
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
	dplyr::inner_join(subset(ref, select=c(ID,REF,ALT,CHROM,POS,MAF)), by=c("ID", "REF", "ALT"))

save(gwas_h, file=args[["out"]])

