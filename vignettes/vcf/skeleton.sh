#!/bin/bash



# BCF file with:
# - correct REF and ALT alleles
# - AF field in INFO that is the ALT allele frequency
ref="../../../reference/1000g/1kg_v3_nomult.bcf"


# GWAS file:
gwas="~/mr-eve/gwas-instrument-subsets/studies/2/elastic.gz"


# 0. Clean GWAS



# STILL TO DO



# 1. check for name merges against sqlite
# This file has rs ID merges
# https://www.ncbi.nlm.nih.gov/projects/SNP/snp_db_table_description.cgi?t=RsMergeArch
# Step 1 is to update any rs IDs in the GWAS based on this file



# STILL TO DO



# 2. Convert any chr:pos SNPs in the GWAS to rs IDs
# Do this by extracting variants that are missing rs IDs in GWAS
# and find the rs ID in the reference
# and update the GWAS file



# STILL TO DO



# 3. Get subset of reference in tab format


gunzip -c elastic.gz | cut -f 1 > snplist.txt
wc -l snplist.txt

time bcftools view -i'ID=@temp' $ref | bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' | sed '1 i\
CHROM\tPOS\tID\tREF\tALT\tAF
' | gzip -c > ref_extract.txt.gz



# 4. Harmonise the GWAS file against the reference
# This needs to retain indels so use the TwoSampleMR::harmonise_data
# function. It will: 
# - switch effect alleles
# - handle sequence coded indels
# - convert D/I indels to sequence coding (as in the reference)
# - check for forward strand and flip if necessary
# - tries to harmonise with only effect allele if other allele not available
# 4b. Write out to bcf format
# After harmonising can use the TwoSampleMR::write_vcf function
# It will create file based on extension and index.



Rscript harmonise_against_ref.r \
--ref-file ref_extract.txt.gz \
--ref-build b37 \
--gwas-file $gwas \
--gwas-header FALSE \
--gwas-snp 1 \
--gwas-ref 3 \
--gwas-alt 2 \
--gwas-af 4 \
--gwas-beta 5 \
--gwas-se 6 \
--gwas-pval 7 \
--gwas-n0 8 \
--gwas-n1 NA \
--out harmonised.bcf



# 5. Create report and json document of harmonising stats
# This could be included above

