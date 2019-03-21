#!/bin/bash

# Download 1000 genomes

wget ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
wget ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi
wget ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
wget ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi
wget ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz
wget ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz.tbi


for f in *vcf.gz
do
	echo $f
	i=`echo $f | cut -d "." -f 2`
	echo $i
	bcftools view -G -Ob $f > 1kg_v3_$i.bcf
	bcftools index 1kg_v3_$i.bcf
done

bcftools concat 1kg_v3_*.bcf -Ob > 1kg_v3.bcf

# Multi allelic SNPs are problematic for harmonisation because we would have to update TwoSampleMR functions to allow this.
# Most GWASs drop multi-allelic SNPs so we can probably try to ignore, though remains to be seen how much we lose due to this strategy
# Split and keep first biallelic version of each variant only
# Is there a better way to handle multi-allelic SNPs?
bcftools norm -m- 1kg_v3.bcf | bcftools norm -d all -Ob > 1kg_v3_nomult.bcf

# Note we can also use dbSNP which has ALT allele frequencies
# ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20180418.vcf.gz

# Or dbSNP filter of only common variants
# ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz
# Though this is slightly prooblematic because they have a CAF column instead of a standard AF column, 
# so handling multiallelic SNPs won't automatically handle this column
