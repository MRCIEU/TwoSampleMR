#!/bin/bash
set -euo pipefail

module load apps/gatk-4.0.8.1

#validate input VCF
gatk ValidateVariants \
-R ~/db/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$1" \
--dbsnp ~/db/gatk/2.8/b37/dbsnp_138.b37.vcf
