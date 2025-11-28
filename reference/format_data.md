# Read exposure or outcome data

Reads in exposure data. Checks and organises columns for use with MR or
enrichment tests. Infers p-values when possible from beta and se.

## Usage

``` r
format_data(
  dat,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  z_col = "z",
  info_col = "info",
  chr_col = "chr",
  pos_col = "pos",
  log_pval = FALSE
)
```

## Arguments

- dat:

  Data frame. Must have header with at least SNP column present.

- type:

  Is this the exposure or the outcome data that is being read in? The
  default is `"exposure"`.

- snps:

  SNPs to extract. If NULL then doesn't extract any and keeps all. The
  default is `NULL`.

- header:

  The default is `TRUE`.

- phenotype_col:

  Optional column name for the column with phenotype name corresponding
  the the SNP. If not present then will be created with the value
  `"Outcome"`. The default is `"Phenotype"`.

- snp_col:

  Required name of column with SNP rs IDs. The default is `"SNP"`.

- beta_col:

  Required for MR. Name of column with effect sizes. The default is
  `"beta"`.

- se_col:

  Required for MR. Name of column with standard errors. The default is
  `"se"`.

- eaf_col:

  Required for MR. Name of column with effect allele frequency. The
  default is `"eaf"`.

- effect_allele_col:

  Required for MR. Name of column with effect allele. Must contain only
  the characters "A", "C", "T" or "G". The default is `"effect_allele"`.

- other_allele_col:

  Required for MR. Name of column with non effect allele. Must contain
  only the characters "A", "C", "T" or "G". The default is
  `"other_allele"`.

- pval_col:

  Required for enrichment tests. Name of column with p-value. The
  default is `"pval"`.

- units_col:

  Optional column name for units. The default is `"units"`.

- ncase_col:

  Optional column name for number of cases. The default is `"ncase"`.

- ncontrol_col:

  Optional column name for number of controls. The default is
  `"ncontrol"`.

- samplesize_col:

  Optional column name for sample size. The default is `"samplesize"`.

- gene_col:

  Optional column name for gene name. The default is `"gene"`.

- id_col:

  The default is `"id"`.

- min_pval:

  Minimum allowed p-value. The default is `1e-200`.

- z_col:

  The default is `"z"`.

- info_col:

  The default is `"info_col"`.

- chr_col:

  The default is `"chr_col"`.

- pos_col:

  The default is `"pos"`.

- log_pval:

  The pval is -log10(P). The default is `FALSE`.

## Value

data frame
