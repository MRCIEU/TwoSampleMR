# Read exposure data

Reads in exposure data. Checks and organises columns for use with MR or
enrichment tests. Infers p-values when possible from beta and se.

## Usage

``` r
read_exposure_data(
  filename,
  clump = FALSE,
  sep = " ",
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
  log_pval = FALSE,
  chr_col = "chr",
  pos_col = "pos",
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  pop = "EUR",
  bfile = NULL,
  plink_bin = NULL
)
```

## Arguments

- filename:

  Filename. Must have header with at least SNP column present.

- clump:

  Whether to perform LD clumping with
  [`clump_data()`](https://mrcieu.github.io/TwoSampleMR/reference/clump_data.md)
  on the exposure data. The default is `FALSE`.

- sep:

  Specify delimiter in file. The default is a space, i.e. `" "`.

- phenotype_col:

  Optional column name for the column with phenotype name corresponding
  the the SNP. If not present then will be created with the value
  "Outcome". The default is `"Phenotype"`.

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

  Required for MR. Name of column with effect allele. Must be "A", "C",
  "T" or "G". The default is `"effect_allele"`.

- other_allele_col:

  Required for MR. Name of column with non effect allele. Must be "A",
  "C", "T" or "G". The default is `"other_allele"`.

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

  Optional column name to give the dataset an ID. Will be generated
  automatically if not provided for every trait / unit combination. The
  default is `"id"`.

- min_pval:

  Minimum allowed p-value. The default is `1e-200`.

- log_pval:

  The p-value is -log10(P). The default is `FALSE`.

- chr_col:

  Optional column name for chromosome. Default is `"chr"`.

- pos_col:

  Optional column name for genetic position Default is `"pos"`.

- clump_kb:

  Clumping window, default is `10000`.

- clump_r2:

  Clumping r2 cutoff. Note that this default value has recently changed
  from `0.01` to `0.001`.

- clump_p1:

  Clumping sig level for index SNPs, default is `1`.

- pop:

  Super-population to use as reference panel. Default = `"EUR"`. Options
  are `"EUR"`, `"SAS"`, `"EAS"`, `"AFR"`, `"AMR"`. `'legacy'` also
  available - which is a previously used version of the EUR panel with a
  slightly different set of markers

- bfile:

  If this is provided then will use the API. Default = `NULL`

- plink_bin:

  If `NULL` and `bfile` is not `NULL` then will detect packaged plink
  binary for specific OS. Otherwise specify path to plink binary.
  Default = `NULL`

## Value

data frame
