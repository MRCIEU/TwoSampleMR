# Attempt to perform MVMR using local data

Allows you to read in summary data from text files to format the
multivariable exposure dataset.

## Usage

``` r
mv_extract_exposures_local(
  filenames_exposure,
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
  pval_threshold = 5e-08,
  plink_bin = NULL,
  bfile = NULL,
  clump_r2 = 0.001,
  clump_kb = 10000,
  pop = "EUR",
  harmonise_strictness = 2
)
```

## Arguments

- filenames_exposure:

  Filenames for each exposure dataset. Must have header with at least
  SNP column present. Following arguments are used for determining how
  to read the filename and clumping etc.

- sep:

  Specify delimiter in file. The default is space, i.e. `sep=" "`. If
  length is 1 it will use the same `sep` value for each exposure
  dataset. You can provide a vector of values, one for each exposure
  dataset, if the values are different across datasets. The same applies
  to all dataset-formatting options listed below.

- phenotype_col:

  Optional column name for the column with phenotype name corresponding
  the the SNP. If not present then will be created with the value
  `"Outcome"`. Default is `"Phenotype"`.

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

  The pval is -log10(P). The default is `FALSE`.

- pval_threshold:

  Default=`5e-8` for clumping

- plink_bin:

  If `NULL` and `bfile` is not `NULL` then will detect packaged plink
  binary for specific OS. Otherwise specify path to plink binary.
  Default = `NULL`

- bfile:

  If this is provided then will use the API. Default = `NULL`

- clump_r2:

  Default=`0.001` for clumping

- clump_kb:

  Default=`10000` for clumping

- pop:

  Which 1000 genomes super population to use for clumping when using the
  server

- harmonise_strictness:

  See action argument in
  [`harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.md).
  Default=`2`

## Value

List

## Details

Note that you can provide an array of column names for each column,
which is of length `filenames_exposure`
