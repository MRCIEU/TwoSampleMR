# Exposure data

``` r
library(TwoSampleMR)
```

## Introduction

A data frame of the instruments for an exposure is required. Each line
has the information for one variant for one exposure. The minimum
information required for MR analysis is the following:

- `SNP` - rs ID
- `beta` - The effect size. If the trait is binary then log(OR) should
  be used
- `se` - The standard error of the effect size
- `effect_allele` - The allele of the SNP which has the effect marked in
  `beta`

Other information that is useful for MR can also be provided:

- `other_allele` - The non-effect allele
- `eaf` - The effect allele frequency
- `Phenotype` - The name of the phenotype for which the SNP has an
  effect

You can also provide the following extra information:

- `chr` - Physical position of variant (chromosome)
- `position` - Physical position of variant (position)
- `samplesize` - Sample size for estimating the effect size
- `ncase` - Number of cases
- `ncontrol` - Number of controls
- `pval` - The P-value for the SNP’s association with the exposure
- `units` - The units in which the effects are presented
- `gene` - The gene or other annotation for the the SNP

## Reading in from a file

The data can be read in from a text file using the `read_exposure_data`
function. The file must have a header with column names corresponding to
the columns described above.

### Example 1: The default column names are used

An example of a text file with the default column names is provided as
part of the package, the first few rows look like this:

    Phenotype SNP beta se effect_allele other_allele eaf pval units gene samplesize
    BMI rs10767664 0.19 0.0306122448979592 A T 0.78 5e-26 kg/m2 BDNF 225238
    BMI rs13078807 0.1 0.0204081632653061 G A 0.2 4e-11 kg/m2 CADM2 221431
    BMI rs1514175 0.07 0.0204081632653061 A G 0.43 8e-14 kg/m2 TNNI3K 207641
    BMI rs1558902 0.39 0.0204081632653061 A T 0.42 5e-120 kg/m2 FTO 222476
    BMI rs10968576 0.11 0.0204081632653061 G A 0.31 3e-13 kg/m2 LRRN6C 247166
    BMI rs2241423 0.13 0.0204081632653061 G A 0.78 1e-18 kg/m2 LBXCOR1 227886

The exact path to the file will be different on everyone’s computer, but
it can be located like this:

``` r
bmi_file <- system.file("extdata", "bmi.txt", package = "TwoSampleMR")
```

You can read the data in like this:

``` r
bmi_exp_dat <- read_exposure_data(bmi_file)
head(bmi_exp_dat)
#>          SNP beta.exposure se.exposure effect_allele.exposure
#> 1 rs10767664          0.19  0.03061224                      A
#> 2 rs13078807          0.10  0.02040816                      G
#> 3  rs1514175          0.07  0.02040816                      A
#> 4  rs1558902          0.39  0.02040816                      A
#> 5 rs10968576          0.11  0.02040816                      G
#> 6  rs2241423          0.13  0.02040816                      G
#>   other_allele.exposure eaf.exposure pval.exposure units.exposure gene.exposure
#> 1                     T         0.78         5e-26          kg/m2          BDNF
#> 2                     A         0.20         4e-11          kg/m2         CADM2
#> 3                     G         0.43         8e-14          kg/m2        TNNI3K
#> 4                     T         0.42        5e-120          kg/m2           FTO
#> 5                     A         0.31         3e-13          kg/m2        LRRN6C
#> 6                     A         0.78         1e-18          kg/m2       LBXCOR1
#>   samplesize.exposure exposure mr_keep.exposure pval_origin.exposure
#> 1              225238      BMI             TRUE             reported
#> 2              221431      BMI             TRUE             reported
#> 3              207641      BMI             TRUE             reported
#> 4              222476      BMI             TRUE             reported
#> 5              247166      BMI             TRUE             reported
#> 6              227886      BMI             TRUE             reported
#>   units.exposure_dat id.exposure data_source.exposure
#> 1              kg/m2      ImbABK             textfile
#> 2              kg/m2      ImbABK             textfile
#> 3              kg/m2      ImbABK             textfile
#> 4              kg/m2      ImbABK             textfile
#> 5              kg/m2      ImbABK             textfile
#> 6              kg/m2      ImbABK             textfile
```

The output from this function is a new data frame with standardised
column names:

- `SNP`
- `exposure`
- `beta.exposure`
- `se.exposure`
- `effect_allele.exposure`
- `other_allele.exposure`
- `eaf.exposure`
- `mr_keep.exposure`
- `pval.exposure`
- `pval_origin.exposure`
- `id.exposure`
- `data_source.exposure`
- `units.exposure`
- `gene.exposure`
- `samplesize.exposure`

The function attempts to match the columns to the ones it expects. It
also checks that the data type is as expected.

If the required data for MR to be performed is not present (SNP name,
effect size, standard error, effect allele) for a particular SNP, then
the column `mr_keep.exposure` will be `FALSE`.

### Example 2: The text file has non-default column names

If the text file does not have default column names, this can still be
read in as follows. Here are the first few rows of an example:

    rsid,effect,SE,a1,a2,a1_freq,p-value,Units,Gene,n
    rs10767664,0.19,0.030612245,A,T,0.78,5.00E-26,kg/m2,BDNF,225238
    rs13078807,0.1,0.020408163,G,A,0.2,4.00E-11,kg/m2,CADM2,221431
    rs1514175,0.07,0.020408163,A,G,0.43,8.00E-14,kg/m2,TNNI3K,207641
    rs1558902,0.39,0.020408163,A,T,0.42,5.00E-120,kg/m2,FTO,222476

Note that this is a CSV file, with commas separating fields. The file is
located here:

``` r
bmi2_file <- system.file("extdata/bmi.csv", package = "TwoSampleMR")
```

To read in this data:

``` r
bmi_exp_dat <- read_exposure_data(
  filename = bmi2_file,
  sep = ",",
  snp_col = "rsid",
  beta_col = "effect",
  se_col = "SE",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  eaf_col = "a1_freq",
  pval_col = "p-value",
  units_col = "Units",
  gene_col = "Gene",
  samplesize_col = "n"
)
#> No phenotype name specified, defaulting to 'exposure'.
head(bmi_exp_dat)
#>          SNP beta.exposure se.exposure effect_allele.exposure
#> 1 rs10767664          0.19  0.03061224                      A
#> 2 rs13078807          0.10  0.02040816                      G
#> 3  rs1514175          0.07  0.02040816                      A
#> 4  rs1558902          0.39  0.02040816                      A
#> 5 rs10968576          0.11  0.02040816                      G
#> 6  rs2241423          0.13  0.02040816                      G
#>   other_allele.exposure eaf.exposure pval.exposure units.exposure gene.exposure
#> 1                     T         0.78         5e-26          kg/m2          BDNF
#> 2                     A         0.20         4e-11          kg/m2         CADM2
#> 3                     G         0.43         8e-14          kg/m2        TNNI3K
#> 4                     T         0.42        5e-120          kg/m2           FTO
#> 5                     A         0.31         3e-13          kg/m2        LRRN6C
#> 6                     A         0.78         1e-18          kg/m2       LBXCOR1
#>   samplesize.exposure exposure mr_keep.exposure pval_origin.exposure
#> 1              225238 exposure             TRUE             reported
#> 2              221431 exposure             TRUE             reported
#> 3              207641 exposure             TRUE             reported
#> 4              222476 exposure             TRUE             reported
#> 5              247166 exposure             TRUE             reported
#> 6              227886 exposure             TRUE             reported
#>   units.exposure_dat id.exposure data_source.exposure
#> 1              kg/m2      Ku3B84             textfile
#> 2              kg/m2      Ku3B84             textfile
#> 3              kg/m2      Ku3B84             textfile
#> 4              kg/m2      Ku3B84             textfile
#> 5              kg/m2      Ku3B84             textfile
#> 6              kg/m2      Ku3B84             textfile
```

If the `Phenotype` column is not provided (as is the case in this
example) then it will assume that the phenotype’s name is simply
“exposure”. This is entered in the `exposure` column. It can be renamed
manually:

``` r
bmi_exp_dat$exposure <- "BMI"
```

## Using an existing data frame

If the data already exists as a data frame in R then it can be converted
into the correct format using the
[`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md)
function. For example, here is some randomly created data:

``` r
random_df <- data.frame(
  SNP = c("rs1", "rs2"),
  beta = c(1, 2),
  se = c(1, 2),
  effect_allele = c("A", "T")
)
random_df
#>   SNP beta se effect_allele
#> 1 rs1    1  1             A
#> 2 rs2    2  2             T
```

This can be formatted like so:

``` r
random_exp_dat <- format_data(random_df, type = "exposure")
#> No phenotype name specified, defaulting to 'exposure'.
#> Warning in format_data(random_df, type = "exposure"): The following columns are not present but are helpful for harmonisation
#> other_alleleeaf
#> Inferring p-values
random_exp_dat
#>   SNP beta.exposure se.exposure effect_allele.exposure exposure
#> 1 rs1             1           1                      A exposure
#> 2 rs2             2           2                      T exposure
#>   mr_keep.exposure pval.exposure pval_origin.exposure id.exposure
#> 1             TRUE     0.3173105             inferred      4nec1Z
#> 2             TRUE     0.3173105             inferred      4nec1Z
#>   other_allele.exposure eaf.exposure
#> 1                    NA           NA
#> 2                    NA           NA
```

## Obtaining instruments from existing catalogues

A number of sources of instruments have already been curated and are
available for use. They are provided as data objects in the
`MRInstruments` package. To install:

``` r
remotes::install_github("MRCIEU/MRInstruments")
```

This package contains a number of data.frames, each of which is a
repository of SNP-trait associations. How to access the data frames is
detailed below:

### GWAS catalog

The NHGRI-EBI GWAS catalog contains a catalog of significant
associations obtained from GWASs. This version of the data is filtered
and harmonised to contain associations that have the required data to
perform MR, to ensure that the units used to report effect sizes from a
particular study are all the same, and other data cleaning operations.

To use the GWAS catalog:

``` r
library(MRInstruments)
data(gwas_catalog)
head(gwas_catalog)
#>                                                 Phenotype_simple
#> 1                           Eosinophil percentage of white cells
#> 2                                              Eosinophil counts
#> 3 Medication use (agents acting on the renin-angiotensin system)
#> 4                                       Post bronchodilator FEV1
#> 5                         DNA methylation variation (age effect)
#> 6                                         Ankylosing spondylitis
#>                                                MAPPED_TRAIT_EFO
#> 1                           eosinophil percentage of leukocytes
#> 2                                              eosinophil count
#> 3 Agents acting on the renin-angiotensin system use measurement
#> 4          forced expiratory volume, response to bronchodilator
#> 5                                               DNA methylation
#> 6                                        ankylosing spondylitis
#>                                                              MAPPED_TRAIT_EFO_URI
#> 1                                            http://www.ebi.ac.uk/efo/EFO_0007991
#> 2                                            http://www.ebi.ac.uk/efo/EFO_0004842
#> 3                                            http://www.ebi.ac.uk/efo/EFO_0009931
#> 4 http://www.ebi.ac.uk/efo/EFO_0004314, http://purl.obolibrary.org/obo/GO_0097366
#> 5                                       http://purl.obolibrary.org/obo/GO_0006306
#> 6                                            http://www.ebi.ac.uk/efo/EFO_0003898
#>                                                                                                                                                Initial_sample_description
#> 1                                                                                                                                   172,378 European ancestry individuals
#> 2                                                                                                                                   172,275 European ancestry individuals
#> 3                                                                                                      62,752 European ancestry cases, 174,778 European ancestry controls
#> 4 10,094 European ancestry current and former smoker individuals, 3,260 African American current and former smoker individuals, 178 current and former smoker individuals
#> 5                                                                                                                                                   Up to 954 individuals
#> 6                                                    921 Turkish ancestry cases, 907 Turkish ancestry controls, 422 Iranian ancestry cases, 754 Iranian ancestry controls
#>   Replication_sample_description STUDY.ACCESSION
#> 1                           <NA>      GCST004600
#> 2                           <NA>      GCST004606
#> 3                           <NA>      GCST007930
#> 4                           <NA>      GCST003262
#> 5                           <NA>      GCST006660
#> 6                           <NA>      GCST007844
#>                                                        Phenotype Phenotype_info
#> 1                           Eosinophil percentage of white cells               
#> 2                                              Eosinophil counts               
#> 3 Medication use (agents acting on the renin-angiotensin system)               
#> 4                                       Post bronchodilator FEV1               
#> 5                         DNA methylation variation (age effect)               
#> 6                                         Ankylosing spondylitis               
#>   PubmedID   Author Year        SNP chr bp_ens_GRCh38   Region       gene
#> 1 27863252 Astle WJ 2016  rs1000005  21      33060745 21q22.11 AP000282.2
#> 2 27863252 Astle WJ 2016  rs1000005  21      33060745 21q22.11 AP000282.2
#> 3 31015401     Wu Y 2019  rs1000010   3      11562645   3p25.3      VGLL4
#> 4 26634245  Lutz SM 2015 rs10000225   4     144312789  4q31.21 Intergenic
#> 5 30348214  Zhang Q 2018 rs10000513   4     160334994   4q32.1         NR
#> 6 30946743     Li Z 2019 rs10000518   4      11502867  4p15.33     HS3ST1
#>               Gene_ens effect_allele other_allele        beta          se  pval
#> 1 AP000282.2,LINC00945             C            G -0.02652552 0.003826531 2e-13
#> 2 AP000282.2,LINC00945             C            G -0.02481715 0.003571429 7e-12
#> 3                                  G            A -0.03724189 0.006377551 6e-09
#> 4           Intergenic             A            T -0.04400000 0.009420188 3e-06
#> 5                   NR          <NA>         <NA>          NA          NA 4e-08
#> 6                                  G            A  0.73396926          NA 6e-06
#>              units      eaf date_added_to_MRBASE
#> 1    unit decrease 0.589400           2019-08-29
#> 2    unit decrease 0.589400           2019-08-29
#> 3    unit decrease 0.351806           2019-08-29
#> 4 NR unit decrease 0.350000           2019-08-29
#> 5             <NA>       NA           2019-08-29
#> 6             <NA>       NA           2019-08-29
```

For example, to obtain instruments for body mass index using the
Speliotes et al 2010 study:

``` r
bmi_gwas <-
  subset(gwas_catalog,
         grepl("Speliotes", Author) &
           Phenotype == "Body mass index")
bmi_exp_dat <- format_data(bmi_gwas)
```

### Metabolites

Independent top hits from GWASs on 121 metabolites in whole blood are
stored in the `metab_qtls` data object. Use
[`?metab_qtls`](https://rdrr.io/pkg/MRInstruments/man/metab_qtls.html)
to get more information.

``` r
data(metab_qtls)
head(metab_qtls)
#>   phenotype chromosome  position       SNP effect_allele other_allele      eaf
#> 1     AcAce          8   9181395 rs2169387             G            A 0.870251
#> 2     AcAce         11 116648917  rs964184             C            G 0.857715
#> 3       Ace          6  12042473 rs6933521             C            T 0.120256
#> 4       Ala          2  27730940 rs1260326             C            T 0.638817
#> 5       Ala          2  65220910 rs2160387             C            T 0.403170
#> 6       Ala         12  47201814 rs4554975             G            A 0.644059
#>        beta       se     pval n_studies     n
#> 1  0.085630 0.015451 3.61e-08        11 19257
#> 2 -0.096027 0.014624 6.71e-11        11 19261
#> 3 -0.091667 0.015885 8.10e-09        14 24742
#> 4 -0.104582 0.009940 7.40e-26        13 22569
#> 5 -0.071001 0.009603 1.49e-13        14 24793
#> 6 -0.069135 0.009598 6.12e-13        14 24792
```

For example, to obtain instruments for Alanine:

``` r
ala_exp_dat <- format_metab_qtls(subset(metab_qtls, phenotype == "Ala"))
```

### Proteins

Independent top hits from GWASs on 47 protein levels in whole blood are
stored in the `proteomic_qtls` data object. Use
[`?proteomic_qtls`](https://rdrr.io/pkg/MRInstruments/man/proteomic_qtls.html)
to get more information.

``` r
data(proteomic_qtls)
head(proteomic_qtls)
#>   analyte chr  position        SNP  gene location annotation other_allele
#> 1   CFHR1   1 196698945 rs12144939   CFH      cis   missense            T
#> 2    IL6r   1 154425456 rs12126142  IL6R      cis   missense            A
#> 3   ApoA4  11 116677723  rs1263167 APOA4      cis intergenic            G
#> 4    SELE   9 136149399   rs507666   ABO    trans   intronic            A
#> 5 FetuinA   3 186335941  rs2070633  AHSG      cis   missense            T
#> 6     ACE  17  61566031     rs4343   ACE      cis synonymous            A
#>   effect_allele   eaf   maf      pval   beta         se
#> 1             G 0.643 0.357 8.99e-143 -1.108 0.04355258
#> 2             G 0.608 0.392 1.81e-106  0.850 0.03878364
#> 3             A 0.803 0.197  2.64e-54 -0.919 0.05922332
#> 4             G 0.809 0.191  1.01e-52 -0.882 0.05771545
#> 5             C 0.676 0.324  2.88e-44 -0.629 0.04506925
#> 6             G 0.508 0.492  6.66e-44  0.493 0.03547679
```

For example, to obtain instruments for the ApoH protein:

``` r
apoh_exp_dat <-
  format_proteomic_qtls(subset(proteomic_qtls, analyte == "ApoH"))
```

### Gene expression levels

Independent top hits from GWASs on 32432 gene identifiers and in 44
tissues are available from the GTEX study in `gtex_eqtl`. Use
[`?gtex_eqtl`](https://rdrr.io/pkg/MRInstruments/man/gtex_eqtl.html) to
get more information.

``` r
data(gtex_eqtl)
head(gtex_eqtl)
#>                 tissue     gene_name gene_start         SNP snp_position
#> 1 Adipose Subcutaneous RP4-669L17.10   1:317720   rs2519065     1:787151
#> 2 Adipose Subcutaneous RP11-206L10.1   1:661611  rs11804171     1:723819
#> 3 Adipose Subcutaneous RP11-206L10.3   1:677193 rs149110718     1:759227
#> 4 Adipose Subcutaneous RP11-206L10.2   1:700306 rs148649543     1:752796
#> 5 Adipose Subcutaneous RP11-206L10.9   1:714150  rs12184279     1:717485
#> 6 Adipose Subcutaneous RP11-206L10.8   1:736259  rs10454454     1:754954
#>   effect_allele other_allele      beta        se        pval   n
#> 1             A            G  0.551788 0.0747180 2.14627e-12 298
#> 2             A            T -0.917475 0.1150060 4.99967e-14 298
#> 3             T            C  0.807571 0.1776530 8.44694e-06 298
#> 4             T            C  0.745393 0.0958531 1.82660e-13 298
#> 5             A            C  1.927250 0.2247390 9.55098e-16 298
#> 6             A            G  1.000400 0.1787470 5.61079e-08 298
```

For example, to obtain instruments for the IRAK1BP1 gene expression
levels in subcutaneous adipose tissue:

``` r
irak1bp1_exp_dat <-
  format_gtex_eqtl(subset(
    gtex_eqtl,
    gene_name == "IRAK1BP1" & tissue == "Adipose Subcutaneous"
  ))
#> Warning in format_data(gtex_eqtl_subset, type = type, phenotype_col = type, : The following columns are not present but are helpful for harmonisation
#> eaf
#> Inferring p-values
```

### DNA methylation levels

Independent top hits from GWASs on 0 DNA methylation levels in whole
blood across 5 time points are available from the ARIES study in
`aries_mqtl`. Use
[`?aries_mqtl`](https://rdrr.io/pkg/MRInstruments/man/aries_mqtl.html)
to get more information.

``` r
data(aries_mqtl)
head(aries_mqtl)
#>          SNP timepoint        cpg    beta        pval         se snp_chr
#> 1 esv2656832         1 cg21826606  0.3459 1.60408e-26 0.03265336       1
#> 2 esv2658098         1 cg22681495 -0.6263 1.55765e-66 0.03643240      15
#> 3 esv2660043         1 cg24276624 -0.5772 3.16370e-26 0.05481823      11
#> 4 esv2660043         1 cg11157765 -0.5423 1.33928e-22 0.05583777      11
#> 5 esv2660673         1 cg05832925 -0.5919 2.88011e-50 0.03982467      11
#> 6 esv2660769         1 cg05859533 -0.6224 1.49085e-58 0.03868158      16
#>    snp_pos effect_allele other_allele    eaf   sex   age    units
#> 1 25591901             I            R 0.3974 mixed Birth SD units
#> 2 86057007             D            R 0.2076 mixed Birth SD units
#> 3 69982552             D            R 0.1450 mixed Birth SD units
#> 4 69982552             D            R 0.1450 mixed Birth SD units
#> 5 74024905             D            R 0.1671 mixed Birth SD units
#> 6 57725395             D            R 0.2136 mixed Birth SD units
#>   island_location cpg_chr  cpg_pos    gene gene_location cis_trans
#> 1         N_Shore       1 25593055                             cis
#> 2                      15 86058755  AKAP13          Body       cis
#> 3                      11 69982941    ANO1          Body       cis
#> 4                      11 69982996    ANO1          Body       cis
#> 5         S_Shelf      11 74026371                             cis
#> 6                      16 57727230 CCDC135       TSS1500       cis
```

For example, to obtain instruments for cg25212131 CpG DNA methylation
levels in at birth:

``` r
cg25212131_exp_dat <-
  format_aries_mqtl(subset(aries_mqtl, cpg == "cg25212131" &
                             age == "Birth"))
```

### IEU OpenGWAS database

The IEU OpenGWAS database contains the entire summary statistics for
thousands of GWASs. You can browse them here:
<https://gwas.mrcieu.ac.uk/>

You can use this database to define the instruments for a particular
exposure. You can also use this database to obtain the effects for
constructing polygenic risk scores using different p-value thresholds.

You can check the status of the API:

``` r
ieugwasr::api_status()
```

To obtain a list and details about the available GWASs do the following:

``` r
ao <- available_outcomes()
head(ao)
```

For information about authentication see
<https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication>.

The
[`available_outcomes()`](https://mrcieu.github.io/TwoSampleMR/reference/available_outcomes.md)
function returns a table of all the available studies in the database.
Each study has a unique ID. e.g., You might obtain

``` r
head(subset(ao, select = c(trait, id)))
#>           trait         id
#> 1 Schizophrenia ieu-b-5103
#> 2 Schizophrenia ieu-b-5102
#> 3 Schizophrenia ieu-b-5101
#> 4 Schizophrenia ieu-b-5100
#> 5 Schizophrenia ieu-b-5099
#> 6 Schizophrenia ieu-b-5098
```

To extract instruments for a particular trait using a particular study,
for example to obtain SNPs for body mass index using the Locke et
al. 2015 GIANT study, you specify the study ID as follows:

``` r
bmi2014_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')
```

``` r
str(bmi2014_exp_dat)
#> 'data.frame':    79 obs. of  15 variables:
#>  $ pval.exposure         : num  2.18e-08 4.57e-11 5.06e-14 5.45e-10 1.88e-28 ...
#>  $ samplesize.exposure   : num  339152 339065 313621 338768 338123 ...
#>  $ chr.exposure          : chr  "1" "1" "1" "1" ...
#>  $ se.exposure           : num  0.003 0.0031 0.0087 0.0029 0.003 0.0037 0.0031 0.003 0.0038 0.003 ...
#>  $ beta.exposure         : num  -0.0168 0.0201 0.0659 0.0181 0.0331 0.0497 -0.0227 0.0221 0.0209 0.0175 ...
#>  $ pos.exposure          : int  47684677 78048331 110082886 201784287 72837239 177889480 49589847 96924097 164567689 181550962 ...
#>  $ id.exposure           : chr  "ieu-a-2" "ieu-a-2" "ieu-a-2" "ieu-a-2" ...
#>  $ SNP                   : chr  "rs977747" "rs17381664" "rs7550711" "rs2820292" ...
#>  $ effect_allele.exposure: chr  "G" "C" "T" "C" ...
#>  $ other_allele.exposure : chr  "T" "T" "C" "A" ...
#>  $ eaf.exposure          : num  0.5333 0.425 0.0339 0.5083 0.6083 ...
#>  $ exposure              : chr  "Body mass index || id:ieu-a-2" "Body mass index || id:ieu-a-2" "Body mass index || id:ieu-a-2" "Body mass index || id:ieu-a-2" ...
#>  $ mr_keep.exposure      : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ pval_origin.exposure  : chr  "reported" "reported" "reported" "reported" ...
#>  $ data_source.exposure  : chr  "igd" "igd" "igd" "igd" ...
```

This returns a set of LD clumped SNPs that are GWAS significant for BMI.
You can specify various parameters for this function:

- `p1` = P-value threshold for keeping a SNP
- `clump` = Whether or not to return independent SNPs only (default is
  `TRUE`)
- `r2` = The maximum LD R-square allowed between returned SNPs
- `kb` = The distance in which to search for LD R-square values

By changing changing the `p1` parameter it is possible to obtain SNP
effects for constructing polygenic risk scores.

## Clumping

For standard two sample MR it is important to ensure that the
instruments for the exposure are independent. Once instruments have been
identified for an exposure variable, the IEU OpenGWAS database can be
used to perform clumping.

You can provide a list of SNP IDs, the SNPs will be extracted from 1000
genomes data, LD calculated between them, and amongst those SNPs that
have LD R-square above the specified threshold only the SNP with the
lowest P-value will be retained. To do this, use the following command:

``` r
bmi_exp_dat <- clump_data(bmi2014_exp_dat)
```

``` r
str(bmi_exp_dat)
#> 'data.frame':    30 obs. of  16 variables:
#>  $ SNP                   : chr  "rs10767664" "rs13078807" "rs1514175" "rs1558902" ...
#>  $ beta.exposure         : num  0.19 0.1 0.07 0.39 0.11 0.13 0.06 0.09 0.13 0.06 ...
#>  $ se.exposure           : num  0.0306 0.0204 0.0204 0.0204 0.0204 ...
#>  $ effect_allele.exposure: chr  "A" "G" "A" "A" ...
#>  $ other_allele.exposure : chr  "T" "A" "G" "T" ...
#>  $ eaf.exposure          : num  0.78 0.2 0.43 0.42 0.31 0.78 0.41 0.24 0.21 0.21 ...
#>  $ pval.exposure         : num  5e-26 4e-11 8e-14 5e-120 3e-13 ...
#>  $ units.exposure        : chr  "kg/m2" "kg/m2" "kg/m2" "kg/m2" ...
#>  $ gene.exposure         : chr  "BDNF" "CADM2" "TNNI3K" "FTO" ...
#>  $ samplesize.exposure   : int  225238 221431 207641 222476 247166 227886 209051 218439 209849 220081 ...
#>  $ exposure              : chr  "BMI" "BMI" "BMI" "BMI" ...
#>  $ mr_keep.exposure      : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ pval_origin.exposure  : chr  "reported" "reported" "reported" "reported" ...
#>  $ units.exposure_dat    : chr  "kg/m2" "kg/m2" "kg/m2" "kg/m2" ...
#>  $ id.exposure           : chr  "FXhiAH" "FXhiAH" "FXhiAH" "FXhiAH" ...
#>  $ data_source.exposure  : chr  "textfile" "textfile" "textfile" "textfile" ...
```

The
[`clump_data()`](https://mrcieu.github.io/TwoSampleMR/reference/clump_data.md)
function takes any data frame that has been formatted to be an exposure
data type of data frame. Note that for the instruments in the
MRInstruments package the SNPs are already LD clumped.

**Note:** The LD reference panel only includes SNPs (no INDELs). There
are five super-populations from which LD can be calculated, by default
European samples are used. Only SNPs with MAF \> 0.01 within-population
are available.

**NOTE:** If a variant is dropped from your unclumped data it could be
because it is absent from the reference panel. For more flexibility,
including using your own LD reference data, see here:
<https://mrcieu.github.io/ieugwasr/>
