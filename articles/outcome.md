# Outcome data

``` r
library(TwoSampleMR)
```

Once instruments for the exposure trait have been specified, those
variants need to be extracted from the outcome trait.

## Available studies in IEU GWAS database

The IEU GWAS database (IGD) contains complete GWAS summary statistics
from a large number of studies. You can browse them here:

<https://gwas.mrcieu.ac.uk/>

To obtain details about the available GWASs programmatically do the
following:

``` r
ao <- available_outcomes()
```

``` r
head(ao)
#>           id         trait ncase group_name year       author consortium
#> 1 ieu-b-5103 Schizophrenia  1234     public 2022 Trubetskoy V        PGC
#> 2 ieu-b-5102 Schizophrenia 52017     public 2022 Trubetskoy V        PGC
#> 3 ieu-b-5101 Schizophrenia 12305     public 2022 Trubetskoy V        PGC
#> 4 ieu-b-5100 Schizophrenia 64322     public 2022 Trubetskoy V        PGC
#> 5 ieu-b-5099 Schizophrenia 76755     public 2022 Trubetskoy V        PGC
#> 6 ieu-b-5098 Schizophrenia  5998     public 2022 Trubetskoy V        PGC
#>                 sex     pmid                         population  unit
#> 1 Males and Females 35396580         Hispanic or Latin American logOR
#> 2 Males and Females 35396580                           European logOR
#> 3 Males and Females 35396580                         East Asian logOR
#> 4 Males and Females 35396580                              Mixed logOR
#> 5 Males and Females 35396580                              Mixed logOR
#> 6 Males and Females 35396580 African American or Afro-Caribbean logOR
#>   sample_size       build ncontrol category subcategory      ontology
#> 1        4324 HG19/GRCh37     3090  Disease          NA MONDO:0005090
#> 2      127906 HG19/GRCh37    75889  Disease          NA MONDO:0005090
#> 3       27363 HG19/GRCh37    15058  Disease          NA MONDO:0005090
#> 4      155269 HG19/GRCh37    90947  Disease          NA MONDO:0005090
#> 5      320404 HG19/GRCh37   243649  Disease          NA MONDO:0005090
#> 6        9824 HG19/GRCh37     3826  Disease          NA MONDO:0005090
#>                                                                      note mr
#> 1                                                                    <NA> NA
#> 2                                                                    <NA> NA
#> 3                                                                    <NA> NA
#> 4                            Core - East Asian and European meta analysis NA
#> 5 Primary - meta analysis of Eur, East Asian, African American and Latino NA
#> 6                                                                    <NA> NA
#>   nsnp  doi coverage study_design priority sd
#> 1   NA <NA>     <NA>         <NA>       NA NA
#> 2   NA <NA>     <NA>         <NA>       NA NA
#> 3   NA <NA>     <NA>         <NA>       NA NA
#> 4   NA <NA>     <NA>         <NA>       NA NA
#> 5   NA <NA>     <NA>         <NA>       NA NA
#> 6   NA <NA>     <NA>         <NA>       NA NA
```

For information about authentication see
<https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication>.

The `available_outcomes` function returns a table of all the available
studies in the database. Each study has a unique ID. e.g.

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

## Extracting particular SNPs from particular studies

If we want to perform MR of BMI against coronary heart disease, we need
to identify the SNPs that influence the BMI, and then extract those SNPs
from a GWAS on coronary heart disease.

Let’s get the Locke et al 2014 instruments for BMI as an example:

``` r
bmi_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')
```

``` r
head(bmi_exp_dat)
#>   pval.exposure samplesize.exposure chr.exposure se.exposure beta.exposure
#> 1   2.18198e-08              339152            1      0.0030       -0.0168
#> 2   4.56773e-11              339065            1      0.0031        0.0201
#> 3   5.05941e-14              313621            1      0.0087        0.0659
#> 4   5.45205e-10              338768            1      0.0029        0.0181
#> 5   1.88018e-28              338123            1      0.0030        0.0331
#> 6   2.28718e-40              339078            1      0.0037        0.0497
#>   pos.exposure id.exposure        SNP effect_allele.exposure
#> 1     47684677     ieu-a-2   rs977747                      G
#> 2     78048331     ieu-a-2 rs17381664                      C
#> 3    110082886     ieu-a-2  rs7550711                      T
#> 4    201784287     ieu-a-2  rs2820292                      C
#> 5     72837239     ieu-a-2  rs7531118                      C
#> 6    177889480     ieu-a-2   rs543874                      G
#>   other_allele.exposure eaf.exposure                      exposure
#> 1                     T       0.5333 Body mass index || id:ieu-a-2
#> 2                     T       0.4250 Body mass index || id:ieu-a-2
#> 3                     C       0.0339 Body mass index || id:ieu-a-2
#> 4                     A       0.5083 Body mass index || id:ieu-a-2
#> 5                     T       0.6083 Body mass index || id:ieu-a-2
#> 6                     A       0.2667 Body mass index || id:ieu-a-2
#>   mr_keep.exposure pval_origin.exposure data_source.exposure
#> 1             TRUE             reported                  igd
#> 2             TRUE             reported                  igd
#> 3             TRUE             reported                  igd
#> 4             TRUE             reported                  igd
#> 5             TRUE             reported                  igd
#> 6             TRUE             reported                  igd
```

We now need to find a suitable GWAS for coronary heart disease. We can
search the available studies:

``` r
ao[grepl("heart disease", ao$trait), ]
#>                                 id
#> 7958                 finn-b-I9_CHD
#> 10829                   ukb-b-3983
#> 12268                ukb-e-I25_AFR
#> 14221                   ukb-b-2205
#> 14897                      ieu-a-7
#> 15028 finn-b-I9_SECONDRIGHT_EXNONE
#> 16883                   ukb-b-7436
#> 18046                    ukb-a-534
#> 18048                ukb-e-I25_CSA
#> 20207           finn-b-I9_OTHHEART
#> 22249         finn-b-I9_VHD_EXNONE
#> 22725          finn-b-I9_PULMHEART
#> 23614                      ieu-a-9
#> 24603           finn-b-FG_OTHHEART
#> 27414             ebi-a-GCST000998
#> 33264          finn-b-FG_PULMHEART
#> 35220                 ukb-d-I9_CHD
#> 35231        finn-b-I9_OTHILLHEART
#> 36931                   ukb-b-8184
#> 37302 finn-b-I9_OTHILLHEART_EXNONE
#> 37894                   ukb-b-1668
#> 38319                finn-b-I9_IHD
#> 38448            finn-b-I9_RHEUFEV
#> 38602                      ieu-a-8
#> 41994          finn-b-I9_ISCHHEART
#> 42095                 ukb-d-I9_IHD
#> 42882        finn-b-I9_SECONDRIGHT
#> 43747           ukb-d-I9_CHD_NOREV
#> 45294                      ieu-a-6
#> 46126                finn-b-I9_VHD
#> 46689                  ukb-b-16606
#>                                                                                                                           trait
#> 7958                                                                                         Major coronary heart disease event
#> 10829                                                Diagnoses - main ICD10: I25.9 Chronic ischaemic heart disease, unspecified
#> 12268                                                                                       I25 Chronic ischaemic heart disease
#> 14221 Diagnoses - secondary ICD10: Z82.4 Family history of ischaemic heart disease and other diseases of the circulatory system
#> 14897                                                                                                    Coronary heart disease
#> 15028                                                                      Secondary right heart disease (no controls excluded)
#> 16883                                                          Diagnoses - secondary ICD10: I25.1 Atherosclerotic heart disease
#> 18046                                                               Diagnoses - main ICD10: I25 Chronic ischaemic heart disease
#> 18048                                                                                       I25 Chronic ischaemic heart disease
#> 20207                                                                                        Other heart diseases (I9_OTHHEART)
#> 22249                                                   Valvular heart disease including rheumatic fever (no controls excluded)
#> 22725                                                                Pulmonary heart disease, diseases of pulmonary circulation
#> 23614                                                                                                    Coronary heart disease
#> 24603                                                                                        Other heart diseases (FG_OTHHEART)
#> 27414                                                                                                    Coronary heart disease
#> 33264                                                                                                   Pulmonary heart disease
#> 35220                                                                                        Major coronary heart disease event
#> 35231                                                                                       Other or ill-defined heart diseases
#> 36931                                           Diagnoses - secondary ICD10: I25.9 Chronic ischaemic heart disease, unspecified
#> 37302                                                                Other or ill-defined heart diseases (no controls excluded)
#> 37894                                                               Diagnoses - main ICD10: I25.1 Atherosclerotic heart disease
#> 38319                                                                                  Ischaemic heart disease, wide definition
#> 38448                                                                                        Rheumatic fever incl heart disease
#> 38602                                                                                                    Coronary heart disease
#> 41994                                                                                                   Ischemic heart diseases
#> 42095                                                                                  Ischaemic heart disease, wide definition
#> 42882                                                                                             Secondary right heart disease
#> 43747                                                           Major coronary heart disease event excluding revascularizations
#> 45294                                                                                                    Coronary heart disease
#> 46126                                                                          Valvular heart disease including rheumatic fever
#> 46689                                         Diagnoses - secondary ICD10: I25.8 Other forms of chronic ischaemic heart disease
#>       ncase group_name year       author        consortium               sex
#> 7958  21012     public 2021           NA                NA Males and Females
#> 10829  1195     public 2018 Ben Elsworth           MRC-IEU Males and Females
#> 12268   302     public 2020 Pan-UKB team                NA Males and Females
#> 14221  9330     public 2018 Ben Elsworth           MRC-IEU Males and Females
#> 14897 60801     public 2015       Nikpay CARDIoGRAMplusC4D Males and Females
#> 15028   428     public 2021           NA                NA Males and Females
#> 16883  5771     public 2018 Ben Elsworth           MRC-IEU Males and Females
#> 18046  8755     public 2017        Neale         Neale Lab Males and Females
#> 18048  1205     public 2020 Pan-UKB team                NA Males and Females
#> 20207 62081     public 2021           NA                NA Males and Females
#> 22249 38209     public 2021           NA                NA Males and Females
#> 22725  4564     public 2021           NA                NA Males and Females
#> 23614 63746     public 2013     Deloukas CARDIoGRAMplusC4D Males and Females
#> 24603 58173     public 2021           NA                NA Males and Females
#> 27414 22233     public 2011  Schunkert H                NA                NA
#> 33264  4185     public 2021           NA                NA Males and Females
#> 35220 10157     public 2018    Neale lab                NA Males and Females
#> 35231   713     public 2021           NA                NA Males and Females
#> 36931  5861     public 2018 Ben Elsworth           MRC-IEU Males and Females
#> 37302   713     public 2021           NA                NA Males and Females
#> 37894 12171     public 2018 Ben Elsworth           MRC-IEU Males and Females
#> 38319 31640     public 2021           NA                NA Males and Females
#> 38448   573     public 2021           NA                NA Males and Females
#> 38602 22233     public 2011  Schunkert H        CARDIoGRAM Males and Females
#> 41994 30952     public 2021           NA                NA Males and Females
#> 42095 20857     public 2018    Neale lab                NA Males and Females
#> 42882   428     public 2021           NA                NA Males and Females
#> 43747 10157     public 2018    Neale lab                NA Males and Females
#> 45294 15420     public 2011        Peden               C4D Males and Females
#> 46126 38209     public 2021           NA                NA Males and Females
#> 46689  5738     public 2018 Ben Elsworth           MRC-IEU Males and Females
#>           pmid                         population     unit sample_size
#> 7958        NA                           European       NA          NA
#> 10829       NA                           European       SD      463010
#> 12268       NA African American or Afro-Caribbean       NA        6636
#> 14221       NA                           European       SD      463010
#> 14897 26343387                              Mixed log odds      184305
#> 15028       NA                           European       NA          NA
#> 16883       NA                           European       SD      463010
#> 18046       NA                           European       SD      337199
#> 18048       NA                        South Asian       NA        8876
#> 20207       NA                           European       NA          NA
#> 22249       NA                           European       NA          NA
#> 22725       NA                           European       NA          NA
#> 23614 23202125                              Mixed log odds      194427
#> 24603       NA                           European       NA          NA
#> 27414 21378990                           European    logOR       86995
#> 33264       NA                           European       NA          NA
#> 35220       NA                           European       NA      361194
#> 35231       NA                           European       NA          NA
#> 36931       NA                           European       SD      463010
#> 37302       NA                           European       NA          NA
#> 37894       NA                           European       SD      463010
#> 38319       NA                           European       NA          NA
#> 38448       NA                           European       NA          NA
#> 38602 21378990                           European log odds       86995
#> 41994       NA                           European       NA          NA
#> 42095       NA                           European       NA      361194
#> 42882       NA                           European       NA          NA
#> 43747       NA                           European       NA      361194
#> 45294 21378988                              Mixed log odds       30482
#> 46126       NA                           European       NA          NA
#> 46689       NA                           European       SD      463010
#>             build ncontrol category    subcategory ontology
#> 7958  HG19/GRCh37   197780   Binary             NA       NA
#> 10829 HG19/GRCh37   461815   Binary             NA       NA
#> 12268 HG19/GRCh37     6334   Binary             NA       NA
#> 14221 HG19/GRCh37   453680   Binary             NA       NA
#> 14897 HG19/GRCh37   123504  Disease Cardiovascular       NA
#> 15028 HG19/GRCh37   218364   Binary             NA       NA
#> 16883 HG19/GRCh37   457239   Binary             NA       NA
#> 18046 HG19/GRCh37   328444       NA             NA       NA
#> 18048 HG19/GRCh37     7671   Binary             NA       NA
#> 20207 HG19/GRCh37   156711   Binary             NA       NA
#> 22249 HG19/GRCh37   180583   Binary             NA       NA
#> 22725 HG19/GRCh37   214228   Binary             NA       NA
#> 23614 HG19/GRCh37   130681  Disease Cardiovascular       NA
#> 24603 HG19/GRCh37   160619   Binary             NA       NA
#> 27414 HG19/GRCh37    64762       NA             NA       NA
#> 33264 HG19/GRCh37   214607   Binary             NA       NA
#> 35220 HG19/GRCh37   351037   Binary             NA       NA
#> 35231 HG19/GRCh37   156711   Binary             NA       NA
#> 36931 HG19/GRCh37   457149   Binary             NA       NA
#> 37302 HG19/GRCh37   218079   Binary             NA       NA
#> 37894 HG19/GRCh37   450839   Binary             NA       NA
#> 38319 HG19/GRCh37   187152   Binary             NA       NA
#> 38448 HG19/GRCh37   218219   Binary             NA       NA
#> 38602 HG19/GRCh37    64762  Disease Cardiovascular       NA
#> 41994 HG19/GRCh37   187840   Binary             NA       NA
#> 42095 HG19/GRCh37   340337   Binary             NA       NA
#> 42882 HG19/GRCh37   214228   Binary             NA       NA
#> 43747 HG19/GRCh37   351037   Binary             NA       NA
#> 45294 HG19/GRCh37    15062  Disease Cardiovascular       NA
#> 46126 HG19/GRCh37   156711   Binary             NA       NA
#> 46689 HG19/GRCh37   457272   Binary             NA       NA
#>                                                                                       note
#> 7958                                                                                I9_CHD
#> 10829 41202#I259: Output from GWAS pipeline using Phesant derived variables from UKBiobank
#> 12268                                                                                   NA
#> 14221 41204#Z824: Output from GWAS pipeline using Phesant derived variables from UKBiobank
#> 14897                                                                                 <NA>
#> 15028                                                                I9_SECONDRIGHT_EXNONE
#> 16883 41204#I251: Output from GWAS pipeline using Phesant derived variables from UKBiobank
#> 18046                                                                                   NA
#> 18048                                                                                   NA
#> 20207                                                                          I9_OTHHEART
#> 22249                                                                        I9_VHD_EXNONE
#> 22725                                                                         I9_PULMHEART
#> 23614                                                                                 <NA>
#> 24603                                                                          FG_OTHHEART
#> 27414                                                                                   NA
#> 33264                                                                         FG_PULMHEART
#> 35220                                                                                   NA
#> 35231                                                                       I9_OTHILLHEART
#> 36931 41204#I259: Output from GWAS pipeline using Phesant derived variables from UKBiobank
#> 37302                                                                I9_OTHILLHEART_EXNONE
#> 37894 41202#I251: Output from GWAS pipeline using Phesant derived variables from UKBiobank
#> 38319                                                                               I9_IHD
#> 38448                                                                           I9_RHEUFEV
#> 38602                                                                                 <NA>
#> 41994                                                                         I9_ISCHHEART
#> 42095                                                                                   NA
#> 42882                                                                       I9_SECONDRIGHT
#> 43747                                                                                   NA
#> 45294                                                                                 <NA>
#> 46126                                                                               I9_VHD
#> 46689 41204#I258: Output from GWAS pipeline using Phesant derived variables from UKBiobank
#>       mr     nsnp  doi coverage study_design priority sd
#> 7958   1 16380466 <NA>     <NA>         <NA>        0 NA
#> 10829  1  9851867 <NA>     <NA>         <NA>        1 NA
#> 12268  1 15478580 <NA>     <NA>         <NA>        0 NA
#> 14221  1  9851867 <NA>     <NA>         <NA>        1 NA
#> 14897  1  9455779 <NA>     <NA>         <NA>        1 NA
#> 15028  1 16380466 <NA>     <NA>         <NA>        0 NA
#> 16883  1  9851867 <NA>     <NA>         <NA>        1 NA
#> 18046  1 10894596 <NA>     <NA>         <NA>        1 NA
#> 18048  1  9811287 <NA>     <NA>         <NA>        0 NA
#> 20207  1 16380466 <NA>     <NA>         <NA>        0 NA
#> 22249  1 16380466 <NA>     <NA>         <NA>        0 NA
#> 22725  1 16380466 <NA>     <NA>         <NA>        0 NA
#> 23614  1    79129 <NA>     <NA>         <NA>        1 NA
#> 24603  1 16380466 <NA>     <NA>         <NA>        0 NA
#> 27414  1  2415020 <NA>     <NA>         <NA>        0 NA
#> 33264  1 16380466 <NA>     <NA>         <NA>        0 NA
#> 35220  1 13295130 <NA>     <NA>         <NA>        0 NA
#> 35231  1 16380177 <NA>     <NA>         <NA>        0 NA
#> 36931  1  9851867 <NA>     <NA>         <NA>        1 NA
#> 37302  1 16380466 <NA>     <NA>         <NA>        0 NA
#> 37894  1  9851867 <NA>     <NA>         <NA>        1 NA
#> 38319  1 16380466 <NA>     <NA>         <NA>        0 NA
#> 38448  1 16380466 <NA>     <NA>         <NA>        0 NA
#> 38602  1  2420361 <NA>     <NA>         <NA>        2 NA
#> 41994  1 16380466 <NA>     <NA>         <NA>        0 NA
#> 42095  1 13586589 <NA>     <NA>         <NA>        0 NA
#> 42882  1 16380459 <NA>     <NA>         <NA>        0 NA
#> 43747  1 13295130 <NA>     <NA>         <NA>        0 NA
#> 45294  1   540233 <NA>     <NA>         <NA>        3 NA
#> 46126  1 16380358 <NA>     <NA>         <NA>        0 NA
#> 46689  1  9851867 <NA>     <NA>         <NA>        1 NA
```

The most recent CARDIOGRAM GWAS is ID number `ieu-a-7`. We can extract
the BMI SNPs from this GWAS as follows:

``` r
chd_out_dat1 <- extract_outcome_data(
  snps = bmi_exp_dat$SNP,
  outcomes = 'ieu-a-7'
)
```

The
[`extract_outcome_data()`](https://mrcieu.github.io/TwoSampleMR/reference/extract_outcome_data.md)
function is flexible. The `snps` argument only requires an array of
rsIDs, and the `outcomes` argument can be a vector of outcomes, e.g.

``` r
chd_out_dat2 <- extract_outcome_data(
  snps = c("rs234", "rs17097147"),
  outcomes = c('ieu-a-2', 'ieu-a-7')
)
```

will extract the two SNPs from each of the outcomes `ieu-a-2` and
`ieu-a-7`.

## LD proxies

By default if a particular requested SNP is not present in the outcome
GWAS then a SNP (proxy) that is in LD with the requested SNP (target)
will be searched for instead. LD proxies are defined using 1000 genomes
European sample data. The effect of the proxy SNP on the outcome is
returned, along with the proxy SNP, the effect allele of the proxy SNP,
and the corresponding allele (in phase) for the target SNP.

The parameters for handling LD proxies are as follows:

- `proxies` = TRUE or FALSE (TRUE by default)
- `rsq` = numeric value of minimum rsq to find a proxy. Default is 0.8,
  minimum is 0.6
- `palindromes` = Allow palindromic SNPs? Default is 1 (yes)
- `maf_threshold` = If palindromes allowed then what is the maximum
  minor allele frequency of palindromes allowed? Default is 0.3.

## Using local GWAS summary data

If you have GWAS summary data that is not present in IEU GWAS database,
this can still be used to perform analysis.

Supposing there is a GWAS summary file called “gwas_summary.csv” with
e.g. 2 million rows and it looks like this:

    rsid,effect,SE,a1,a2,a1_freq,p-value,Units,Gene,n
    rs10767664,0.19,0.030612245,A,T,0.78,5.00E-26,kg/m2,BDNF,225238
    rs13078807,0.1,0.020408163,G,A,0.2,4.00E-11,kg/m2,CADM2,221431
    rs1514175,0.07,0.020408163,A,G,0.43,8.00E-14,kg/m2,TNNI3K,207641
    rs1558902,0.39,0.020408163,A,T,0.42,5.00E-120,kg/m2,FTO,222476
    ...
    ...

To extract the exposure SNPs from this data, we would use the following
command:

``` r
outcome_dat <- read_outcome_data(
  snps = bmi_exp_dat$SNP,
  filename = "gwas_summary.csv",
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
```

This returns an outcome data frame with only the SNPs that were
requested (if those SNPs were present in the “gwas_summary.csv” file).

## Outcome data format

The `extract_outcome_data` function returns a table of SNP effects for
the requested SNPs on the requested outcomes. The format of the data is
similar to the exposure data format, except the main columns are as
follows:

- `SNP`
- `beta.outcome`
- `se.outcome`
- `samplesize.outcome`
- `ncase.outcome`
- `ncontrol.outcome`
- `pval.outcome`
- `eaf.outcome`
- `effect_allele.outcom`
- `other_allele.outcome`
- `units.outcome`
- `outcome`
- `consortium.outcome`
- `year.outcome`
- `pmid.outcome`
- `id.outcome`
- `originalname.outcome`
- `proxy.outcome`
- `target_snp.outcome`
- `proxy_snp.outcome`
- `target_a1.outcome`
- `target_a2.outcome`
- `proxy_a1.outcome`
- `proxy_a2.outcome`
- `mr_keep.outcome`
- `data_source.outcome`

## More advanced use of local data

We have developed a summary data format called “GWAS VCF”, which is
designed to store GWAS results in a strict and performant way. It is
possible to use this format with the TwoSampleMR package. Going down
this avenue also allows you to use LD proxy functionality using your own
LD reference files (or ones that we provide). For more details, see this
package that explains the format and how to query it in R:

<https://github.com/mrcieu/gwasvcf>

and this package for how to connect the data to other packages including
TwoSampleMR

<https://github.com/MRCIEU/gwasglue>
