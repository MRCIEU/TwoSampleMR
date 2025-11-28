# Harmonise data

``` r
library(TwoSampleMR)
```

## Introduction

The exposure data and outcome data are now obtained, e.g.:

``` r
bmi_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')
chd_out_dat <- extract_outcome_data(snps = bmi_exp_dat$SNP, outcomes = 'ieu-a-7')
```

but it is important to harmonise the effects. This means that the effect
of a SNP on the exposure and the effect of that SNP on the outcome must
each correspond to the same allele.

**Note:** The IEU GWAS database contains data that is already
harmonised, meaning that the non-effect allele is aligned to the human
genome reference sequence (build 37). It’s still recommended to
harmonise, but in principle everything should be on the forward strand
and effect alleles always relating to the same allele. Some
discrepancies could arise if there are multi-allelic variants that are
represented as different bi-allelic variants in different studies.

To harmonise the exposure and outcome data, do the following:

``` r
dat <- harmonise_data(
  exposure_dat = bmi_exp_dat,
  outcome_dat = chd_out_dat
)
#> Harmonising Body mass index || id:ieu-a-2 (ieu-a-2) and Coronary heart disease || id:ieu-a-7 (ieu-a-7)
```

This creates a new data frame that has the exposure data and outcome
data combined.

If there were 3 exposure traits and 3 outcome traits then there will be
9 sets of harmonisations being performed - harmonising the SNP effects
of exposure trait 1 against outcome trait 1; exposure trait 1 against
outcome trait 2; and so on.

## Dealing with strand issues

Recent GWASs typically present the effects of a SNP in reference to the
allele on the forward strand. But as reference panels are updated the
forward strand sometimes changes, and GWASs from a few years ago aren’t
guaranteed to be using forward strand conventions.

Some examples are shown below:

### Correct, unambiguous

    exposure effect = 0.5
    effect allele = A
    other allele = G

    outcome effect = 0.05
    effect allele = A
    other allele = G

Here the effect allele on the exposure and the outcome is the same

### Incorrect reference, unambiguous

    exposure effect = 0.5
    effect allele = A
    other allele = G

    outcome effect = -0.05
    effect allele = C
    other allele = T

Here the outcome GWAS is presenting the effect for the alternate allele
on the reverse strand. We need to flip the outcome effect to 0.05 to
correspond to the same allele as the exposure GWAS on the forward
strand.

### Ambiguous

    exposure effect = 0.5
    effect allele = A
    other allele = G

    outcome effect = -0.05
    effect allele = A
    other allele = C

Here the alleles do not correspond for the same SNP, so this SNP will be
discarded from the analysis.

### Palindromic SNP, inferrable

    exposure effect = 0.5
    effect allele = A
    other allele = T
    effect allele frequency = 0.11

    outcome effect = -0.05
    effect allele = A
    other allele = T
    effect allele frequency = 0.91

Here the alleles correspond, but it is a palindromic SNP, such that the
alleles on the forward strand are the same as on the reverse strand (A/T
on forward is T/A on the reverse). However, the allele frequency of the
effect allele gives us information - if the outcome effect allele (A)
were on the forward strand we would expect it to have a low allele
frequency, but given it has a high frequency (0.91) we infer that the
outcome GWAS is presenting the effect on the reverse strand for the
alternative allele. We would flip the effect to 0.05 for the outcome
GWAS.

### Palindromic SNP, not inferrable

    exposure effect = 0.5
    effect allele = A
    other allele = T
    effect allele frequency = 0.50

    outcome effect = -0.05
    effect allele = A
    other allele = T
    effect allele frequency = 0.50

This is similar to the above, except the allele frequency no longer
gives us information about the strand. We would discard this SNP. This
is done for any palindromic SNPs that have minor allele frequency above
0.42.

### Options

There are three options to harmonising the data.

1.  Assume all alleles are presented on the forward strand
2.  Try to infer the forward strand alleles using allele frequency
    information
3.  Correct the strand for non-palindromic SNPs, but drop all
    palindromic SNPs

By default, the `harmonise_data` function uses option 2, but this can be
modified using the `action` argument,
e.g. `harmonise_data(exposure_dat, outcome_dat, action = 3)`.

## Drop duplicate exposure-outcome summary sets

After data harmonisation, users may find that their dataset contains
duplicate exposure-outcome summary sets. This can arise, for example,
when a GWAS consortium has released multiple results from separate GWAS
analyses for the same trait. For example, there are multiple GWAS
summary datasets for body mass index and coronary heart disease:

``` r
ao <- available_outcomes()
```

``` r
ao[ao$trait == "Body mass index", c("trait", "id", "pmid", "author", "sample_size", "nsnp")]
#>                 trait                 id     pmid                    author
#> 3958  Body mass index ebi-a-GCST90103751 35051171                   Wong HS
#> 4015  Body mass index ebi-a-GCST90095039 35399580 Fern<U+00E1>ndez-Rhodes L
#> 4020  Body mass index ebi-a-GCST90095034 35399580 Fern<U+00E1>ndez-Rhodes L
#> 6032  Body mass index ebi-a-GCST90029007 29892013                    Loh PR
#> 6821  Body mass index ebi-a-GCST90025994 34226706                 Barton AR
#> 7045  Body mass index ebi-a-GCST90018947 34594039                  Sakaue S
#> 7259  Body mass index ebi-a-GCST90018727 34594039                  Sakaue S
#> 10738 Body mass index           ieu-a-94 23754948                Randall JC
#> 12717 Body mass index            ieu-a-2 25673413                  Locke AE
#> 14254 Body mass index           ieu-a-95 23754948                Randall JC
#> 16676 Body mass index          ieu-a-974 25673413                  Locke AE
#> 19343 Body mass index            bbj-a-3 28892062                Ishigaki K
#> 26475 Body mass index   ebi-a-GCST006368 30108127               Hoffmann TJ
#> 28065 Body mass index         ieu-b-4815       NA                   Howe LJ
#> 28419 Body mass index            bbj-a-2 28892062                Ishigaki K
#> 32371 Body mass index         ieu-b-4816       NA                   Howe LJ
#> 33869 Body mass index          ieu-a-785 25673413                  Locke AE
#> 39217 Body mass index   ebi-a-GCST002783 25673413                  Locke AE
#> 40065 Body mass index            bbj-a-1 28892062                Ishigaki K
#> 43262 Body mass index   ebi-a-GCST004904 28892062                 Akiyama M
#> 43743 Body mass index   ebi-a-GCST006802 26961502                   Wood AR
#> 47894 Body mass index          ieu-a-835 25673413                  Locke AE
#> 48917 Body mass index   ebi-a-GCST008025 31217584                 Wojcik GL
#> 49208 Body mass index         ieu-a-1089 26961502                      Wood
#>       sample_size     nsnp
#> 3958        21930  6370138
#> 4015       330793  2401077
#> 4020        56161  8764141
#> 6032       532396 11973091
#> 6821       457756  4238669
#> 7045       359983 19066885
#> 7259       163835 12502877
#> 10738       60586  2736876
#> 12717      339224  2555511
#> 14254       73137  2736876
#> 16676      171977  2494613
#> 19343       72390  6108953
#> 26475      315347 27854527
#> 28065       51852       NA
#> 28419       85894  6108953
#> 32371       99998  7191606
#> 33869      152893  2477659
#> 39217      236781  2529499
#> 40065      158284  5961600
#> 43262      158284  5952516
#> 43743      119688  8580466
#> 47894      322154  2554668
#> 48917       21955 34343880
#> 49208      120286  8654252
ao[ao$trait == "Coronary heart disease", c("trait", "id", "pmid", "author", "ncase", "ncontrol", "nsnp")]
#>                        trait               id     pmid      author ncase
#> 14897 Coronary heart disease          ieu-a-7 26343387      Nikpay 60801
#> 23614 Coronary heart disease          ieu-a-9 23202125    Deloukas 63746
#> 27414 Coronary heart disease ebi-a-GCST000998 21378990 Schunkert H 22233
#> 38602 Coronary heart disease          ieu-a-8 21378990 Schunkert H 22233
#> 45294 Coronary heart disease          ieu-a-6 21378988       Peden 15420
#>       ncontrol    nsnp
#> 14897   123504 9455779
#> 23614   130681   79129
#> 27414    64762 2415020
#> 38602    64762 2420361
#> 45294    15062  540233
```

There are therefore multiple potential combinations of body mass index
and coronary heart disease, which would likely lead to duplicate MR
analyses. We recommend that users prune their datasets so that only the
exposure-outcome combination with the highest expected power is
retained. This can be done by selecting the exposure-outcome summary set
with the largest sample size for the outcome, using the power_prune
function:

``` r
dat <- power_prune(dat, method = 1, dist.outcome = "binary")
```

This drops the duplicate exposure-outcome sets with the smaller outcome
sample size (number of cases for binary outcomes). Remaining duplicates
are then dropped on the basis of the exposure sample size. However, if
there are a large number of SNPs available to instrument an exposure,
the outcome GWAS with the better SNP coverage may provide better power
than the outcome GWAS with the larger sample size. This can occur, for
example, if the larger outcome GWAS has used a targeted genotyping
array. In such instances, it may be better to prune studies on the basis
of instrument strength (i.e. variation in exposure explained by the
instrumental SNPs) as well as sample size. This can be done by setting
the method argument to 2:

``` r
dat <- power_prune(dat, method = 2, dist.outcome = "binary")
```

This procedure drops duplicate exposure-outcome sets on the basis of
instrument strength and sample size, and assumes that the SNP-exposure
effects correspond to a continuous trait with a normal distribution
(i.e. exposure should not be binary). The SNP-outcome effects can
correspond to either a binary or continuous trait (default behaviour is
to assume a binary distribution). If the exposure is binary then method
1 should be used.
