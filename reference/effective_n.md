# Estimate the effective sample size in a case control study

Computes the effective sample size as the harmonic mean of the number of
cases and controls, \$\$N\_{eff} = \frac{2}{1/N\_{cases} +
1/N\_{controls}},\$\$ as given by Winkler et al. (2014)
<https://www.nature.com/articles/nprot.2014.071>.

## Usage

``` r
effective_n(ncase, ncontrol)
```

## Arguments

- ncase:

  Vector of number of cases

- ncontrol:

  Vector of number of controls

## Value

Vector of effective sample size

## Details

Note that other software, such as METAL (Willer et al., 2010,
[doi:10.1093/bioinformatics/btq340](https://doi.org/10.1093/bioinformatics/btq340)
), uses \\4/(1/N\_{cases} + 1/N\_{controls})\\, which is twice this
value. With that definition a balanced study (equal numbers of cases and
controls) has an effective sample size equal to its total sample size,
whereas with this function it is equal to half the total sample size,
i.e., the number of cases. If you require the METAL definition multiply
the result of this function by 2.

This function is used by
[`add_rsq()`](https://mrcieu.github.io/TwoSampleMR/reference/add_rsq.md)
to obtain the sample sizes for binary traits (those with units of
`"log odds"`) which are then used by
[`directionality_test()`](https://mrcieu.github.io/TwoSampleMR/reference/directionality_test.md)
and
[`steiger_filtering()`](https://mrcieu.github.io/TwoSampleMR/reference/steiger_filtering.md).
