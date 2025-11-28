# Changelog

## TwoSampleMR v0.6.24

(Release date 2025-10-31)

- We once again obtain **mr.raps** from its default branch. We issue a
  warning if the user has a version less than version 0.4.3 installed.

## TwoSampleMR v0.6.23

(Release date 2025-10-31)

- We now obtain **mr.raps** (version 0.4.3) from the open pull request
  which includes a fix. We will switch back to the default branch once
  that pull request is merged.

## TwoSampleMR v0.6.22

(Release date 2025-08-27)

- Additional reformatting of code base.
- Minor updates to GitHub Actions workflows.

## TwoSampleMR v0.6.21

(Release date 2025-08-11)

- Fixed typos in Steiger related documentation (thanks to
  [@eleanorsanderson](https://github.com/eleanorsanderson) and Kaitlin
  Wade for providing the motivating example for this).
- Some code cleanup to make debugging easier.
- Reduced size of installed package by removing some unused files.

## TwoSampleMR v0.6.20

(Release date 2025-08-01)

- TwoSampleMR now requires **ieugwasr** version 1.1.0 which provides
  improved error messages.
- Changed error handling - now if the API gives an error code it is
  propagated as an error with message in **TwoSampleMR**.
- Headers sent to API are updated.

## TwoSampleMR v0.6.19

(Release date 2025-07-21)

- Further fixes for the forthcoming version 4 of ggplot2.

## TwoSampleMR v0.6.18

(Release date 2025-07-15)

- Fixes to
  [`mr_forest_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_forest_plot.md),
  [`forest_plot_1_to_many()`](https://mrcieu.github.io/TwoSampleMR/reference/forest_plot_1_to_many.md),
  [`mr_leaveoneout_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_leaveoneout_plot.md),
  and
  [`mr_scatter_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_scatter_plot.md)
  for compatibility with the forthcoming release of version 4 of
  **ggplot2**.
- Made some minor improvements in the *Perform MR* vignette including
  adding back several of the forest plots.

## TwoSampleMR v0.6.17

(Release date 2025-06-22)

- Continued amends to the plot for MR-GRIP (thanks
  [@fdudbridge](https://github.com/fdudbridge))

## TwoSampleMR v0.6.16

(Release date 2025-06-05)

- Added weak instrument adjustment estimate and intercept to returned
  output of
  [`mr_grip()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_grip.md)
  (thanks [@fdudbridge](https://github.com/fdudbridge))
- Amended
  [`mr_scatter_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_scatter_plot.md)
  to issue a message about the scale the MR-GRIP estimate is plotted on
  (thanks [@fdudbridge](https://github.com/fdudbridge))

## TwoSampleMR v0.6.15

(Release date 2025-05-01)

- Bumped the minimum required version of R to 4.1.0. This is due to a
  dependency (the scales package) of a dependency (the ggplot2 package),
  which now has this requirement.

## TwoSampleMR v0.6.14

(Release date 2025-03-28)

- Minor amends to the
  [`mr_grip()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_grip.md)
  returned object names (thanks
  [@fdudbridge](https://github.com/fdudbridge))
- Fixed some typos in the helpfiles and vignettes

## TwoSampleMR v0.6.13

(Release date 2025-03-26)

- Added
  [`mr_grip()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_grip.md)
  function which implements the MR-GRIP (modified MR-Egger with the
  Genotype Recoding Invariance Property) method of Dudbridge and Bowden
  et al. (2025). The new method can be accessed by
  `mr(dat, method_list = "mr_grip")` or it can be added to the default
  list of methods with
  `mr(dat, method_list = c(subset(mr_method_list(), use_by_default)$obj, "mr_grip"))`.
- Added Pub Med IDs for more of the methods.
- The
  [`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md)
  function no longer causes a stack overflow when its `dat` argument is
  not a variable (thanks to
  [@DarwinAwardWinner](https://github.com/DarwinAwardWinner))

## TwoSampleMR v0.6.12

(Release date 2025-03-18)

- Fixed a bug in
  [`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md)
  (thanks to [@lea-urpa](https://github.com/lea-urpa))
- For MVMR data extraction in
  [`mv_extract_exposures()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_extract_exposures.md)
  and
  [`mv_extract_exposures_local()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_extract_exposures_local.md),
  SNPs that do not pass harmonisation between exposures in MVMR are now
  dropped (thanks [@yatest](https://github.com/yatest))
- The comments about **mr.raps** in the Perform MR vignette have been
  updated

## TwoSampleMR v0.6.11

(Release date 2025-03-06)

- The tests using **mr.raps** are now run conditionally.
- The **rsnps** package has been removed from Suggests because the
  configuration of **mr.raps** has been corrected.

## TwoSampleMR v0.6.10

(Release date 2025-03-03)

- The **mr.raps** package was archived from CRAN on 2025-03-01. A later
  version (0.4.1) than was on CRAN (0.2) is available from its GitHub
  repo <https://github.com/qingyuanzhao/mr.raps>. We believe this later
  version still works as expected in **TwoSampleMR**. However, version
  0.4.1 depends upon **rsnps** package which is itself now only
  available on its GitHub repo. Hence, we have added both packages to
  our <https://mrcieu.r-universe.dev/> and added that to the
  `Additional_repositories` field in the DESCRIPTION. We have also moved
  mr.raps to the soft dependency list (Suggests list) and added rsnps to
  our Suggests list, both with remotes.
- Fix to
  [`default_parameters()`](https://mrcieu.github.io/TwoSampleMR/reference/default_parameters.md)
  helpfile

## TwoSampleMR v0.6.9

(Release date 2025-02-05)

- Fixed a bug in
  [`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md)
  when the `log_pval` argument was set to `TRUE`. The specified p-value
  column is now used (thanks to
  [@luddeluddis](https://github.com/luddeluddis))
- Amend references to MR-Base to OpenGWAS

## TwoSampleMR v0.6.8

(Release date 2024-09-06)

- Replaced some [`unique()`](https://rdrr.io/r/base/unique.html) calls
  in
  [`power_prune()`](https://mrcieu.github.io/TwoSampleMR/reference/power_prune.md)
  with [`mean()`](https://rdrr.io/r/base/mean.html) to ensure scalar
  `iv.se` values (thanks [@phageghost](https://github.com/phageghost))
- Slightly improved formatting of code base

## TwoSampleMR v0.6.7

(Release date 2024-08-21)

- Update OpenGWAS API URLs
- Minor tweak to `R CMD check` GitHub Actions due to the rjson hard
  dependency of the MendelianRandomization package now requiring R \>=
  4.4.0
- Add dark mode to pkgdown site

## TwoSampleMR v0.6.6

(Release date 2024-07-06)

- Improve a test
- Improve permissions in GitHub Actions workflows
- Bump minimum required version of **ieugwasr** to 1.0.1
- Made some amends to the code to bring it more in line with **lintr**
  recommendations
- Added omitted **tidyr** soft dependency

## TwoSampleMR v0.6.5

(Release date: 2024-06-30)

- Bumped version of **roxygen2** for creating package documentation
- Update the earliest version of R in the `R CMD check` GitHub Actions
  workflow to be 4.3.2. This is because the **meta** dependency depends
  on **lme4**, and the recent 1.1-35.4 release of **lme4** requires
  **Matrix** 1.6-2 which was released a few days after R 4.3.2.
- Made package tests more robust to non-response from the OpenGWAS API

## TwoSampleMR v0.6.4

(Release date: 2024-06-05)

- Update installation instructions in README.md
- Fixed a bug in which the wrong indels recoding function was called
  (thanks [@ruochiz](https://github.com/ruochiz))

## TwoSampleMR v0.6.3

(Release date: 2024-05-23)

- Update package startup message

## TwoSampleMR v0.6.2

(Release date: 2024-05-09)

- [`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md)
  now errors if it detects its `dat` object is of class `'data.table'`
  and issues a message informing the user to make their dat object
  simply a `'data.frame'` (thanks to Si Fang
  [@sifang1678](https://github.com/sifang1678))

## TwoSampleMR v0.6.1

(Release date: 2024-04-30)

- The **MendelianRandomization** package has been moved from a hard
  dependency to a soft dependency. This is because its dependency
  package **Matrix** now requires R \>= 4.4.0. Making
  **MendelianRandomization** a soft dependency means we don’t need to
  make TwoSampleMR have the same requirement.

## TwoSampleMR v0.6.0

(Release date: 2024-04-22)

- TwoSampleMR now uses the CRAN version of the **ieugwasr** package.
  Importantly this includes the new authentication system for the
  OpenGWAS API. Please see
  <https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication>
  for more information about how to set this up.

## TwoSampleMR v0.5.11

(Release date: 2024-03-21)

- In
  [`mr_leaveoneout_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_leaveoneout_plot.md)
  and
  [`mr_forest_plot()`](https://mrcieu.github.io/TwoSampleMR/reference/mr_forest_plot.md)
  amended `size` argument to `linewidth` as per ggplot2 version 3.4.0.
- Add some datasets such that tests, continuous integration services,
  and creation of the vignettes don’t rely on the availability of the
  OpenGWAS server.
- Various improvements to helpfiles.

## TwoSampleMR v0.5.10

(Release date: 2024-02-20)

- Added `bfile` and `plink_bin` arguments to
  [`clump_data()`](https://mrcieu.github.io/TwoSampleMR/reference/clump_data.md)
- Improvements to file reading and dataset formatting capabilities of
  [`mv_extract_exposures_local()`](https://mrcieu.github.io/TwoSampleMR/reference/mv_extract_exposures_local.md)
  to create the multiple exposure dataset

## TwoSampleMR v0.5.9

(Release date: 2024-02-01)

- Fixed a minor issue in
  [`dat_to_RadialMR()`](https://mrcieu.github.io/TwoSampleMR/reference/dat_to_RadialMR.md)
- Minor improvements to
  [`make_dat()`](https://mrcieu.github.io/TwoSampleMR/reference/make_dat.md)
  default arguments and helpfile
- Minor improvements to package tests
- Amendments to GitHub Actions workflows
- Updated several URLs which had changed

## TwoSampleMR v0.5.8

(Release date: 2023-11-16)

- Improved speed of harmonisation using data.table functions (thanks
  [@nicksunderland](https://github.com/nicksunderland))
- Updated URL to R-CMD-check README badge
- Updates to GitHub Actions workflows

## TwoSampleMR v0.5.7

(Release date: 2023-05-29)

- Move car package to Suggests to allow TwoSampleMR to install on R
  between versions 4.0.0 and 4.1.0
- In DESCRIPTION use pkgdepends syntax for MRPRESSO package due its
  repository name being different to the package name so that installing
  TwoSampleMR under pak continues to work
- Various minor code tweaks to fix 2 R CMD check notes
- Add Cairo package to Suggests list (thanks
  [@hdraisma](https://github.com/hdraisma))
- Fix error in outcome data vignette (thanks
  [@hdraisma](https://github.com/hdraisma))
- Some p-values that should have been ~0 were being stored as 1 in the
  elasticsearch database. This has now been fixed and those datasets
  have been clumped again to re-define the tophits. A full list of
  affected GWAS is available here:
  <https://github.com/MRCIEU/opengwas-infpval-fix>
- Updated steiger filtering to use effective sample size for case
  control studies (thanks to [@niekverw](https://github.com/niekverw))
- Fixed issue with tri-allelic SNPs in harmonisation. Credit to Clare
  Horscroft ([@chorscroft](https://github.com/chorscroft)) for spotting
  the error and fixing
- Fixed an issue with experimental version of local multivariable MR
  method. Credit to Mischa Lundberg
  ([@MischaLundberg](https://github.com/MischaLundberg)).
- Catching edge cases for retrieving sample size meta data
- Updating default rsq estimation function to use beta and standard
  error instead of p-value, should improve numerical stability
- Allow chr and pos to be read in from local summary data files
- When reading in local data without p-values, editing the inferred
  p-value method to be two-sided
- All images in the vignettes (and hence also in the rendered pkgdown
  website) now have accompanying alt text descriptions
- The accompanying website for the package now uses Bootstrap 5, which
  means a search facility is enabled
- The NAMESPACE has been simplified, hence the package load time is very
  slightly improved

## TwoSampleMR v0.5.6

(Release date: 2021-03-25)

- Fix to scatter plot (thanks to Yossi Farjoun
  [@yfarjoun](https://github.com/yfarjoun))
- Update to mr.raps parameters (thanks to Qingyuan Zhao
  [@qingyuanzhao](https://github.com/qingyuanzhao))
- Bug fix to MVMR (thanks to Conor Judge
  [@conorjudge](https://github.com/conorjudge))
- Fix to harmonise_data (thanks to Leland Taylor
  [@letaylor](https://github.com/letaylor))
- Documentation (thanks to
  [@jinghuazhao](https://github.com/jinghuazhao))

## TwoSampleMR v0.5.5

(Release date: 2020-08-09)

- Updating `clump_data` function to operate on outcome datasets in the
  same way as it operates on exposure datasets. Credit goes to Marina
  Vabistsevits for spotting this and suggesting a solution.
- Removing ios function, this has now moved to mr.ios package here:
  <https://github.com/universe77/mr.ios>
- Temporarily removing some studies because the reported effect allele
  may have been incorrect, will reinstate after this has been further
  investigated. A list of studies quarantined below:
  - ieu-a-756
  - ieu-a-757
  - ieu-a-758
  - ieu-a-759
  - ieu-a-760
  - ieu-a-761
  - ieu-a-762
  - ieu-a-763
  - ieu-a-764
  - ieu-a-765
  - ieu-a-766
  - ieu-a-767
  - ieu-a-768
  - ieu-a-769
  - ieu-a-770
  - ieu-a-771
  - ieu-a-772
  - ieu-a-773
  - ieu-a-774
  - ieu-a-775
  - ieu-a-776
  - ieu-a-777
  - ieu-a-778
  - ieu-a-779
  - bbj-a-64
  - bbj-a-65
  - bbj-a-66
  - bbj-a-67
  - bbj-a-68
  - bbj-a-69
  - ebi-a-GCST004364
  - ebi-a-GCST005215
  - ebi-a-GCST005216
  - ebi-a-GCST005221
  - ebi-a-GCST005222
  - ieu-a-1086
  - ieu-a-761
  - ieu-a-762
  - ieu-a-763
  - ieu-a-767
  - ieu-a-777
  - ieu-a-779

## TwoSampleMR v0.5.4

(Release date: 2020-05-10)

- All datasets now re-instated
- Added options for different populations in LD operations
- When converting to MRInput format and supplying an LD matrix, it is
  possible that multi-allelic variants will be represented differently
  on in the GWAS and the LD reference panel. Ambiguous alignments were
  not being removed, now fixed. Credit goes to Mona Almramhi for
  spotting and fixing this issue.

## TwoSampleMR v0.5.3

(Release date: 2020-04-02)

- When converting to MRInput format and supplying an LD matrix, the LD
  matrix SNP order was not matching the summary data order. Credit goes
  to Mona Almramhi for spotting and fixing this issue.
- Reinstating all datasets that were previously disabled (ukb-a, ukb-d,
  ubm-a)
- Fixed bug with mr_wrapper. Thanks to Gunn-Helen Moen for this.

## TwoSampleMR v0.5.2

(Release date: 2020-03-11)

- No longer marking LD functions as deprecated for now. Thanks to Jonas
  Bovijn for discussions on this.
- Various fixes for `R CMD check` warnings and notes.

## TwoSampleMR v0.5.1

(Release date: 2020-02-14)

- A number of datasets have been found to have issues since 0.5.0. These
  include:
- A minority of non-effect alleles being incorrect in the ieu-a batch.
  The consequence of this is harmonisation may have thrown out some SNPs
  due to harmonisation mismatches. Error arose in 0.5.0 and now fixed
- p-value issues with the ubm-a batch. This would have led to fewer
  top-hits being identified than they should have. Error arose in 0.5.0
  and currently disabled
- Effect allele frequency issues with the ukb-a batch, potentially due
  to misreported effect alleles. Error arose in 0.5.0 and currently
  disabled

## TwoSampleMR v0.5.0

(Release date: 2020-01-01)

- Major update, details here:
  <https://mrcieu.github.io/TwoSampleMR/articles/gwas2020.html>

## TwoSampleMR v0.4.26

(Release date: 2019-12-01)

- Improved precision of low p-values in steiger tests. Thanks to Hannah
  V Meyer and Tom Palmer for this.
- Improved instrument extraction for new datasets

## TwoSampleMR v0.4.25

(Release date: 2019-09-12)

- Changes in googleAuthR package break authentication. Added
  interception to install older version while this is being fixed.
  Please use
  `devtools::install_github("MarkEdmondson1234/googleAuthR@v0.8.1")`

## TwoSampleMR v0.4.24

(Release date: 2019-09-10)

- Bug found in extract_instruments when requesting non-default
  parameters. Thanks to Shantanu Bafna for pointing this out.

## TwoSampleMR v0.4.23

(Release date: 2019-08-12)

- Forcing server to `extract_instruments` when pre-computed outcomes are
  not present by default. The old behaviour is still possible by setting
  `extract_instruments(force_server_if_empty=FALSE)`

## TwoSampleMR v0.4.22

(Release date: 2019-02-22)

- Changing default API address in preparation for moving to version
  0.5.0 which will use the new API

## TwoSampleMR v0.4.21

(Release date: 2019-02-19)

- Updated mixture of experts

## TwoSampleMR v0.4.20

(Release date: 2019-01-31)

- The harmonise function now returns a summary of the harmonisation
  procedure e.g. number of SNPs removed etc. Access via attr(obj, “log”)
- Note that this has been tested and shown to give the same results as
  previously but there is a chance that it might lead to slightly
  different behaviour. Please install the previous version if you would
  prefer to avoid possibilities of changed behaviour -
  devtools::install_github(“<MRCIEU/TwoSampleMR@0.4.18>”)

## TwoSampleMR v0.4.19

(Release date: 2019-01-31)

- Fixed a bug in mr_heterogeneity that would have impacted a minority of
  cases. If the method list was being specified then the order of the
  results didn’t always match the method (MR Egger and IVW were mixed
  up). This did not affect default usage. Thanks to Anna Guyatt for
  pointing this out.
- Added index of suspicion functionality, and penalised mode estimator
- Added transformation function to scale effect estimate units to SD
  scale
- Starting to write change log again!

## TwoSampleMR v0.4.18

(Release date: 2018-12-03)

- Improved performance of harmonisation

## TwoSampleMR v0.4.17

(Release date: 2018-12-03)

- Added facility to harmonise indels
- Documentation and options added to multivariable MR

## TwoSampleMR v0.3.4

(Release date: 2017-11-30)

- Moving over to elastic search database so the request batching is
  changing from 50 SNPs per chunk to 10000. This can be modified through
  extract_outcome_data(splitsize=50)
- Changing harmonise_data behaviour - now does not discard the bad SNPs
  but retains them with the mr_keep column indicating whether or not
  they will be used by the mr analysis functions
- Fixed issue with oauth token
- Updated scatter plot to register the mr_keep column.

## TwoSampleMR v0.3.3

(Release date: 2017-11-23)

- Fixed bug in singlesnp and leaveoneout analyses

## TwoSampleMR v0.3.2

(Release date: 2017-11-22)

- Added function to check for latest version on package load

## TwoSampleMR v0.3.1

(Release date: 2017-11-22)

- One of the external packages that TwoSampleMR depends upon had
  changed, making the authorisation behaviour change. The authorisation
  was timing out after an hour and it was not refreshing after its
  timeout. This has now been fixed - the authorisation token will
  refresh after an hour.
- The authorisation token used to be stored in a hidden file called
  .httr-oauth. This has now been changed - it will be stored in a
  visible file called ‘mrbase.oauth’.
