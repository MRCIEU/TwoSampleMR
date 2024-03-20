# TwoSampleMR v0.5.10

(Release date: 2024-02-20)

* Added `bfile` and `plink_bin` arguments to `clump_data()`
* Improvements to file reading and dataset formatting capabilities of `mv_extract_exposures_local()` to create the multiple exposure dataset

# TwoSampleMR v0.5.9

(Release date: 2024-02-01)

* Fixed a minor issue in `dat_to_RadialMR()`
* Minor improvements to `make_dat()` default arguments and helpfile
* Minor improvements to package tests
* Amendments to GitHub Actions workflows
* Updated several URLs which had changed

# TwoSampleMR v0.5.8

(Release date: 2023-11-16)

* Improved speed of harmonisation using data.table functions (thanks @nicksunderland)
* Updated URL to R-CMD-check README badge
* Updates to GitHub Actions workflows

TwoSampleMR v0.5.7
==============

(Release date: 2023-05-29)

* Move car package to Suggests to allow TwoSampleMR to install on R between versions 4.0.0 and 4.1.0
* In DESCRIPTION use pkgdepends syntax for MRPRESSO package due its repository name being different to the package name so that installing TwoSampleMR under pak continues to work
* Various minor code tweaks to fix 2 R CMD check notes
* Add Cairo package to Suggests list (thanks @hdraisma)
* Fix error in outcome data vignette (thanks @hdraisma)
* Some p-values that should have been ~0 were being stored as 1 in the elasticsearch database. This has now been fixed and those datasets have been clumped again to re-define the tophits. A full list of affected GWAS is available here: https://github.com/MRCIEU/opengwas-infpval-fix
* Updated steiger filtering to use effective sample size for case control studies (thanks to @niekverw)
* Fixed issue with tri-allelic SNPs in harmonisation. Credit to Clare Horscroft (@chorscroft) for spotting the error and fixing
* Fixed an issue with experimental version of local multivariable MR method. Credit to Mischa Lundberg (@MischaLundberg).
* Catching edge cases for retrieving sample size meta data
* Updating default rsq estimation function to use beta and standard error instead of p-value, should improve numerical stability
* Allow chr and pos to be read in from local summary data files
* When reading in local data without p-values, editing the inferred p-value method to be two-sided
* All images in the vignettes (and hence also in the rendered pkgdown website) now have accompanying alt text descriptions
* The accompanying website for the package now uses Bootstrap 5, which means a search facility is enabled
* The NAMESPACE has been simplified, hence the package load time is very slightly improved

TwoSampleMR v0.5.6
==============

(Release date: 2021-03-25)

* Fix to scatter plot (thanks to Yossi Farjoun @yfarjoun)
* Update to mr.raps parameters (thanks to Qingyuan Zhao @qingyuanzhao)
* Bug fix to MVMR (thanks to Conor Judge @conorjudge)
* Fix to harmonise_data (thanks to Leland Taylor @letaylor)
* Documentation (thanks to @jinghuazhao)

TwoSampleMR v0.5.5
==============

(Release date: 2020-08-09)

* Updating `clump_data` function to operate on outcome datasets in the same way as it operates on exposure datasets. Credit goes to Marina Vabistsevits for spotting this and suggesting a solution.
* Removing ios function, this has now moved to mr.ios package here: https://github.com/universe77/mr.ios
* Temporarily removing some studies because the reported effect allele may have been incorrect, will reinstate after this has been further investigated. A list of studies quarantined below:
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

TwoSampleMR v0.5.4
==============

(Release date: 2020-05-10)

* All datasets now re-instated
* Added options for different populations in LD operations
* When converting to MRInput format and supplying an LD matrix, it is possible that multi-allelic variants will be represented differently on in the GWAS and the LD reference panel. Ambiguous alignments were not being removed, now fixed. Credit goes to Mona Almramhi for spotting and fixing this issue.

TwoSampleMR v0.5.3
==============

(Release date: 2020-04-02)

* When converting to MRInput format and supplying an LD matrix, the LD matrix SNP order was not matching the summary data order. Credit goes to Mona Almramhi for spotting and fixing this issue.
* Reinstating all datasets that were previously disabled (ukb-a, ukb-d, ubm-a)
* Fixed bug with mr_wrapper. Thanks to Gunn-Helen Moen for this.

TwoSampleMR v0.5.2
==============

(Release date: 2020-03-11)

* No longer marking LD functions as deprecated for now. Thanks to Jonas Bovijn for discussions on this.
* Various fixes for `R CMD check` warnings and notes.

TwoSampleMR v0.5.1
==============

(Release date: 2020-02-14)

* A number of datasets have been found to have issues since 0.5.0. These include:
* A minority of non-effect alleles being incorrect in the ieu-a batch. The consequence of this is harmonisation may have thrown out some SNPs due to harmonisation mismatches. Error arose in 0.5.0 and now fixed
* p-value issues with the ubm-a batch. This would have led to fewer top-hits being identified than they should have. Error arose in 0.5.0 and currently disabled
* Effect allele frequency issues with the ukb-a batch, potentially due to misreported effect alleles. Error arose in 0.5.0 and currently disabled

TwoSampleMR v0.5.0
==============

(Release date: 2020-01-01)

* Major update, details here: https://mrcieu.github.io/TwoSampleMR/articles/gwas2020.html

TwoSampleMR v0.4.26
==============

(Release date: 2019-12-01)

* Improved precision of low p-values in steiger tests. Thanks to Hannah V Meyer and Tom Palmer for this.
* Improved instrument extraction for new datasets

TwoSampleMR v0.4.25
==============

(Release date: 2019-09-12)

* Changes in googleAuthR package break authentication. Added interception to install older version while this is being fixed. Please use `devtools::install_github("MarkEdmondson1234/googleAuthR@v0.8.1")`

TwoSampleMR v0.4.24
==============

(Release date: 2019-09-10)

* Bug found in extract_instruments when requesting non-default parameters. Thanks to Shantanu Bafna for pointing this out.

TwoSampleMR v0.4.23
==============

(Release date: 2019-08-12)

* Forcing server to `extract_instruments` when pre-computed outcomes are not present by default. The old behaviour is still possible by setting `extract_instruments(force_server_if_empty=FALSE)`

TwoSampleMR v0.4.22
==============

(Release date: 2019-02-22)

* Changing default API address in preparation for moving to version 0.5.0 which will use the new API

TwoSampleMR v0.4.21
==============

(Release date: 2019-02-19)

* Updated mixture of experts

TwoSampleMR v0.4.20
==============

(Release date: 2019-01-31)

* The harmonise function now returns a summary of the harmonisation procedure e.g. number of SNPs removed etc. Access via attr(obj, "log")
* Note that this has been tested and shown to give the same results as previously but there is a chance that it might lead to slightly different behaviour. Please install the previous version if you would prefer to avoid possibilities of changed behaviour - devtools::install_github("MRCIEU/TwoSampleMR@0.4.18")

TwoSampleMR v0.4.19
==============

(Release date: 2019-01-31)

* Fixed a bug in mr_heterogeneity that would have impacted a minority of cases. If the method list was being specified then the order of the results didn't always match the method (MR Egger and IVW were mixed up). This did not affect default usage. Thanks to Anna Guyatt for pointing this out.
* Added index of suspicion functionality, and penalised mode estimator
* Added transformation function to scale effect estimate units to SD scale
* Starting to write change log again!

TwoSampleMR v0.4.18
==============

(Release date: 2018-12-03)

* Improved performance of harmonisation

TwoSampleMR v0.4.17
==============

(Release date: 2018-12-03)

* Added facility to harmonise indels
* Documentation and options added to multivariable MR

TwoSampleMR v0.3.4
==============

(Release date: 2017-11-30)

* Moving over to elastic search database so the request batching is changing from 50 SNPs per chunk to 10000. This can be modified through extract_outcome_data(splitsize=50)
* Changing harmonise_data behaviour - now does not discard the bad SNPs but retains them with the mr_keep column indicating whether or not they will be used by the mr analysis functions
* Fixed issue with oauth token
* Updated scatter plot to register the mr_keep column.

TwoSampleMR v0.3.3
==============

(Release date: 2017-11-23)

* Fixed bug in singlesnp and leaveoneout analyses

TwoSampleMR v0.3.2
==============

(Release date: 2017-11-22)

* Added function to check for latest version on package load

TwoSampleMR v0.3.1
==============

(Release date: 2017-11-22)

* One of the external packages that TwoSampleMR depends upon had changed, making the authorisation behaviour change. The authorisation was timing out after an hour and it was not refreshing after its timeout. This has now been fixed - the authorisation token will refresh after an hour.

* The authorisation token used to be stored in a hidden file called .httr-oauth. This has now been changed - it will be stored in a visible file called 'mrbase.oauth'.
