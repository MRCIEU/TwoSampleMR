# Perform LD clumping on SNP data

Uses PLINK clumping method, where SNPs in LD within a particular window
will be pruned. The SNP with the lowest p-value is retained.

## Usage

``` r
clump_data(
  dat,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR",
  bfile = NULL,
  plink_bin = NULL
)
```

## Arguments

- dat:

  Output from
  [`format_data()`](https://mrcieu.github.io/TwoSampleMR/reference/format_data.md).
  Must have a SNP name column (`SNP`), SNP chromosome column
  (`chr_name`), SNP position column (`chrom_start`). If `id.exposure` or
  `pval.exposure` not present they will be generated.

- clump_kb:

  Clumping window, default is `10000`.

- clump_r2:

  Clumping r2 cutoff. Note that this default value has recently changed
  from `0.01` to `0.001`.

- clump_p1:

  Clumping sig level for index SNPs, default is `1`.

- clump_p2:

  Clumping sig level for secondary SNPs, default is `1`.

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

Data frame

## Details

This function interacts with the OpenGWAS API, which houses LD reference
panels for the 5 super-populations in the 1000 genomes reference panel.
It includes only bi-allelic SNPs with MAF \> 0.01, so it's quite
possible that a variant you want to include in the clumping process will
be absent. If it is absent, it will be automatically excluded from the
results.

You can check if your variants are present in the LD reference panel
using
[`ieugwasr::ld_reflookup()`](https://mrcieu.github.io/ieugwasr/reference/ld_reflookup.html).

This function does put load on the OpenGWAS servers, which makes life
more difficult for other users. We have implemented a method and made
available the LD reference panels to perform clumping locally, see
[`ieugwasr::ld_clump()`](https://mrcieu.github.io/ieugwasr/reference/ld_clump.html)
and related vignettes for details.
