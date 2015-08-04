# Two Sample MR


## Part 1
- Read in 
- Get effects and standard errors for exposure and outcome relative to the same allele
- Get standard errors


## Part 2

### Perform 2 sample MR using
- Inverse variance weighted approach (wald ratio with standard error estimated by delta method)
- ML approach
- Heterogeneity tests

### Sensitivity analysis
- Egger regression
- Weighted median function

### Plots
- Funnel plot
- Forest plot
- Scatter plots (slopes with and without intercept)



## To add:

- Allow SNPs to be chr_pos or rsid


## To do:

- Setup 1000 genomes data for LD pruning
- Code by chr:pos
- Keep only Europeans
- Make sure no repeated SNP IDs


## Wish list

- multiple ethnicities
- LD 



## Using the database


    # Read in exposure data
    # This is required to have the following columns with headers:
    # SNP, beta, se, eaf, effect_allele, other_allele
    exposure <- read_exposure_data("filename.txt", "trait name")

    # (Recommended) Prune on LD
    exposure <- ld_pruning(exposure)

    # Identify the available GWAS summary stats
    # Select the ones you want to test against
    ss_list <- summary_stats_list()

    # Perform lookup
    # This extracts the summary stats for your SNPs and harmonises the effect sizes
    lookup <- extract_outcome_data(exposure)

    # Perform MR analysis
    mr_res <- run_mr(lookup)
    summary(mr_res)

    # 
