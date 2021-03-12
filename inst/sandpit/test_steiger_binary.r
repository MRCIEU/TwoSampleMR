# Test binary steiger

devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)

# Extract data for BMI vs CHD
bmi <- extract_instruments(2)
chd <- extract_outcome_data(bmi$SNP, 7)
dat <- harmonise_data(bmi, chd)

# Calculate variance explained in BMI (using pvals and sample size)
dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)

# Calculate variance explained in CHD (using log OR and allele frequencies)
dat$r.outcome <- get_r_from_lor(dat$beta.outcome, dat$eaf.outcome, dat$ncase.outcome, dat$ncontrol.outcome, prevalence=0.4)

# Test causal direction
mr_steiger2(dat$r.exposure, dat$r.outcome, dat$samplesize.exposure, dat$samplesize.outcome)
