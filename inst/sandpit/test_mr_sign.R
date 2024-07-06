
# MR sign concordance - how often are the SNP-exposure and SNP-outcome effects in the same direction?

# Use binomial test assuming that on average should be 50% under the null hypothesis

# What is the minumum number of SNPs required to achieve p-value of 0.05?

param <- expand.grid(n = 1:100, x=0:100)
param <- subset(param, x <= n)
for (i in seq_len(nrow(param)))
{
	param$pval[i] <- binom.test(x=param$x[i], n=param$n[i], p=0.5)$p.value
}
min(param$pval)

library(ggplot2)
ggplot(param, aes(x=x, y=n)) +
geom_point(aes(colour=pval < 0.05))

library(dplyr)
group_by(param, n) %>% summarise(minp = min(pval))

# Looks like 6 is the smallest number of SNPs that can achieve p < 0.05
library(TwoSampleMR)
a <- extract_instruments(2)
b <- extract_outcome_data(a$SNP,7)
dat <- harmonise_data(a, b)
mr(dat, method_list=c("mr_ivw", "mr_sign"))
mr_scatter_plot(mr(dat, method_list=c("mr_ivw", "mr_sign")), dat)

# Could consider doing a parametric bootstrap type analysis
