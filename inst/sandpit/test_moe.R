library(devtools)
load_all()

a <- extract_instruments(2)
b <- extract_outcome_data(a$SNP, 1001)
dat <- harmonise_data(a,b)
load("~/Dropbox/rf.rdata")
out <- mr_moe(dat, rf)

d <- dat_to_MRInput(dat)[[1]]

MendelianRandomization::mr_ivw(d, model="fixed")@StdError
o <- summary(lm(beta.outcome ~ 0 + beta.exposure, weight=1/se.outcome^2, data=dat))
o$sigma
o$coefficients[1,2] / o$sigma

o <- summary(lm(beta.outcome ~ 0 + beta.exposure, weight=1/se.outcome^2, data=dat))
o$sigma
o$coefficients[1,2] / o$sigma

library(magrittr)
mr_rucker(dat)[[1]]$rucker
