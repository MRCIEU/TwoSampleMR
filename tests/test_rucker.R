source("~/repo/instrument-directionality/scripts/simulation-functions.r")
library(devtools)
load_all()


n1 <- 50000
n2 <- 50000
nsnp <- 50

# models


effs <- make_effs(ninst1=nsnp, var_xy=0.5, var_g1x=0.5, var_g1y=0, mu_g1y=0)
pop1 <- make_pop(effs, n1)
pop2 <- make_pop(effs, n2)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G1)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G1)
datA <- recode_dat(make_dat(dat1, dat2))
run_rucker(datA)
plot(beta.outcome ~ beta.exposure, datA)

with(datA, mr_mode(beta.exposure, beta.outcome, se.exposure, se.outcome))

effs <- make_effs(ninst1=nsnp, var_xy=0.5, var_g1x=0.5, var_g1y=0.1, mu_g1y=0)
pop1 <- make_pop(effs, n1)
pop2 <- make_pop(effs, n2)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G1)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G1)
datB <- recode_dat(make_dat(dat1, dat2))
run_rucker(datB)
plot(beta.outcome ~ beta.exposure, datB)

mr_pleiotropy_test(datB)

effs <- make_effs(ninst1=nsnp, var_xy=0.5, var_g1x=0.5, var_g1y=0, mu_g1y=0.1)
pop1 <- make_pop(effs, n1)
pop2 <- make_pop(effs, n2)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G1)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G1)
datC <- recode_dat(make_dat(dat1, dat2))
run_rucker(datC)
plot(beta.outcome ~ beta.exposure, datC)


effs <- make_effs(ninst1=nsnp, var_xy=0.5, var_g1x=0.5, var_g1y=0.1, mu_g1y=0.1)
pop1 <- make_pop(effs, n1)
pop2 <- make_pop(effs, n2)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G1)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G1)
datD <- recode_dat(make_dat(dat1, dat2))
run_rucker(datD)
plot(beta.outcome ~ beta.exposure, datD)


param <- expand.grid(
	nsim = c(1:50),
	nsnp = c(5,20,50),
	nid1 = 50000,
	nid2 = 50000,
	var_xy = c(0, 0.05, 0.1, 0.5),
	var_g1x = c(0, 0.025, 0.1),
	mu_g1y = c(-0.1, 0, 0.1)
)

dim(param)

out <- list()
for(i in 1:nrow(param))
{
	effs <- make_effs(ninst1=param$nsnp[i], var_xy=param$var_xy[i], var_g1x=param$var_g1x[i], mu_g1y=param$mu_g1y[i])
	pop1 <- make_pop(effs, param$nid1[i])
	pop2 <- make_pop(effs, param$nid2[i])
	dat1 <- get_effs(pop1$x, pop1$y, pop1$G1)
	dat2 <- get_effs(pop2$x, pop2$y, pop2$G1)
	dat <- recode_dat(make_dat(dat1, dat2))
	out <- run_rucker(dat)
	out$param <- param[i,]
	res[[i]] <- out
}



# https://www.ncbi.nlm.nih.gov/pubmed/21473747
# http://amstat.tandfonline.com/doi/abs/10.1080/00031305.2016.1165735


