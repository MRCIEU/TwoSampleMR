source("~/repo/instrument-directionality/scripts/simulation-functions.r")
library(devtools)
load_all()


## Minimal example
devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
bp <- read.table(system.file(package="TwoSampleMR", "data/DebbieData_2.txt"))
names(bp) <- c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
bp$mr_keep <- TRUE
bp$id.exposure <- bp$id.outcome <- bp$exposure <- bp$outcome <- 1
bp$SNP <- paste0("SNP", 1:nrow(bp))

r1 <- mr_rucker(bp)
r2 <- mr_rucker_cooksdistance(bp)

r1$selected
r2$selected

a1 <- mr_rucker(bp)
a <- mr_rucker_cooksdistance(bp)
b <- mr_rucker_bootstrap(bp)
bp2 <- subset(bp, ! SNP %in% a$removed_snps)
c <- mr_rucker_bootstrap(bp2)

pdf("rucker_bootstrap.pdf", width=10, height=10)
gridExtra::grid.arrange(
	b$q_plot + labs(title="All variants Rucker bootstrap"),
	c$q_plot + labs(title="Excluding Cooks D > 4/n Rucker bootstrap"),
	b$e_plot + labs(title="All variants Rucker bootstrap"),
	c$e_plot + labs(title="Excluding Cooks D > 4/n Rucker bootstrap"),
	ncol=2
)
dev.off()

# permutations
res1 <- array(0, 100)
res2 <- array(0, 100)
res3 <- array(0, 100)
for(i in 1:100)
{
	message(i)
	bp3 <- bp
	index <- sample(1:nrow(bp3), replace=FALSE)
	bp3$beta.outcome <- bp3$beta.outcome[index]
	bp3$se.outcome <- bp3$se.outcome[index]
	r1 <- mr_rucker_cooksdistance(bp3)
	r2 <- rucker_bootstrap(bp3)
	r3 <- rucker_bootstrap(subset(bp3, ! SNP %in% r1$removed_snps))
	res1[i] <- r1$selected$P
	res2[i] <- r2$res$P[5]
	res3[i] <- r3$res$P[5]
}

res4 <- array(0, 100)
for(i in 1:100)
{
	message(i)
	bp3 <- bp
	index <- sample(1:nrow(bp3), replace=FALSE)
	bp3$beta.outcome <- bp3$beta.outcome[index]
	bp3$se.outcome <- bp3$se.outcome[index]
	res4[i] <- with(bp3, mr_ivw(beta.exposure, beta.outcome, se.exposure, se.outcome))$pval
}


a <- influence.measures(mod)
cooks.distance(mod) > 4/nrow(bp)
max(cooks.distance(mod))

a <- mr_all(bp)

mr_scatter_plot(mr(bp), bp)
mr_scatter_plot(mr(bp2), bp2)

a$rucker$rucker

n1 <- 50000
n2 <- 50000
nsnp <- 50

# models


effs <- make_effs(ninst1=nsnp, var_xy=0.2, var_g1x=0.5, var_g1y=0, mu_g1y=0)
pop1 <- make_pop(effs, n1)
pop2 <- make_pop(effs, n2)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G1)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G1)
datA <- recode_dat(make_dat(dat1, dat2))
a <- mr_all(datA)
datAp <- subset(datA, pval.exposure < 5e-8)
dim(datAp)

b <- mr_all(datAp)
plot(beta.outcome ~ beta.exposure, datA)
plot(beta.outcome ~ beta.exposure, datAp)

a$rucker$q_plot
b$rucker$q_plot

a$rucker$e_plot
a$res

with(datA, mr_mode(beta.exposure, beta.outcome, se.exposure, se.outcome))

effs <- make_effs(ninst1=nsnp, var_xy=0.05, var_g1x=0.2, var_g1y=0.01, mu_g1y=0)
pop1 <- make_pop(effs, n1)
pop2 <- make_pop(effs, n2)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G1)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G1)
datB <- recode_dat(make_dat(dat1, dat2))
b <- mr_all(datB)
b$rucker$q_plot
b$rucker$e_plot

r <- rucker_bootstrap(datB)
plot(beta.outcome ~ beta.exposure, datB)
r$q_plot

mr_pleiotropy_test(datB)

effs <- make_effs(ninst1=nsnp, var_xy=0.5, var_g1x=0.5, var_g1y=0, mu_g1y=0.1)
pop1 <- make_pop(effs, n1)
pop2 <- make_pop(effs, n2)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G1)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G1)
datC <- recode_dat(make_dat(dat1, dat2))
run_rucker(datC)
plot(beta.outcome ~ beta.exposure, datC)


effs <- make_effs(ninst1=nsnp, var_xy=0.5, var_g1x=0.1, var_g1y=0.002, mu_g1y=0.1)
pop1 <- make_pop(effs, n1)
pop2 <- make_pop(effs, n2)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G1)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G1)
datD <- recode_dat(make_dat(dat1, dat2))
plot(beta.outcome ~ beta.exposure, datD)


d <- mr_all(datD)
d$rucker$q_plot

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




## Irwin's data [Williams 1987]
xi <- 1:5
yi <- c(0,2,14,19,30)    # number of mice responding to dose xi
mi <- rep(40, 5)         # number of mice exposed
summary(lmI <- glm(cbind(yi, mi -yi) ~ xi, family = binomial))
signif(cooks.distance(lmI), 3)   # ~= Ci in Table 3, p.184
(imI <- influence.measures(lmI))



for(i in 1:BootSim){
BXG = rnorm(length(BetaXG),BetaXG,seBetaXG)
BYG = rnorm(length(BetaYG),BetaYG,seBetaYG)

if(weights==1){W = BXG^2/seBetaYG^2}
if(weights==2){W = 1/(seBetaYG^2/BXG^2 + (BYG^2)*seBetaXG^2/BXG^4)}

            BIVw = BIV*sqrt(W)
            sW   = sqrt(W)
IR          = lm(BIVw ~ -1+sW);IVW[i] = IR$coef[1]
MR          = lm(BIVw ~    sW);E[i]   = MR$coef[2]

DF1      = length(BetaYG)-1
phi_IVW  = summary(IR)$sigma^2
QQ[i]    = DF1*phi_IVW
DF2      = length(BetaYG)-2
phi_E    = summary(MR)$sigma^2
QQd[i]   = DF2*phi_E

Qp       = 1-pchisq(Q,DF1)

if(QQ[i] <= qchisq(1-alpha,DF1)){Mod[i]=1}
if(QQ[i] >= qchisq(1-alpha,DF1)){Mod[i]=2}
if(QQ[i] >= qchisq(1-alpha,DF1) & QQ[i] - QQd[i] >= qchisq(1-alpha,1)){Mod[i]=3}
if(QQ[i] >= qchisq(1-alpha,DF1)& QQ[i] - QQd[i] >= qchisq(1-alpha,1)& QQd[i] >=qchisq(1-alpha,DF2)){Mod[i]=4}

}
