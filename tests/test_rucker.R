library(devtools)
load_all()


n1 <- 50000
n2 <- 50000
nsnp <- 50
bxy <- 0.5
bgx <- runif(nsnp, 0.5, 1)
bgy <- rnorm(nsnp)


effs <- make_effs(nsnp, var_gu=0, var_ux=0, var_uy=0, var_xy=0.5, var_gx=0.5, var_gy=0, mu_gy=0)
pop1 <- make_pop(nsnp, n1, effs)
pop2 <- make_pop(nsnp, n2, effs)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G)
datA <- recode_dat(make_dat(dat1, dat2))
resA <- analyse_simulation(datA, pop2)



effs <- make_effs(ninst1 = nsnp, var_xy = 0.4, var_g1x = 0.3, ninst2 = nsnp, var_g2y = 0.2)
pop1 <- make_pop(effs, n1)
dat1 <- recode_dat(get_effs(pop1$y, pop1$x, pop1$G1))
dat2 <- recode_dat(get_effs(pop1$y, pop1$x, pop1$G2))
dat <- rbind(dat1, dat2)
dats <- subset(dat, pval.exposure < 5e-8)

mr(dat1)
mr(dat2)
mr(dat)
mr(dats)


mr_scatter_plot(mr(dat1), dat1)[[1]]
mr_scatter_plot(mr(dat2), dat2)[[1]]
mr_scatter_plot(mr(dat), dat)[[1]]
mr_scatter_plot(mr(dats), dats)[[1]]


mr_scatter_plot(mr(dats), dats)[[1]] + geom_point(data=dats[index,], colour="red")

dats <- direct_vs_reverse(dats)



effs <- make_effs(ninst1 = nsnp, var_xy = 0.4, var_g1x = 0.3, ninst2 = nsnp, var_g2y = 0.2, var_g1y = 0.05, var_g2x = 0.05)
pop1 <- make_pop(effs, n1)
dat1 <- recode_dat(get_effs(pop1$y, pop1$x, pop1$G1))
dat1$inst <- "x"
dat2 <- recode_dat(get_effs(pop1$y, pop1$x, pop1$G2))
dat2$inst <- "y"
dat <- rbind(dat1, dat2)
dats <- subset(dat, pval.exposure < 5e-8)
dats <- direct_vs_reverse(dats)

mr_scatter_plot(mr(dats), dats)[[1]] + geom_point(data=dats[index,], colour="red")

mr_scatter_plot(mr(dats), dats)[[1]] + geom_point(data=dats, aes(colour=inst))
mr(dats)
mr(dat[dat$pval.exposure < 5e-8, ])






effs <- make_effs(ninst1 = nsnp, var_xy = 0.4, var_g1x = 0.3, ninst2 = nsnp, var_g2y = 0.2, ninstu = 50, var_guu = 0.4, var_ux=0.2, var_uy=0.2)
pop1 <- make_pop(effs, n1)

dat1 <- recode_dat(get_effs(pop1$y, pop1$x, pop1$G1))
dat2 <- recode_dat(get_effs(pop1$y, pop1$x, pop1$G2))
dat3 <- recode_dat(get_effs(pop1$y, pop1$x, pop1$Gu))
dat1$inst <- "x"
dat2$inst <- "y"
dat3$inst <- "u"

dat <- rbind(dat1, dat2, dat3)
dats <- subset(dat, pval.exposure < 5e-8)
table(dats$inst)

table(dats$inst, dats$mr_keep)

mr_scatter_plot(mr(dats), dats)[[1]] + geom_point(data=dats, aes(colour=inst))
mr(dats[dats$inst == "u",])

mr(dat1)
mr(dat[dat$pval.exposure < 5e-8,])
mr_scatter_plot(mr(dat), dat[dat$mr_keep,])[[1]] + geom_point(data=dats, aes(colour=inst))


d1 <- gwas(pop1$y, cbind(pop1$G1, pop1$G2, pop1$Gu))
d2 <- gwas(pop1$x, cbind(pop1$G1, pop1$G2, pop1$Gu))
d3 <- gwas(pop1$u, cbind(pop1$G1, pop1$G2, pop1$Gu))

d <- find_invalid_instruments(d1, d2, d3)
keep <- (d$sig_inst & !d$remove_reverse & !d$remove_confounder)

table(dat$inst, keep)


# models

effs <- make_effs(nsnp, var_gu=0, var_ux=0, var_uy=0, var_xy=0.5, var_gx=0.5, var_gy=0.1, mu_gy=0)
pop1 <- make_pop(nsnp, n1, effs)
pop2 <- make_pop(nsnp, n2, effs)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G)
datB <- recode_dat(make_dat(dat1, dat2))
resB <- analyse_simulation(datB, pop2)


effs <- make_effs(nsnp, var_gu=0, var_ux=0, var_uy=0, var_xy=0.5, var_gx=0.5, var_gy=0, mu_gy=0.1)
pop1 <- make_pop(nsnp, n1, effs)
pop2 <- make_pop(nsnp, n2, effs)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G)
datC <- recode_dat(make_dat(dat1, dat2))
resC <- analyse_simulation(datC, pop2)


effs <- make_effs(nsnp, var_gu=0, var_ux=0, var_uy=0, var_xy=0.5, var_gx=0.5, var_gy=0.1, mu_gy=0.1)
pop1 <- make_pop(nsnp, n1, effs)
pop2 <- make_pop(nsnp, n2, effs)
dat1 <- get_effs(pop1$x, pop1$y, pop1$G)
dat2 <- get_effs(pop2$x, pop2$y, pop2$G)
datD <- recode_dat(make_dat(dat1, dat2))
resD <- analyse_simulation(datD, pop2)


plot(beta.outcome ~ beta.exposure, datA)
plot(beta.outcome ~ beta.exposure, datB)
plot(beta.outcome ~ beta.exposure, datC)
plot(beta.outcome ~ beta.exposure, datD)


rucker(datA)
rucker(datB)
rucker(datC)
rucker(datD)




rucker_ivw_fe <- function(b_exp, b_out, se_exp, se_out)
{
	y <- b_out / se_out
	x <- b_exp / se_out
	mod <- summary(lm(y ~ 0 + x))
	b <- coefficients(mod)[1,1]
	se <- coefficients(mod)[1,2]
	pval <- coefficients(mod)[1,4]
	Q <- sum((y - x*b)^2)
	Q_df <- length(b_exp) - 1
	Q_pval <- pchisq(Q, Q_df, low = FALSE)
	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}

rucker_ivw_re <- function(b_exp, b_out, se_exp, se_out)
{
	out <- rucker_ivw_fe(b_exp, b_out, se_exp, se_out)
	b <- out$b
	w <- b_exp^2 / se_out^2
	phi <- out$Q / (out$nsnp - 1)
	se <- phi / sum(w)
	pval <- pt(abs(b/se), out$nsnp-1, lower.tail=FALSE) * 2
	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), Q = out$Q, Q_df = out$Q_df, Q_pval = out$Q_pval))
}

run_rucker <- function(dat)
{
	nsnp <- nrow(dat)
	b_exp <- dat$beta.exposure
	b_out <- dat$beta.outcome
	se_exp <- dat$se.exposure
	se_out <- dat$se.outcome
	w <- b_exp^2 / se_out^2


	# IVW FE
	y <- b_out / se_out
	x <- b_exp / se_out
	mod_ivw_fe <- summary(lm(y ~ 0 + x))
	b_ivw_fe <- coefficients(mod_ivw_fe)[1,1]
	se_ivw_fe <- coefficients(mod_ivw_fe)[1,2]
	pval_ivw_fe <- coefficients(mod_ivw_fe)[1,4]
	Q_ivw <- sum((y - x*b_ivw_fe)^2)
	Q_df_ivw <- length(b_exp) - 1
	Q_pval_ivw <- pchisq(Q_ivw, Q_df_ivw, low = FALSE)

	# IVW MRE
	b_ivw_re <- b_ivw_fe
	phi_ivw <- Q_ivw / (nsnp - 1)
	se_ivw_re <- phi_ivw / sum(w)
	pval_ivw_re <- pt(abs(b_ivw_re/se_ivw_re), nsnp-1, lower.tail=FALSE) * 2


	# Egger MRE
	i <- 1 / se_out
	mod_egger <- summary(lm(y ~ 0 + i + x))
	b_egger_fe <- coefficients(mod_egger)[2,1]
	se_egger_fe <- coefficients(mod_egger)[2,1]

	se_egger_fe


	# 


	dat <- data.frame(
		method = c("IVW fixed effects", "IVW random effects"),
		b = c(b_ivw_fe, b_ivw_re),
		se = c(se_ivw_fe, se_ivw_re),
		pval = c(pval_ivw_fe, pval_ivw_re),
		nsnp = nsnp
	)

	return(dat)
}

run_rucker(datB)





# https://www.ncbi.nlm.nih.gov/pubmed/21473747
# http://amstat.tandfonline.com/doi/abs/10.1080/00031305.2016.1165735


