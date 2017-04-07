library(tidyverse)


makePhen <- function(effs, indep, vy=1, vx=rep(1, length(effs)), my=0)
{
	if(is.null(dim(indep))) indep <- cbind(indep)
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	stopifnot(sum(cors^2) <= 1)
	cors <- c(cors, sqrt(1-sum(cors^2)))
	indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
	y <- drop(scale(rowSums(indep)) * sqrt(vy)) + my
	return(y)
}

chooseEffects <- function(nsnp, totvar, sqrt=TRUE, mua=0)
{
	eff <- rnorm(nsnp)
	eff <- sign(eff) * eff^2
	aeff <- abs(eff)
	sc <- sum(aeff) / totvar
	out <- eff / sc
	if(sqrt)
	{
		out <- sqrt(abs(out)) * sign(out)
	}
	return(out + mua)
}

fastAssoc <- function(y, x)
{
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	# fitted <- ahat + x * bhat
	# residuals <- y - fitted
	# SSR <- sum((residuals - mean(residuals))^2)
	# SSF <- sum((fitted - mean(fitted))^2)

	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
	se <- abs(bhat / tval)

	# Fval <- (SSF) / (SSR/(n-2))
	# pval <- pf(Fval, 1, n-2, lowe=F)
	p <- pf(fval, 1, n-2, lowe=F)
	return(list(
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p
	))
}

gwas <- function(y, g)
{
	out <- matrix(0, ncol(g), 5)
	for(i in 1:ncol(g))
	{
		o <- fastAssoc(y, g[,i])
		out[i, ] <- unlist(o)
	}
	out <- as.data.frame(out)
	names(out) <- names(o)
	return(out)
}

get_effs <- function(x, y, g)
{
	gwasx <- gwas(x, g)
	gwasy <- gwas(y, g)

	dat <- data.frame(
		exposure="X",
		id.exposure="X",
		outcome="Y",
		id.outcome="Y",
		beta.exposure=gwasx$bhat,
		beta.outcome=gwasy$bhat,
		se.exposure=gwasx$se,
		se.outcome=gwasy$se,
		pval.exposure=gwasx$pval,
		pval.outcome=gwasy$pval,
		mr_keep=TRUE
	)
	return(dat)
}

dat_sign <- function(dat)
{
	sign0 <- function(x) {
		x[x == 0] <- 1
		return(sign(x))
	}
	index <- sign0(dat$beta.exposure) == -1
	dat$beta.exposure <- abs(dat$beta.exposure)
	dat$beta.outcome[index] <- dat$beta.outcome[index] * -1
	return(dat)
}

recode_dat <- function(dat)
{
	a <- lm(beta.outcome ~ beta.exposure, dat)$coefficients[1]
	index <- dat$beta.exposure < 0
	dat$beta.exposure[index] <- dat$beta.exposure[index] * -1
	dat$beta.outcome[index] <- dat$beta.outcome[index] * -1 + 2 * a
	dat$index <- index
	return(dat)
}

make_effs <- function(ninst, var_gu, var_ux, var_uy, var_xy, var_gx, var_gy, mu_gy=0)
{
	# 1 SNP influences one confounder
	eff_gu <- rep(var_gu, ninst)
	eff_ux <- chooseEffects(ninst, var_ux)
	eff_uy <- chooseEffects(ninst, var_uy)
	eff_xy <- var_xy
	eff_gx <- chooseEffects(ninst, var_gx)
	eff_gy <- chooseEffects(ninst, var_gy, mua=mu_gy)
	return(list(
		eff_gu = eff_gu,
		eff_ux = eff_ux,
		eff_uy = eff_uy,
		eff_xy = eff_xy,
		eff_gx = eff_gx,
		eff_gy = eff_gy
	))
}

make_pop <- function(ninst, nid, effs)
{

	G <- matrix(rbinom(nid * ninst, 2, 0.5), nid, ninst)
	U <- matrix(0, nid, ninst)
	for(i in 1:ninst)
	{
		U[,i] <- makePhen(effs$eff_gu[i], G[,i])
	}

	x <- makePhen(c(effs$eff_gx, effs$eff_ux), cbind(G, U))
	y <- makePhen(c(effs$eff_xy, effs$eff_uy, effs$eff_gy), cbind(x, U, G))
	return(list(
		x=x,
		y=y,
		U=U,
		G=G
	))
}

make_dat <- function(exposure, outcome)
{
	dat <- cbind(
		exposure[,grepl("exposure", names(exposure))],
		outcome[,grepl("outcome", names(outcome))]
	)
	dat$mr_keep <- TRUE
	return(dat)
}

analyse_simulation <- function(dat, pop2)
{
	prs <- pop2$G %*% sign(dat$beta.exposure)
	wprs <- pop2$G %*% dat$beta.exposure
	mod1 <- summary(lm(pop2$y ~ prs)) %>% coefficients
	mod2 <- summary(lm(pop2$y ~ wprs)) %>% coefficients
	scoreres <- data.frame(
		id.exposure = "X", id.outcome = "Y", exposure = "X", outcome = "Y",
		method = c("Basic PRS", "Weighted PRS"),
		nsnp = nrow(dat),
		b = c(mod1[2,1], mod2[2,1]),
		se = c(mod1[2,2], mod2[2,2]),
		pval = c(mod1[2,4], mod2[2,4])
	)
	ivres <- mr(dat, method_list=mr_method_list()$obj[c(10, 6, 8)])
	res <- rbind(ivres, scoreres)
	return(res)
}
