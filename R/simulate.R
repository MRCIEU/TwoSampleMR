
makePhen <- function(effs, indep, vy=1, vx=rep(1, length(effs)), my=0)
{
	if(is.null(dim(indep))) indep <- cbind(indep)
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	sc <- sum(cors^2)
	if(sc >= 1)
	{
		print(sc)
		stop("effects explain more than 100% of variance")
	}
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
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p, n=n
	))
}

gwas <- function(y, g)
{
	out <- matrix(0, ncol(g), 6)
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
		samplesize.exposure=gwasx$n,
		samplesize.outcome=gwasy$n,
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

make_effs <- function(ninst1, var_g1u=0, var_g1x, var_g1y=0, mu_g1y=0, var_xy, var_ux=0, var_uy=0, ninst2=0, var_g2u=0, var_g2x=0, var_g2y=0, mu_g2x=0, ninstu=0, var_guu=0)
{
	# 1 SNP influences one confounder
	eff_g1u <- chooseEffects(ninst1, var_g1u)
	eff_g1x <- chooseEffects(ninst1, var_g1x)
	eff_g1y <- chooseEffects(ninst1, var_g1y, mua=mu_g1y)

	eff_g2u <- chooseEffects(ninst2, var_g2u)
	eff_g2x <- chooseEffects(ninst2, var_g2x, mua=mu_g2x)
	eff_g2y <- chooseEffects(ninst2, var_g2y)

	eff_guu <- chooseEffects(ninstu, var_guu)

	eff_xy <- var_xy
	eff_ux <- var_ux
	eff_uy <- var_uy

	return(list(
		eff_g1u = eff_g1u,
		eff_g1x = eff_g1x,
		eff_g1y = eff_g1y,
		eff_g2u = eff_g2u,
		eff_g2x = eff_g2x,
		eff_g2y = eff_g2y,
		eff_ux = eff_ux,
		eff_uy = eff_uy,
		eff_xy = eff_xy,
		eff_guu = eff_guu
	))
}

make_pop <- function(effs, nid)
{
	ninst1 <- length(effs$eff_g1x)
	ninst2 <- length(effs$eff_g2x)
	ninstu <- length(effs$eff_guu)
	G1 <- matrix(rbinom(nid * ninst1, 2, 0.5), nid, ninst1)
	G2 <- matrix(rbinom(nid * ninst2, 2, 0.5), nid, ninst2)
	Gu <- matrix(rbinom(nid * ninstu, 2, 0.5), nid, ninstu)

	u <- makePhen(c(effs$eff_g1u, effs$eff_g2u, effs$eff_guu), cbind(G1, G2, Gu))
	x <- makePhen(c(effs$eff_g1x, effs$eff_g2x, effs$eff_ux), cbind(G1, G2, u))
	y <- makePhen(c(effs$eff_xy, effs$eff_uy, effs$eff_g1y, effs$eff_g2y), cbind(x, u, G1, G2))
	return(list(
		x=x,
		y=y,
		u=u,
		G1=G1,
		G2=G2,
		Gu=Gu
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
	mod1 <- coefficients(summary(lm(pop2$y ~ prs)))
	mod2 <- coefficients(summary(lm(pop2$y ~ wprs)))
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
