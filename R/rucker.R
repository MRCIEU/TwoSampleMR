library(simex)

walds <- function(b_exp, b_out, se_exp, se_out)
{
	bj <- b_out / b_exp
	wj_1 <- b_exp^2 / se_out^2
	wj_2 <- 1 / (se_out^2 / b_exp^2 + b_out^2 * se_exp^2 / b_exp^4)
	b_ivw <- sum(wj_2 * bj) / sum(wj_2)
	wj_3 <- 1 / (se_out^2 / b_exp^2 + b_ivw^2 * se_exp^2 / b_exp^2)
	l <- list(
		"First order" = list(bj=bj, wj=wj_1),
		"Second order" = list(bj=bj, wj=wj_2),
		"Modified second order" = list(bj=bj, wj=wj_3)
	)
	return(l)
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

Isq <- function(y,s)
{
	k <- length(y)
	w <- 1/s^2
	sum.w <- sum(w)
	mu.hat <- sum(y*w)/sum.w
	Q <- sum(w*(y-mu.hat)^2)
	Isq <- (Q - (k-1))/Q
	Isq <- max(0,Isq)
	return(Isq)
}


# mr_ivw_fe <- function (b_exp, b_out, se_exp, se_out)
# {
# 	lapply(x, function(x){
# 		wj <- x$wj
# 		bj <- x$bj
# 		b <- sum(wj * bj) / sum(wj)
# 		se <- sqrt(1 / sum(wj))
# 		pval <- 2 * pnorm(abs(b/se), low = FALSE)
# 		Q <- sum(wj * (bj - b)^2)
# 		Q_df <- length(b_exp) - 1
# 		Q_pval <- pchisq(Q, Q_df, low = FALSE)
# 		data.frame(
# 			b1=b, se1=se, pval1=pval,
# 		)
# 		return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
# 	})
# }



mr_ivw_fe <- function (b_exp, b_out, se_exp, se_out)
{
	bj <- b_out / b_exp
	wj <- b_exp^2 / se_out^2
	b <- sum(wj * bj) / sum(wj)
	se <- sqrt(1 / sum(wj))
	pval <- 2 * pnorm(abs(b/se), low = FALSE)
	Q <- sum(wj * (bj - b)^2)
	Q_df <- length(b_exp) - 1
	Q_pval <- pchisq(Q, Q_df, low = FALSE)
	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}


mr_ivw_are <- function (b_exp, b_out, se_exp, se_out, parameters)
{
	if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2) {
		return(list(b = NA, se = NA, pval = NA, nsnp = NA, Q = NA, Q_df = NA, Q_pval = NA, tau_sq=NA))
	}

	bj <- b_out/b_exp
	sej <- sqrt((se_out^2/b_exp^2))

	res <- meta::metagen(bj, sej)

	# var(alpha)
	tau_sq <- res$tau^2

	b <- res$TE.random
	se <- res$seTE.random
	pval <- res$pval.random
	Q_pval <- pchisq(res$Q, res$df.Q, low = FALSE)
	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), Q = res$Q, Q_df = res$df.Q, Q_pval = Q_pval, tau_sq=tau_sq))
}

mr_ivw_mre <- function (b_exp, b_out, se_exp, se_out, parameters) 
{
	if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2) {
		return(list(b = NA, se = NA, pval = NA, nsnp = NA, Q = NA, Q_df = NA, Q_pval = NA))
	}

	wj <- b_exp^2 / se_out^2
	bj <- b_out / b_exp

	b <- sum(wj * bj) / sum(wj)

	Q <- sum(wj * (bj - b)^2)
	phi <- Q / (length(b_exp) - 1)

	se <- sqrt(phi / sum(wj))

	pval <- 2 * pnorm(abs(b/se), low = FALSE)

	Q_df <- length(b_exp) - 1
	Q_pval <- pchisq(Q, Q_df, low = FALSE)
	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}

# mr_egger_fe <- function (b_exp, b_out, se_exp, se_out, parameters)
# {
# 	if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3) {
# 		return(list(b = NA, se = NA, pval = NA, nsnp = NA, Q = NA, Q_df = NA, Q_pval = NA))
# 	}

# 	Z <- cbind(1 / se_out, b_exp / se_out)

# 	G <- matrix(b_out / se_out, length(b_out))
# 	p1 <- solve(t(Z) %*% Z)
# 	b <- p1 %*% t(Z) %*% G
# 	se <- sqrt(diag(p1))

# 	wj <- b_exp^2 / se_out^2
# 	bj <- b_out / b_exp

# 	Q <- sum(wj * (bj - b[1]/b_exp - b[2])^2)

# 	pval <- 2 * pnorm(abs(b/se), low = FALSE)

# 	Q_df <- length(b_exp) - 2
# 	Q_pval <- pchisq(Q, Q_df, low = FALSE)
# 	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
# }

mr_egger_fe <- function(b_exp, b_out, se_exp, se_out, isq_threshold = 0.95)
{
	res <- mr_egger_mre(b_exp, b_out, se_exp, se_out, isq_threshold=0)
	res$se <- res$se / res$phi
	res$se_i <- res$se_i / res$phi
	res$pval <- 2 * pt(abs(res$b/res$se), length(b_exp) - 2, lower.tail = FALSE)
	res$pval_i <- 2 * pt(abs(res$b_i/res$se_i), length(b_exp) - 2, lower.tail = FALSE)

	if(res$isq < isq_threshold)
	{
		mod.sim <- simex(
			res$mod,
			B=1000,
			measurement.error = se_exp,
			SIMEXvariable = "beta.exposure",
			fitting.method = "quad",
			asymptotic = "FALSE"
		)
		smod.sim <- summary(mod.sim)
		res$simex_mod <- mod.sim
		res$simex_smod <- smod.sim
		res$simex.b <- smod.sim$coefficients$jackknife[2,1]
		res$simex.b_i <- smod.sim$coefficients$jackknife[1,1]
		res$simex.se <- smod.sim$coefficients$jackknife[2,2] / res$phi
		res$simex.se_i <- smod.sim$coefficients$jackknife[1,2] / res$phi
		res$pval <- 2 * pt(abs(res$simex.b/res$simex.se), length(b_exp) - 2, lower.tail = FALSE)
		res$pval_i <- 2 * pt(abs(res$simex.b_i/res$simex.se_i), length(b_exp) - 2, lower.tail = FALSE)
	}
	return(res)
}

# mr_egger_mre <- function (b_exp, b_out, se_exp, se_out, parameters)
# {
# 	if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3) {
# 		return(list(b = NA, se = NA, pval = NA, nsnp = NA, Q = NA, Q_df = NA, Q_pval = NA))
# 	}

# 	Z <- cbind(1 / se_out, b_exp / se_out)

# 	G <- matrix(b_out / se_out, length(b_out))
# 	p1 <- solve(t(Z) %*% Z)
# 	b <- p1 %*% t(Z) %*% G
# 	se <- sqrt(diag(p1))

# 	wj <- b_exp^2 / se_out^2
# 	bj <- b_out / b_exp

# 	Q <- sum(wj * (bj - b[1]/b_exp - b[2])^2)
# 	phi <- Q / (length(b_exp) - 2)

# 	se <- se * phi
# 	pval <- 2 * pnorm(abs(b/se), low = FALSE)

# 	Q_df <- length(b_exp) - 2
# 	Q_pval <- pchisq(Q, Q_df, low = FALSE)
# 	return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
# }


mr_egger_mre <- function(b_exp, b_out, se_exp, se_out, isq_threshold = 0.95)
{
	stopifnot(length(b_exp) == length(b_out))
	stopifnot(length(se_exp) == length(se_out))
	stopifnot(length(b_exp) == length(se_out))
	nulllist <- list(b = NA, se = NA, pval = NA, nsnp = NA, b_i = NA, se_i = NA, pval_i = NA, Q = NA, Q_df = NA, Q_pval = NA, mod = NA, smod = NA, dat = NA)

	if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3) 
	{
		return(nulllist)
	}

	
	mod <- lm(beta.outcome ~ beta.exposure, data=dat, weights = 1/dat$se.outcome^2, x=TRUE, y=TRUE)

	to_flip <- b_exp < 0
	b_exp[to_flip] <- b_exp[to_flip] * -1
	b_out[to_flip] = b_out[to_flip] * -1 + 2 * mod$coefficients[1]

	dat <- data.frame(beta.exposure = b_exp, beta.outcome = b_out, se.exposure = se_exp, se.outcome = se_out, flipped = to_flip)

	mod <- lm(beta.outcome ~ beta.exposure, data=dat, weights = 1/dat$se.outcome^2, x=TRUE, y=TRUE)
	smod <- summary(mod)

	
	if (nrow(coefficients(smod)) <= 1) 
	{
		warning("Collinearities in MR Egger, try LD pruning the exposure variables.")
		return(nulllist)
	}

	b <- coefficients(smod)[2, 1]
	se <- coefficients(smod)[2, 2]/min(1, smod$sigma)
	pval <- 2 * pt(abs(b/se), length(b_exp) - 2, lower.tail = FALSE)
	b_i <- coefficients(smod)[1, 1]
	se_i <- coefficients(smod)[1, 2]/min(1, smod$sigma)
	pval_i <- 2 * pt(abs(b_i/se_i), length(b_exp) - 2, lower.tail = FALSE)
	# Q <- smod$sigma^2 * (length(b_exp) - 2)
	w <- b_exp^2 / se_out^2
	Q <- sum(w * (b_out / b_exp - b_i / b_exp - b)^2)
	Q_df <- length(b_exp) - 2
	Q_pval <- pchisq(Q, Q_df, low = FALSE)
	phi <- Q / (length(b_exp) - 2)
	isq <- Isq(b_exp/se_out, se_exp/se_out)

	l <- list(b = b, se = se, pval = pval, nsnp = length(b_exp), b_i = b_i, se_i = se_i, pval_i = pval_i, Q = Q, Q_df = Q_df, Q_pval = Q_pval, mod = mod, smod=smod, dat = dat, isq = isq, phi = phi)

	if(isq < isq_threshold)
	{
		mod.sim <- simex(
			mod,
			B=1000,
			measurement.error = se_exp,
			SIMEXvariable = "beta.exposure",
			fitting.method="quad",
			asymptotic = "FALSE"
		)
		smod.sim <- summary(mod.sim)
		l$simex_mod <- mod.sim
		l$simex_smod <- smod.sim
		l$simex.b <- smod.sim$coefficients$jackknife[2,1]
		l$simex.b_i <- smod.sim$coefficients$jackknife[1,1]
		l$simex.se <- smod.sim$coefficients$jackknife[2,2]
		l$simex.se_i <- smod.sim$coefficients$jackknife[1,2]
		l$simex.pval <- smod.sim$coefficients$jackknife[2,4]
		l$simex.pval_i <- smod.sim$coefficients$jackknife[1,4]

	}

	return(l)
}



rucker <- function(dat, Qthresh = 0.05, Igxthresh = 0.9)
{
	A <- mr_ivw_fe(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
	B <- mr_ivw_mre(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
	C <- mr_egger_fe(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, Igxthresh)
	D <- mr_egger_mre(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, Igxthresh)

	mods <- data.frame(
		mod = LETTERS[1:4],
		b = sapply(list(A,B,C,D), function(x) x$b),
		se = sapply(list(A,B,C,D), function(x) x$se),
		pval = sapply(list(A,B,C,D), function(x) x$pval),
		Q = sapply(list(A,B,C,D), function(x) x$Q),
		Q_pval = sapply(list(A,B,C,D), function(x) x$Q_pval)
	)

	Qdiff <- max(0, A$Q - C$Q)
	Qdiff_p <- pchisq(Qdiff, 1, lower.tail=FALSE)


	if(A$Q_pval <= Qthresh)
	{
		if(Qdiff_p <= Qthresh)
		{
			if(C$Q_pval <= Qthresh)
			{
				res <- D
				res$model <- "D"				
			} else {
				res <- C
				res$model <- "C"
			}
		} else {
			res <- B
			res$model <- "B"			
		}
	} else {
		res <- A
		res$model <- "A"
	}

	res$Qdiff <- Qdiff
	res$Qdiff_p <- Qdiff_p
	res$Qr <- C$Q / A$Q

	info <- data.frame(
		key = c("best model", "isq", "egger intercept", "Qdiff", "Qr"),
		value = c(res$model, C$isq, C$b_i, Qdiff, C$Q / A$Q),
		pval = c(NA, NA, C$pval_i, Qdiff_p, NA)
	)

	return(list(res=res, mods=mods, info=info))
}



rucker_ivw_fe <- function(b_exp, b_out, se_exp, se_out)
{
	mod <- summary(lm(b_out / se_out ~ -1 + b_exp/se_out))
	return(mod)
}

