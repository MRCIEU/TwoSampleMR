#' Evaluate the Steiger test's sensitivity to measurement error
#'
#' @param rgx_o Observed variance of exposure explained by SNPs
#' @param rgy_o Observed variance of outcome explained by SNPs
#' @param ... Further arguments to be passed to [lattice::wireframe()]
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{vz}{Total volume of the error parameter space}
#' \item{vz0}{Volume of the parameter space that gives the incorrect answer}
#' \item{vz1}{Volume of the paramtere space that gives the correct answer}
#' \item{sensitivity_ratio}{Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error}
#' \item{pl}{plot of parameter space}
#' }
steiger_sensitivity <- function(rgx_o, rgy_o, ...)
{
	if(rgy_o > rgx_o)
	{
		a <- rgy_o
		b <- rgx_o
	} else {
		a <- rgx_o
		b <- rgy_o
	}

	d <- expand.grid(rxx_o=seq(rgx_o,1,length.out=50), ryy_o=seq(rgy_o,1,length.out=50), type=c("A","B"))
	d$rgy <- rgy_o / d$ryy_o
	d$rgx <- rgx_o / d$rxx_o
	d$z <- d$rgy - d$rgx
	d$z[d$type=="A"] <- 0
	mycolors.trans = grDevices::rgb(c(255,0), c(0,0), 
               c(0,255),alpha = c(70,255), maxColorValue = 255) 

	temp <- lattice::wireframe(
		z ~ rxx_o * ryy_o, 
		groups=type, 
		data=d, 
		scales=list(arrows=FALSE), 
		col.groups = mycolors.trans, 
		drape=FALSE, 
		ylab=expression(rho[xx[o]]), 
		xlab=expression(rho[yy[o]]),
		zlab=expression(rho[gy]-rho[gx]),
		par.settings = list(axis.line=list(col="transparent")),
		...
	)

	vz <- a * log(a) - b * log(b) + a*b*(log(b)-log(a))
	vz0 <- -2*b - b * log(a) - a*b*log(a) + 2*a*b

	vz1 <- abs(vz - vz0)

	sensitivity <- vz0 / (2 * vz0 + abs(vz))
	sensitivity_ratio <- vz1 / vz0

	return(list(
		vz = vz,
		vz0 = vz0,
		vz1 = vz1,
		# sensitivity = sensitivity,
		sensitivity_ratio = sensitivity_ratio,
		pl = temp
	))
}

ss_conf_calcs <- function(bxy, bgx, bux, buy, vg, vu, vex, vey) {
    args <- as.list(environment())
    bxyo <- ((bgx^2*bxy*vg + bux^2*bxy*vu + bux*buy*vu + bxy*vex) / (vg*bgx^2 + vu*bux^2 + vex))
    vx <- bgx^2 * vg + bux^2 * vu + vex
    vy <- bxy^2*bgx^2*vg + (bxy*bux+buy)^2*vu + bxy^2*(vex) + vey
    conf <- bux * vu * buy / (vg * bgx^2 + vu * bux^2 + vex)
    rsqxyo <- bxyo^2 * vx / vy
    rsqxyos <- rsqxyo * sign(bxyo)
    rsqxy <- bxy^2 * vx / vy
    rsqxys <- rsqxy * sign(bxy)
    rsqgx <- bgx^2*vg / (bgx^2 * vg + bux^2 * vu + vex)
    rsqgy <- bgx^2*bxy^2*vg / (bxy^2*bgx^2*vg + (bxy*bux+buy)^2*vu + bxy^2*(vex) + vey)
    rsqux <- bux^2*vu / (bgx^2 * vg + bux^2 * vu + vex)
    rsquy <- (buy + bux * bxy)^2 * vu / (bxy^2*bgx^2*vg + (bxy*bux+buy)^2*vu + bxy^2*(vex) + vey)
    rsquxs <- rsqux * sign(bux)
    rsquys <- rsquy * sign(buy)
    return(c(args, list(
        vx=vx,
        vy=vy,
        bxyo=bxyo, 
        conf=conf, 
        rsqgx=rsqgx, 
        rsqgy=rsqgy, 
        rsqux=rsqux, 
        rsquy=rsquy, 
        rsqxy=rsqxy, 
        rsqxyo=rsqxyo, 
        rsqxyos=rsqxyos, 
        rsquxs=rsquxs, 
        rsquys=rsquys, 
        rsqxys=rsqxys
    )))
}

ss_conf_1d <- function(bxy=0.1, bxyo=0.2, bgx=0.5, vx=1, vy=1, vu=1, vg=0.5, simsize=100) {
    # vx <- bgx^2 * vg + p$bux_vec^2 * vu + vex
    bux_lim <- sqrt((vx - bgx^2 * vg)/vu)
    bux_vec <- seq(-bux_lim, bux_lim, length.out=simsize)
    # Allow causal effect to vary by +/- 200%
    vex <- vx - bgx^2 * vg - bux_vec^2 * vu
    conf <- bxyo - bxy
    buy_vec <- conf * (bgx^2*vg + bux_vec^2*vu + vex) / (bux_vec * vu)
    vey <- vy - (bxy^2*bgx^2*vg + (bxy*bux_vec+buy_vec)^2*vu + bxy^2*vex)
    # vy <- bxy^2*bgx^2*vg + (bxy*bux_vec+buy_vec)^2*vu + bxy^2*vex + vey
    bux_vec * vu * buy_vec / (vg * bgx^2 + vu * bux_vec^2 + vex)
    res <- ss_conf_calcs(bxy, bgx, bux_vec, buy_vec, vg, vu, vex, vey) %>% 
      dplyr::as_tibble() %>%
      dplyr::mutate(bxy=bxy, bgx=bgx, bux=bux_vec, buy=buy_vec, vg=vg, vu=vu, vex=vex, vey=vey)
    return(res)
}

#' Sensitivity analysis for unmeasured confounding on Steiger inferred causal direction
#' 
#' @description 
#' This function takes known parameters from an MR analysis and determines the range of unmeasured 
#' confounding that would be required to agree with the Steiger test inference, and the range of 
#' unmeasured confounding values that would be required to disagree with the Steiger test inference
#' 
#' @param bxy MR estimate of x -> y
#' @param bxyo Observational estimate of x -> y, if not available can re-run with different values to test sensitivity
#' @param bgx SNP-exposure association
#' @param vx Variance of X
#' @param vy Variance of Y
#' @param vg Variance of SNP (approximately 2*p*(1-p) where p is the allele frequency)
#' @param vu Arbitrary variance of unmeasured confounder, default = 1
#' @param simsize Density of search grid, default=10000
#' @param beta_a Weighting of confounder values, which are beta distributed. Specify 'a' parameter of beta distribution, default=1 implying flat prior
#' @param beta_b Weighting of confounder values, which are beta distributed. Specify 'b' parameter of beta distribution, default=1 implying flat prior
#' @param plot Whether to generate plot. Default=TRUE
#' 
#' @importFrom stats dbeta quantile
#' @importFrom rlang .data
#' @export
#' @return List of results
#' \describe{
#' \item{o}{data frame of possible confounding values and agreement with inferred steiger result}
#' \item{prop}{(weighted) fraction of confounding space that agrees with inferred steiger result}
#' \item{pl}{plot}
#' }
steiger_sensitivity_conf <- function(bxy, bxyo, bgx, vx, vy, vg, vu = 1, simsize=10000, beta_a=1, beta_b=1, plot=TRUE)
{
    o <- dplyr::bind_rows(
        ss_conf_1d(bxy=bxy, bxyo=bxyo, bgx=bgx, vx=vx, vy=vy, vu=vu, vg=vg, simsize=simsize) %>%
            dplyr::mutate(direction="inferred"),
        ss_conf_1d(bxy=1/bxy, bxyo=bxyo * vx / vy, bgx=bgx * bxy, vx=vy, vy=vx, vu=vu, vg=vg, simsize=simsize) %>%
            dplyr::mutate(direction="reverse")
    ) %>%
    dplyr::filter(
        .data$vex >= 0 & 
        .data$vey >= 0 &
        .data$rsquy >= 0 & .data$rsquy <= 1 &
        .data$rsqux >= 0 & .data$rsqux <= 1 &
        .data$rsqgx >= 0 & .data$rsqgx <= 1 &
        .data$rsqgy >= 0 & .data$rsqgy <= 1
    ) %>%
        dplyr::group_by(.data$direction) %>%
        dplyr::do({
            x <- .
            x1 <- x$rsqux[-1]
            x2 <- x$rsqux[-length(x$rsqux)]
            y1 <- x$rsquy[-1]
            y2 <- x$rsquy[-length(x$rsquy)]
            d <- sqrt((x1-x2)^2 + (y1-y2)^2)
            d[d > quantile(d, na.rm=T, probs=0.99)*4] <- NA
            x$d <- c(NA, d)
            x$weight <- dbeta(x$rsqux, shape1=beta_a, shape2=beta_b) * dbeta(x$rsquy, shape1=beta_a, shape2=beta_b)
            x
        })
    
    w <- o$d * o$weight
    w1 <- w[o$direction=="inferred"]
    prop <- sum(w1, na.rm=T) / sum(w, na.rm=T)
    ret <- list(result=o, prop=prop)
    if(plot) {
        ret$pl <- ggplot2::ggplot(o, ggplot2::aes(x=.data$rsquxs, y=.data$rsquys)) +
        ggplot2::geom_point(ggplot2::aes(colour=.data$direction, size=.data$weight))
    }
    return(ret)
}

#' MR Steiger test of directionality
#'
#' A statistical test for whether the assumption that exposure causes outcome is valid
#'
#' @param p_exp Vector of p-values of SNP-exposure
#' @param p_out Vector of p-values of SNP-outcome
#' @param n_exp Sample sizes for p_exp
#' @param n_out Sample sizes for p_out
#' @param r_exp Vector of absolute correlations for SNP-exposure
#' @param r_out Vector of absolute correlations for SNP-outcome
#' @param r_xxo Measurememt precision of exposure
#' @param r_yyo Measurement precision of outcome
#' @param ... Further arguments to be passed to [lattice::wireframe()]
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{r2_exp}{Estimated variance explained in x}
#' \item{r2_out}{Estimated variance explained in y}
#' \item{r2_exp_adj}{Predicted variance explained in x accounting for estimated measurement error}
#' \item{r2_out_adj}{Predicted variance explained in y accounting for estimated measurement error}
#' \item{correct_causal_direction}{`TRUE`/`FALSE`}
#' \item{steiger_test}{p-value for inference of direction}
#' \item{correct_causal_direction_adj}{`TRUE`/`FALSE`, direction of causality for given measurement error parameters}
#' \item{steiger_test_adj}{p-value for inference of direction of causality for given measurement error parameters}
#' \item{vz}{Total volume of the error parameter space}
#' \item{vz0}{Volume of the parameter space that gives the incorrect answer}
#' \item{vz1}{Volume of the paramtere space that gives the correct answer}
#' \item{sensitivity_ratio}{Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error}
#' \item{sensitivity_plot}{Plot of parameter space of causal directions and measurement error}
#' }
mr_steiger <- function(p_exp, p_out, n_exp, n_out, r_exp, r_out, r_xxo = 1, r_yyo=1, ...)
{
	r_exp <- abs(r_exp)
	r_out <- abs(r_out)

	ir_exp <- is.na(r_exp)
	ir_out <- is.na(r_out)

	ip_exp <- is.na(p_exp) | is.na(n_exp)
	ip_out <- is.na(p_out) | is.na(n_out)

	if(any(ir_exp))
	{
		r_exp[ir_exp] <- get_r_from_pn(p_exp[ir_exp & !ip_exp], n_exp[ir_exp & !ip_exp])
	}
	if(any(ir_out))
	{
		r_out[ir_out] <- get_r_from_pn(p_out[ir_out & !ip_out], n_out[ir_out & !ip_out])
	}

	r_exp <- sqrt(sum(r_exp[!is.na(r_exp) | is.na(r_out)]^2))
	r_out <- sqrt(sum(r_out[!is.na(r_exp) | is.na(r_out)]^2))

	stopifnot(r_xxo <= 1 & r_xxo >= 0)
	stopifnot(r_yyo <= 1 & r_yyo >= 0)

	r_exp_adj <- sqrt(r_exp^2 / r_xxo^2)
	r_out_adj <- sqrt(r_out^2 / r_yyo^2)

	sensitivity <- steiger_sensitivity(r_exp, r_out, ...)

	rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
	rtest_adj <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp_adj, r34 = r_out_adj)
	l <- list(
		r2_exp = r_exp^2, 
		r2_out = r_out^2, 
		r2_exp_adj = r_exp_adj^2, 
		r2_out_adj = r_out_adj^2, 
		correct_causal_direction = r_exp > r_out, 
		steiger_test = stats::pnorm(-abs(rtest[["z"]])) * 2,
		correct_causal_direction_adj = r_exp_adj > r_out_adj, 
		steiger_test_adj = stats::pnorm(-abs(rtest_adj[["z"]])) * 2,
		vz = sensitivity$vz,
		vz0 = sensitivity$vz0,
		vz1 = sensitivity$vz1,
		sensitivity_ratio = sensitivity$sensitivity_ratio,
		sensitivity_plot = sensitivity$pl
	)
	return(l)
}


#' Perform MR Steiger test of directionality
#'
#' A statistical test for whether the assumption that exposure causes outcome is valid.
#'
#' @param dat Harmonised exposure and outcome data. Output from [harmonise_data()].
#'
#' @export
#' @return List
directionality_test <- function(dat)
{
	if(! all(c("r.exposure", "r.outcome") %in% names(dat)))
	{
		message("r.exposure and/or r.outcome not present.")
		if(! all(c("pval.exposure", "pval.outcome", "samplesize.exposure", "samplesize.outcome") %in% names(dat)))
		{
			message("Can't calculate approximate SNP-exposure and SNP-outcome correlations without pval.exposure, pval.outcome, samplesize.exposure, samplesize.outcome")
			message("Either supply these values, or supply the r.exposure and r.outcome values")
			message("Note, automated correlations assume quantitative traits. For binary traits please pre-calculate in r.exposure and r.outcome e.g. using get_r_from_lor()")
			return(NULL)
		} else {
			message("Calculating approximate SNP-exposure and/or SNP-outcome correlations, assuming all are quantitative traits. Please pre-calculate r.exposure and/or r.outcome using get_r_from_lor() for any binary traits")
		}
	}
	dtest <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x)
	{
		if(!"r.exposure" %in% names(x)) 
		{
			x$r.exposure <- NA
		}
		if(!"r.outcome" %in% names(x)) 
		{
			x$r.outcome <- NA
		}
		b <- mr_steiger(x$pval.exposure, x$pval.outcome, x$samplesize.exposure, x$samplesize.outcome, x$r.exposure, x$r.outcome)
		a <- data.frame(
			exposure = x$exposure[1],
			outcome = x$outcome[1],
			snp_r2.exposure = b$r2_exp,
			snp_r2.outcome = b$r2_out,
			correct_causal_direction = b$correct_causal_direction,
			steiger_pval = b$steiger_test
		)
		return(a)
	})

	return(dtest)
}



#' MR Steiger test of directionality
#'
#' A statistical test for whether the assumption that exposure causes outcome is valid
#'
#' @param r_exp Vector of correlations of SNP-exposure
#' @param r_out Vector of correlations of SNP-outcome
#' @param n_exp Sample sizes for p_exp
#' @param n_out Sample sizes for p_out
#' @param r_xxo Measurememt precision of exposure
#' @param r_yyo Measurement precision of outcome
#' @param ... Further arguments to be passed to [lattice::wireframe()]
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{r2_exp}{Estimated variance explained in x}
#' \item{r2_out}{Estimated variance explained in y}
#' \item{r2_exp_adj}{Predicted variance explained in x accounting for estimated measurement error}
#' \item{r2_out_adj}{Predicted variance explained in y accounting for estimated measurement error}
#' \item{correct_causal_direction}{`TRUE`/`FALSE`}
#' \item{steiger_test}{p-value for inference of direction}
#' \item{correct_causal_direction_adj}{`TRUE`/`FALSE`, direction of causality for given measurement error parameters}
#' \item{steiger_test_adj}{p-value for inference of direction of causality for given measurement error parameters}
#' \item{vz}{Total volume of the error parameter space}
#' \item{vz0}{Volume of the parameter space that gives the incorrect answer}
#' \item{vz1}{Volume of the paramtere space that gives the correct answer}
#' \item{sensitivity_ratio}{Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error}
#' \item{sensitivity_plot}{Plot of parameter space of causal directions and measurement error}
#' }
mr_steiger2 <- function(r_exp, r_out, n_exp, n_out, r_xxo = 1, r_yyo=1, ...)
{
	index <- any(is.na(r_exp)) | any(is.na(r_out)) | any(is.na(n_exp)) | any(is.na(n_out))
	n_exp <- n_exp[!index]
	n_out <- n_out[!index]

	r_exp <- sqrt(sum(r_exp^2))
	r_out <- sqrt(sum(r_out^2))

	stopifnot(r_xxo <= 1 & r_xxo >= 0)
	stopifnot(r_yyo <= 1 & r_yyo >= 0)

	r_exp_adj <- sqrt(r_exp^2 / r_xxo^2)
	r_out_adj <- sqrt(r_out^2 / r_yyo^2)

	sensitivity <- steiger_sensitivity(r_exp, r_out, ...)

	rtest <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
	rtest_adj <- psych::r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp_adj, r34 = r_out_adj)
	l <- list(
		r2_exp = r_exp^2, 
		r2_out = r_out^2, 
		r2_exp_adj = r_exp_adj^2, 
		r2_out_adj = r_out_adj^2, 
		correct_causal_direction = r_exp > r_out, 
		steiger_test = stats::pnorm(-abs(rtest[["z"]]))*2,
		correct_causal_direction_adj = r_exp_adj > r_out_adj, 
		steiger_test_adj = stats::pnorm(-abs(rtest_adj[["z"]]))*2,
		vz = sensitivity$vz,
		vz0 = sensitivity$vz0,
		vz1 = sensitivity$vz1,
		sensitivity_ratio = sensitivity$sensitivity_ratio,
		sensitivity_plot = sensitivity$pl
	)
	return(l)
}
