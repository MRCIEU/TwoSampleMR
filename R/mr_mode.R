#' MR mode estimators
#'
#' Perform simple, weighted, penalised modes, as well as versions that use the NOME assumption.
#'
#' @md
#' @param dat Output from [`harmonise_data()`].
#' @param parameters List of parameters. The default is `default_parameters()`.
#' @param mode_method The default is `"all"`. The other choices are `'Simple mode'`, `'Weighted mode'`, `'Penalised mode'`, `'Simple mode (NOME)'`, `'Weighted mode (NOME)'`.
#'
#' @export
#' @return data frame
#' @importFrom stats pchisq pt qnorm
mr_mode <- function(dat, parameters=default_parameters(), mode_method="all")
{
	if("mr_keep" %in% names(dat)) dat <- subset(dat, mr_keep)

	if(nrow(dat) < 3) 
	{
		warning("Need at least 3 SNPs")
		return(NULL)
	}

	b_exp <- dat$beta.exposure
	b_out <- dat$beta.outcome
	se_exp <- dat$se.exposure
	se_out <- dat$se.outcome

	#--------------------------------------#
	#Function to compute the point estimate#
	#--------------------------------------#
	#BetaIV.in: ratio estimates
	#seBetaIV.in: standard errors of ratio estimates
	#' @importFrom stats mad sd
	beta <- function(BetaIV.in, seBetaIV.in, phi)
	{
		#Bandwidth rule - modified Silverman's rule proposed by Bickel (2002)
		s <- 0.9*(min(sd(BetaIV.in), mad(BetaIV.in)))/length(BetaIV.in)^(1/5)

		#Standardised weights
		weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)

		beta <- NULL

		for(cur_phi in phi)
		{
			#Define the actual bandwidth
			h <- max(0.00000001, s*cur_phi)
			#Compute the smoothed empirical density function
			densityIV <- density(BetaIV.in, weights=weights, bw=h)
			#Extract the point with the highest density as the point estimate 
			beta[length(beta)+1] <- densityIV$x[densityIV$y==max(densityIV$y)]
		}
		return(beta)
	}

	#------------------------------------------#
	#Function to estimate SEs through bootstrap#
	#------------------------------------------#
	#BetaIV.in: ratio estimates
	#seBetaIV.in: standard errors of ratio estimates
	#beta_Mode.in: point causal effect estimates
	#' @importFrom stats density pchisq rnorm
	boot <- function(BetaIV.in, seBetaIV.in, beta_Mode.in, nboot)
	{
		#Set up a matrix to store the results from each bootstrap iteration
		beta.boot <- matrix(nrow=nboot, ncol=length(beta_Mode.in))

		for(i in 1:nboot) 
		{
			#Re-sample each ratio estimate using SEs derived not assuming NOME
			BetaIV.boot      <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,1])
			#Re-sample each ratio estimate using SEs derived under NOME
			BetaIV.boot_NOME <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,2])

			#Simple mode, not assuming NOME
			beta.boot[i,1:length(phi)] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
			#Weighted mode, not assuming NOME
			beta.boot[i,(length(phi)+1):(2*length(phi))] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=seBetaIV.in[,1], phi=phi)
			#Penalised mode, not assuming NOME
			weights <- 1/seBetaIV.in[,1]^2
			penalty <- pchisq(weights * (BetaIV.boot-beta.boot[i,(length(phi)+1):(2*length(phi))])^2, df=1, lower.tail=FALSE)
			pen.weights <- weights*pmin(1, penalty*parameters$penk)

			beta.boot[i,(2*length(phi)+1):(3*length(phi))] <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=sqrt(1/pen.weights), phi=phi)

			#Simple mode, assuming NOME
			beta.boot[i,(3*length(phi)+1):(4*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)
			#Weighted mode, assuming NOME
			beta.boot[i,(4*length(phi)+1):(5*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=seBetaIV.in[,2], phi=phi)
		}
		return(beta.boot)
	}

	# Parameters
	phi <- parameters$phi
	nboot <- parameters$nboot
	alpha <- parameters$alpha

	#Ratio estimates
	BetaIV   <- b_out/b_exp

	#SEs of ratio estimates
	seBetaIV <- cbind(sqrt((se_out^2)/(b_exp^2) + ((b_out^2)*(se_exp^2))/(b_exp^4)), #SEs NOT assuming NOME
	se_out/abs(b_exp)) #SEs ASSUMING NOME

	#Point causal effect estimate using the simple mode
	beta_SimpleMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=rep(1, length(BetaIV)), phi=phi)

	#Point causal effect estimate using the weighted mode (not asusming NOME)
	beta_WeightedMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,1], phi=phi)
	weights <- 1/seBetaIV[,1]^2
	penalty <- pchisq(weights * (BetaIV-beta_WeightedMode)^2, df=1, lower.tail=FALSE)
	pen.weights <- weights*pmin(1, penalty*parameters$penk) # penalized 

	beta_PenalisedMode <- beta(BetaIV.in=BetaIV, seBetaIV.in=sqrt(1/pen.weights), phi=phi)

	#Point causal effect estimate using the weighted mode (asusming NOME)
	beta_WeightedMode_NOME <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,2], phi=phi)

	#Combine all point effect estimates in a single vector
	beta_Mode <- rep(c(beta_SimpleMode, beta_WeightedMode, beta_PenalisedMode, beta_SimpleMode, beta_WeightedMode_NOME))

	#Compute SEs, confidence intervals and P-value
	beta_Mode.boot <- boot(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV, beta_Mode.in=beta_Mode, nboot=nboot)
	se_Mode <- apply(beta_Mode.boot, 2, mad)

	CIlow_Mode <- beta_Mode-qnorm(1-alpha/2)*se_Mode
	CIupp_Mode <- beta_Mode+qnorm(1-alpha/2)*se_Mode

	P_Mode <- pt(abs(beta_Mode/se_Mode), df=length(b_exp)-1, lower.tail=FALSE)*2

	#Vector to indicate the method referring to each row
	Method <- rep(c('Simple mode', 'Weighted mode', 'Penalised mode', 'Simple mode (NOME)', 'Weighted mode (NOME)'), each=length(phi))

	#Return a data frame containing the results
	id.exposure <- ifelse("id.exposure" %in% names(dat), dat$id.exposure[1], "")
	id.outcome <- ifelse("id.outcome" %in% names(dat), dat$id.outcome[1], "")
	Results <- data.frame(
		id.exposure = id.exposure, 
		id.outcome = id.outcome, 
		method = Method, 
		nsnp = length(b_exp), 
		b = beta_Mode, 
		se = se_Mode, 
		ci_low = CIlow_Mode, 
		ci_upp = CIupp_Mode, 
		pval = P_Mode, 
		stringsAsFactors=FALSE
	)

	if(mode_method == "all")
	{
		return(Results)
	} else {
		stopifnot(all(mode_method %in% Results$method))
		i <- which(Results$method == mode_method)
		return(list(b = Results$b[i], se = Results$se[i], pval=Results$pval[i], nsnp=length(b_exp)))
	}

	return(Results)
}




#' MR weighted mode estimator
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param parameters List containing `phi` - Bandwidth parameter, and `nboot` - number of bootstraps to calculate SE. `default_parameters()` sets `list(phi=1, nboot=1000)`.
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{MR estimate}
#' \item{se}{Standard error}
#' \item{pval}{p-value}
#' }
mr_weighted_mode <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters()) 
{
	index <- !is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)
	if(sum(index) < 3)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	b_exp <- b_exp[index]
	b_out <- b_out[index]
	se_exp <- se_exp[index]
	se_out <- se_out[index]

	return(mr_mode(data.frame(beta.exposure=b_exp, beta.outcome=b_out, se.exposure=se_exp, se.outcome=se_out), parameters=parameters, mode_method="Weighted mode"))
}


#' MR simple mode estimator
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param parameters List containing `phi` - Bandwidth parameter, and `nboot` - number of bootstraps to calculate SE. `default_parameters()` sets `list(phi=1, nboot=1000)`.
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{MR estimate}
#' \item{se}{Standard error}
#' \item{pval}{p-value}
#' }
mr_simple_mode <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters()) 
{
	index <- !is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)
	if(sum(index) < 3)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	b_exp <- b_exp[index]
	b_out <- b_out[index]
	se_exp <- se_exp[index]
	se_out <- se_out[index]

	return(mr_mode(data.frame(beta.exposure=b_exp, beta.outcome=b_out, se.exposure=se_exp, se.outcome=se_out), parameters=parameters, mode_method="Simple mode"))
}


#' MR weighted mode estimator (NOME)
#'
#' Weighted mode estimator
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param parameters List containing `phi` - Bandwidth parameter, and `nboot` - number of bootstraps to calculate SE. `default_parameters()` sets `list(phi=1, nboot=1000)`.
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{MR estimate}
#' \item{se}{Standard error}
#' \item{pval}{p-value}
#' }
mr_weighted_mode_nome <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters()) 
{
	index <- !is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)
	if(sum(index) < 3)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	b_exp <- b_exp[index]
	b_out <- b_out[index]
	se_exp <- se_exp[index]
	se_out <- se_out[index]

	return(mr_mode(data.frame(beta.exposure=b_exp, beta.outcome=b_out, se.exposure=se_exp, se.outcome=se_out), parameters=parameters, mode_method="Weighted mode (NOME)"))
}


#' MR simple mode estimator (NOME)
#' 
#' MR simple mode estimator (NOME).
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param parameters List containing `phi` - Bandwidth parameter, and `nboot` - number of bootstraps to calculate SE. `default_parameters()` sets `list(phi=1, nboot=1000)`.
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{MR estimate}
#' \item{se}{Standard error}
#' \item{pval}{p-value}
#' }
mr_simple_mode_nome <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters()) 
{
	index <- !is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)
	if(sum(index) < 3)
	return(list(b=NA, se=NA, pval=NA, nsnp=NA))

	b_exp <- b_exp[index]
	b_out <- b_out[index]
	se_exp <- se_exp[index]
	se_out <- se_out[index]

	return(mr_mode(data.frame(beta.exposure=b_exp, beta.outcome=b_out, se.exposure=se_exp, se.outcome=se_out), parameters=parameters, mode_method="Simple mode (NOME)"))
}


mr_mode_broken <- function(dat, parameters=default_parameters(), mode_method="all")
{
	res <- plyr::ddply(dat, c("id.exposure", "exposure", "id.outcome", "outcome"), function(x) mr_mode_internal(x, parameters, mode_method="all"))
	if(mode_method != "all")
	{
		return(subset(res, method %in% mode_method))
	} else {
		return(res)
	}
}
