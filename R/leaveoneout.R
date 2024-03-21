#' Leave one out sensitivity analysis
#'
#' @param dat Output from [harmonise_data()].
#' @param method Choose which method to use. The default is `mr_ivw`.
#' @param parameters List of parameters.
#'
#' @export
#' @return List of data frames
mr_leaveoneout <- function(dat, parameters=default_parameters(), method=mr_ivw)
{
	if(!"samplesize.outcome" %in% names(dat))
	{
		dat$samplesize.outcome <- NA
	}

	stopifnot("outcome" %in% names(dat))
	stopifnot("exposure" %in% names(dat))
	stopifnot("beta.exposure" %in% names(dat))
	stopifnot("beta.outcome" %in% names(dat))
	stopifnot("se.exposure" %in% names(dat))
	stopifnot("se.outcome" %in% names(dat))


	res <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(X)
	{
		x <- subset(X, mr_keep)
		nsnp <- nrow(x)
		if(nsnp == 0)
		{
			x <- X[1,]
			d <- data.frame(
				SNP = "All",
				b = NA,
				se = NA,
				p = NA,
				samplesize = NA,
				outcome = x$outcome[1],
				exposure = x$exposure[1]
			)
			return(d)
		}
		if(nsnp > 2)
		{
			l <- lapply(1:nsnp, function(i)
			{
				with(x, method(beta.exposure[-i], beta.outcome[-i], se.exposure[-i], se.outcome[-i], parameters))
			})
			l[[nsnp+1]] <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
			d <- data.frame(
				SNP = c(as.character(x$SNP), "All"),
				b = sapply(l, function(y) y$b),
				se = sapply(l, function(y) y$se),
				p = sapply(l, function(y) y$pval),
				samplesize = x$samplesize.outcome[1]
			)
			d$outcome <- x$outcome[1]
			d$exposure <- x$exposure[1]

		} else {
			a <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
			d <- data.frame(
				SNP = "All",
				b = a$b,
				se = a$se,
				p = a$pval,
				samplesize = x$samplesize.outcome[1]
			)
			d$outcome <- x$outcome[1]
			d$exposure <- x$exposure[1]
		}
		return(d)
	})
	res <- subset(res, select=c(exposure, outcome, id.exposure, id.outcome, samplesize, SNP, b, se, p))
	return(res)
}



#' Plot results from leaveoneout analysis
#' 
#' Plot results from leaveoneout analysis.
#'
#' @param leaveoneout_results Output from [mr_leaveoneout()].
#'
#' @export
#' @return List of plots
mr_leaveoneout_plot <- function(leaveoneout_results)
{
	res <- plyr::dlply(leaveoneout_results, c("id.exposure", "id.outcome"), function(d)
	{
		d <- plyr::mutate(d)
		# Need to have at least 3 SNPs because IVW etc methods can't be performed with fewer than 2 SNPs
		if(sum(!grepl("All", d$SNP)) < 3) {
			return(
				blank_plot("Insufficient number of SNPs")
			)
		}
		d$up <- d$b + 1.96 * d$se
		d$lo <- d$b - 1.96 * d$se
		d$tot <- 1
		d$tot[d$SNP != "All"] <- 0.01
		d$SNP <- as.character(d$SNP)
		nom <- d$SNP[d$SNP != "All"]
		nom <- nom[order(d$b)]
		d <- rbind(d, d[nrow(d),])
		d$SNP[nrow(d)-1] <- ""
		d$b[nrow(d)-1] <- NA
		d$up[nrow(d)-1] <- NA
		d$lo[nrow(d)-1] <- NA
		d$SNP <- ordered(d$SNP, levels=c("All", "", nom))

		ggplot2::ggplot(d, ggplot2::aes(y=SNP, x=b)) +
		ggplot2::geom_vline(xintercept=0, linetype="dotted") +
		# ggplot2::geom_errorbarh(ggplot2::aes(xmin=pmax(lo, min(d$b, na.rm=T)), xmax=pmin(up, max(d$b, na.rm=T)), size=as.factor(tot), colour=as.factor(tot)), height=0) +
		ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo, xmax=up, linewidth=as.factor(tot), colour=as.factor(tot)), height=0) +
		ggplot2::geom_point(ggplot2::aes(colour=as.factor(tot))) +
		ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
		ggplot2::scale_colour_manual(values=c("black", "red")) +
		ggplot2::scale_linewidth_manual(values=c(0.3, 1)) +
		# xlim(c(min(c(0, d$b), na.rm=T), max(c(0, d$b), na.rm=T))) +
		ggplot2::theme(
			legend.position="none", 
			axis.text.y=ggplot2::element_text(size=8), 
			axis.ticks.y=ggplot2::element_line(linewidth=0),
			axis.title.x=ggplot2::element_text(size=8)) +
		ggplot2::labs(y="", x=paste0("MR leave-one-out sensitivity analysis for\n'", d$exposure[1], "' on '", d$outcome[1], "'"))
	})
	res
}
