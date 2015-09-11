mr_forest_plot <- function() {}


#' Create scatter plot with lines showing the causal estimate for different MR tests
#'
#' Requires dev version of ggplot2
#' 
#' @param mr_results Output from \code{mr}
#' @param dat Output from \code{harmonise_exposure_outcome}
#' @export
#' @return List of plots
mr_scatter_plot <- function(mr_results, dat)
{
	mrres <- dlply(dat, .(exposure, outcome), function(d)
	{
		mrres <- subset(mr_results$mr, Exposure == d$exposure[1] & Outcome == d$outcome[1])
		egger <- subset(mr_results$extra, Exposure == d$exposure[1] & Outcome == d$outcome[1])

		mrres$a <- 0
		mrres$a[mrres$Test == "Egger regression"] <- egger$b[1]

		ggplot(data=d, aes(x=beta.exposure, y=beta.outcome)) +
			geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
			geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", width=0) +

			geom_point() +
			geom_abline(data=mrres, aes(intercept=a, slope=b, colour=Test), show.legend=TRUE) +
			scale_colour_brewer(type="qual") +
			labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
			theme(legend.position="top", legend.direction="vertical") +
			guides(colour=guide_legend(ncol=2))
	})
	mrres
}


#' Plot results from leaveoneout analysis
#'
#' @param leaveoneout_results Output from \code{mr_leaveoneout}
#'
#' @export
#' @return List of plots
mr_leaveoneout_plot <- function(leaveoneout_results)
{
	res <- dlply(leaveoneout_results, .(exposure, outcome), function(d)
	{
		if(nrow(d) == 1)
		{
			return(NULL)
		}
		d1 <- subset(d, SNP=="All")
		d2 <- subset(d, SNP!="All")
		d2$up <- d2$b + 1.96 * d2$se
		d2$lo <- d2$b - 1.96 * d2$se
		ggplot(d2, aes(y=SNP, x=b)) +
		geom_errorbarh(aes(xmin=lo, xmax=up), height=0) +
		geom_point() +
		geom_vline(xintercept=0) +
		geom_vline(xintercept=d1$b, linetype="dashed")
	})
	res
}

#' Forest plot
#'
#' @param singlesnp_results from \code{mr_singlesnp}
#'
#' @export
#' @return List of plots
mr_forest_plot <- function(singlesnp_results)
{
	res <- dlply(singlesnp_results, .(exposure, outcome), function(d)
	{
		if(nrow(d) < 3) return(NULL)
		d$up <- d$b + 1.96 * d$se
		d$lo <- d$b - 1.96 * d$se
		d$tot <- 3
		d$tot[d$SNP != "All"] <- 0.1
		d$SNP <- as.character(d$SNP)
		nom <- d$SNP[d$SNP != "All"]
		d <- rbind(d, d[nrow(d),])
		d$SNP[nrow(d)-1] <- ""
		d$b[nrow(d)-1] <- NA
		d$up[nrow(d)-1] <- NA
		d$lo[nrow(d)-1] <- NA
		d$SNP <- ordered(d$SNP, levels=c("All", "", nom))

		ggplot(d, aes(y=SNP, x=b)) +
		geom_vline(xintercept=0, linetype="dotted") +
		geom_errorbarh(aes(xmin=pmax(lo, min(d$b, na.rm=T)), xmax=pmin(up, max(d$b, na.rm=T)), size=tot, colour=as.factor(tot)), height=0) +
		geom_point(aes(colour=as.factor(tot))) +
		geom_hline(aes(yintercept = which(levels(SNP) %in% ""))) +
		scale_colour_manual(values=c("black", "red")) +
		xlim(c(min(c(0, d$b), na.rm=T), max(c(0, d$b), na.rm=T))) +
		theme(legend.position="none")

	})
	res
}





