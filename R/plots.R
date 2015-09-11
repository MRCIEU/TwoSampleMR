mr_funnel_plot <- function(dat) {}


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
#' <full description>
#'
#' @param leaveoneout_results Output from \code{mr_leaveoneout}
#'
#' @export
#' @return List of plots
mr_leaveoneout_plot <- function(leaveoneout_results)
{
	res <- llply(leaveoneout_results, function(d)
	{
		if(is.null(d))
		{
			return(ggplot(NULL))
		}
		d1 <- subset(d, SNP=="All")
		d2 <- subset(d, SNP!="All")
		ggplot(d2, aes(y=SNP, x=b)) +
		geom_errorbarh(aes(xmin=b-se, xmax=b+se), height=0) +
		geom_point() +
		geom_vline(xintercept=0) +
		geom_vline(xintercept=d1$b, linetype="dashed")
	})
	res
}
