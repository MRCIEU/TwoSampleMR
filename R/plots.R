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
			geom_point() +
			geom_abline(data=mrres, aes(intercept=a, slope=b, colour=Test), show.legend=TRUE) +
			scale_colour_brewer(type="qual") +
			labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
			theme(legend.position="top", legend.direction="vertical") +
			guides(colour=guide_legend(ncol=2))
	})
	mrres
}
