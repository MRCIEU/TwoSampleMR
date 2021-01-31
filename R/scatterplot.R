

#' Create scatter plot with lines showing the causal estimate for different MR tests
#'
#' Requires dev version of ggplot2
#' 
#' @param mr_results Output from \code{\link{mr}}.
#' @param dat Output from \code{\link{harmonise_data}}.
#' @export
#' @return List of plots
mr_scatter_plot <- function(mr_results, dat)
{
	# dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
	requireNamespace("ggplot2", quietly=TRUE)
	requireNamespace("plyr", quietly=TRUE)
	mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
	{
		d <- plyr::mutate(d)
		if(nrow(d) < 2 | sum(d$mr_keep) == 0)
		{
			return(blank_plot("Insufficient number of SNPs"))
		}
		d <- subset(d, mr_keep)
		index <- d$beta.exposure < 0
		d$beta.exposure[index] <- d$beta.exposure[index] * -1
		d$beta.outcome[index] <- d$beta.outcome[index] * -1
		mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
		mrres$a <- 0
		if("MR Egger" %in% mrres$method)
		{
			temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
			mrres$a[mrres$method == "MR Egger"] <- temp$b_i
		}

		if("MR Egger (bootstrap)" %in% mrres$method)
		{
			temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
			mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
		}

		ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
			ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
			ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
			ggplot2::geom_point() +
			ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE) +
			ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
			ggplot2::labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
			ggplot2::theme(legend.position="top", legend.direction="vertical") +
			ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))
	})
	mrres
}


blank_plot <- function(message)
{
	requireNamespace("ggplot2", quietly=TRUE)
	ggplot2::ggplot(data.frame(a=0,b=0,n=message)) + ggplot2::geom_text(ggplot2::aes(x=a,y=b,label=n)) + ggplot2::labs(x=NULL,y=NULL) + ggplot2::theme(axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
}


