

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
	dat <- subset(dat, paste(id.outcome, exposure) %in% paste(mr_results$Study.ID, mr_results$Exposure))
	mrres <- dlply(dat, .(exposure, id.outcome), function(d)
	{
		if(nrow(d) < 2)
		{
			return(blank_plot("Insufficient number of SNPs"))
		}
		mrres <- subset(mr_results, Exposure == d$exposure[1] & Study.ID == d$id.outcome[1])
		# egger <- subset(mr_results$extra, Exposure == d$exposure[1] & Outcome == d$outcome[1])

		mrres$a <- 0
		# mrres$a[mrres$Test == "Egger regression"] <- egger$b[1]

		ggplot(data=d, aes(x=beta.exposure, y=beta.outcome)) +
			geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
			geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +

			geom_point() +
			geom_abline(data=mrres, aes(intercept=a, slope=b, colour=Test), show.legend=TRUE) +
			# scale_colour_brewer(type="qual") +
			scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
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
mr_leaveoneout_plot <- function(singlesnp_results)
{
	res <- dlply(singlesnp_results, .(exposure, id.outcome), function(d)
	{
		if(nrow(d) < 3) {
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

		ggplot(d, aes(y=SNP, x=b)) +
		geom_vline(xintercept=0, linetype="dotted") +
		# geom_errorbarh(aes(xmin=pmax(lo, min(d$b, na.rm=T)), xmax=pmin(up, max(d$b, na.rm=T)), size=as.factor(tot), colour=as.factor(tot)), height=0) +
		geom_errorbarh(aes(xmin=lo, xmax=up, size=as.factor(tot), colour=as.factor(tot)), height=0) +
		geom_point(aes(colour=as.factor(tot))) +
		geom_hline(aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
		scale_colour_manual(values=c("black", "red")) +
		scale_size_manual(values=c(0.3, 1)) +
		# xlim(c(min(c(0, d$b), na.rm=T), max(c(0, d$b), na.rm=T))) +
		theme(
			legend.position="none", 
			axis.text.y=element_text(size=8), 
			axis.ticks.y=element_line(size=0),
			axis.title.x=element_text(size=8)) +
		labs(y="", x=paste0("MR leave-one-out sensitivity analysis for\n'", d$exposure[1], "' on '", d$outcome[1], "'"))
	})
	res
}


blank_plot <- function(message)
{
	ggplot(data.frame(a=0,b=0,n=message)) + geom_text(aes(x=a,y=b,label=n)) + labs(x=NULL,y=NULL) + theme(axis.text=element_blank(), axis.ticks=element_blank())
}


#' Forest plot
#'
#' @param singlesnp_results from \code{mr_singlesnp}
#'
#' @export
#' @return List of plots
mr_forest_plot <- function(singlesnp_results)
{
	res <- dlply(singlesnp_results, .(exposure, id.outcome), function(d)
	{
		if(nrow(d) < 3) {
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

		ggplot(d, aes(y=SNP, x=b)) +
		geom_vline(xintercept=0, linetype="dotted") +
		# geom_errorbarh(aes(xmin=pmax(lo, min(d$b, na.rm=T)), xmax=pmin(up, max(d$b, na.rm=T)), size=as.factor(tot), colour=as.factor(tot)), height=0) +
		geom_errorbarh(aes(xmin=lo, xmax=up, size=as.factor(tot), colour=as.factor(tot)), height=0) +
		geom_point(aes(colour=as.factor(tot))) +
		geom_hline(aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
		scale_colour_manual(values=c("black", "red")) +
		scale_size_manual(values=c(0.3, 1)) +
		# xlim(c(min(c(0, d$b), na.rm=T), max(c(0, d$b), na.rm=T))) +
		theme(
			legend.position="none", 
			axis.text.y=element_text(size=8), 
			axis.ticks.y=element_line(size=0),
			axis.title.x=element_text(size=8)) +
		labs(y="", x=paste0("MR effect size for\n'", d$exposure[1], "' on '", d$outcome[1], "'"))
	})
	res
}





