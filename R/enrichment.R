#' Fisher's combined test
#'
#' @param pval Vector of outcome p-values
#' @export
#' @return List with the following elements:
#'         b: MR estimate
#'         se: Standard error
#'         pval: p-value
fishers_combined_test <- function(pval)
{
	pval <- pval[is.finite(pval) & pval <=1 & pval >= 0]
	index <- pval == 0
	if(any(index))
	{
		warning("p-values of 0 are unreliable in Fisher's combined test.")
		pval[index] <- 1e-50
	}
	p <- pchisq(-2 * sum(log(pval)), df=2*length(pval), low=FALSE)
	return(list(pval=p, nsnp=length(pval)))
}





#' Get list of available p-value enrichment methods
#'
#' @export
#' @return Data frame
enrichment_method_list <- function()
{
	a <- list(
		list(
			obj="fishers_combined_test",
			name="Fisher's combined test",
			PubmedID="",
			Description=""
		)
	)
	a <- lapply(a, as.data.frame)
	a <- rbind.fill(a)
	a <- as.data.frame(lapply(a, as.character), stringsAsFactors=FALSE)
	return(a)	
}



enrichment <- function(dat, method_list=enrichment_method_list()$obj)
{
	res <- dlply(subset(dat, mr_keep), .(id.outcome, exposure), function(x)
	{
		message("Performing enrichment analysis of '", x$exposure[1], "' on '", x$displayname.outcome[1], "'")
		res <- lapply(method_list, function(meth)
		{
			get(meth)(x$pval.outcome)
		})
		methl <- enrichment_method_list()
		enrichment_tab <- data.frame(
			Study.ID = x$id.outcome[1],
			Exposure = x$exposure[1],
			Test = methl$name[methl$obj %in% method_list],
			n.SNPs = sapply(res, function(x) x$nsnp),
			pval = sapply(res, function(x) x$pval)
		)
		enrichment_tab <- subset(enrichment_tab, !is.na(pval))
		return(enrichment_tab)
	})
	enrichment_tab <- rbind.fill(lapply(res, function(x) x))

	ao <- available_outcomes()
	ao <- subset(ao, select=c(id, trait, trait_strict, consortium, ethnic, gender, ncase, ncontrol, sample_size, pmid, unit, sd, year, cat, subcat))

	enrichment_tab$ord <- 1:nrow(enrichment_tab)
	enrichment_tab <- merge(enrichment_tab, ao, by.x="Study.ID", by.y="id")
	enrichment_tab <- enrichment_tab[order(enrichment_tab$ord), ]
	enrichment_tab <- subset(enrichment_tab, select=-c(ord))

	return(enrichment_tab)
}
