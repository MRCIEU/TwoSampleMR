#' Fisher's combined test
#'
#' @param pval Vector of outcome p-values
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{MR estimate}
#' \item{se}{Standard error}
#' \item{pval}{p-value}
#' }
fishers_combined_test <- function(pval) {
	pval <- pval[is.finite(pval) & pval <=1 & pval >= 0]
	index <- pval == 0
	if (any(index)) {
		warning("p-values of 0 are unreliable in Fisher's combined test.")
		pval[index] <- 1e-50
	}
	p <- stats::pchisq(-2 * sum(log(pval)), df=2*length(pval), lower.tail=FALSE)
	return(list(pval=p, nsnp=length(pval)))
}





#' Get list of available p-value enrichment methods
#'
#' @export
#' @return Data frame
enrichment_method_list <- function() {
	a <- list(
		list(
			obj="fishers_combined_test",
			name="Fisher's combined test",
			PubmedID="",
			Description=""
		)
	)
	a <- lapply(a, as.data.frame)
	a <- plyr::rbind.fill(a)
	a <- as.data.frame(lapply(a, as.character), stringsAsFactors=FALSE)
	return(a)
}


#' Perform enrichment analysis
#'
#'
#' @param dat Harmonised exposure and outcome data. Output from [harmonise_data()].
#' @param method_list List of methods to use in analysis. Default is \code{enrichment_method_list()$obj}. See [enrichment_method_list()] for details.
#'
#' @export
#' @return data frame
enrichment <- function(dat, method_list=enrichment_method_list()$obj) {
	res <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1) {
		# message("Performing enrichment analysis of '", x$id.exposure[1], "' on '", x$id.outcome[1], "'")

		x <- subset(x1, !is.na(pval.outcome))
		if (nrow(x) == 0) {
			message("No outcome p-values for analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
			return(NULL)
		}
		res <- lapply(method_list, function(meth) {
			get(meth)(x$pval.outcome)
		})
		methl <- enrichment_method_list()
		enrichment_tab <- data.frame(
			outcome = x$outcome[1],
			exposure = x$exposure[1],
			method = methl$name[methl$obj %in% method_list],
			nsnp = sapply(res, function(x) x$nsnp),
			pval = sapply(res, function(x) x$pval)
		)
		enrichment_tab <- subset(enrichment_tab, !is.na(pval))
		return(enrichment_tab)
	})
	return(res)
}
