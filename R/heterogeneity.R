#' Get heterogeneity statistics
#'
#' Get heterogeneity statistics.
#'
#' @param dat Harmonised exposure and outcome data. Output from [harmonise_data()].
#' @param parameters Parameters to be used for various MR methods. Default is output from [default_parameters()].
#' @param method_list List of methods to use in analysis. See [mr_method_list()] for details.
#'
#' @export
#' @return Data frame
mr_heterogeneity <- function(dat, parameters=default_parameters(), method_list = subset(mr_method_list(), heterogeneity_test & use_by_default)$obj) {
	het_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1) {
		# message("Performing MR analysis of '", x$id.exposure[1], "' on '", x$id.outcome[1], "'")
		x <- subset(x1, mr_keep)
		if (nrow(x) < 2) {
			message("Not enough SNPs available for Heterogeneity analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
			return(NULL)
		}
		res <- lapply(method_list, function(meth) {
			get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
		})
		methl <- mr_method_list()
		het_tab <- data.frame(
			outcome = x$outcome[1],
			exposure = x$exposure[1],
			method = methl$name[match(method_list, methl$obj)],
			Q = sapply(res, function(x) x$Q),
			Q_df = sapply(res, function(x) x$Q_df),
			Q_pval = sapply(res, function(x) x$Q_pval)
		)
		het_tab <- subset(het_tab, !(is.na(Q) & is.na(Q_df) & is.na(Q_pval)))
		return(het_tab)
	})

	return(het_tab)
}


#' Test for horizontal pleiotropy in MR analysis
#'
#' Performs MR Egger and returns intercept values.
#'
#' @param dat Harmonised exposure and outcome data. Output from [harmonise_data()].
#'
#' @export
#' @return data frame
mr_pleiotropy_test <- function(dat) {
	ptab <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1) {
		x <- subset(x1, mr_keep)
		if (nrow(x) < 2) {
			message("Not enough SNPs available for pleiotropy analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
			return(NULL)
		}
		res <- mr_egger_regression(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, default_parameters())
		out <- data.frame(
			outcome = x$outcome[1],
			exposure = x$exposure[1],
			egger_intercept = res$b_i,
			se = res$se_i,
			pval = res$pval_i
		)
		return(out)
	})
	return(ptab)
}
