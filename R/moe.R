#' Mixture of experts
#'
#' Based on the method described here https://www.biorxiv.org/content/early/2017/08/23/173682
#' Once all MR methods have been applied to a summary set, you can then use the mixture of experts to predict the method most likely to be the most accurate.
#'
#' @param dat Output from \code{mr_wrapper}. 
#' @param rf The trained random forest for the methods. This is available to download at https://www.dropbox.com/s/5la7y38od95swcf/rf.rdata?dl=0
#'
#' @export
#' @return List
#' @examples
#' 
#' # Example of body mass index on coronary heart disease
#' # Extract and harmonise data
#' a <- extract_instruments(2)
#' b <- extract_outcome_data(a$SNP, 7)
#' ab <- harmonise_data(a,b)
#' 
#' # Apply all MR methods
#' r <- mr_wrapper(ab)
#' 
#' # Update the results with mixture of experts
#' r <- mr_moe(r)
#' 
#' # Now you can view the estimates, and see that they have 
#' # been sorted in order from most likely to least likely to 
#' # be accurate, based on MOE prediction
#' r[[1]]$estimates
mr_moe <- function(res, rf)
{
	require(dplyr)
	require(randomForest)
	lapply(res, function(x)
	{
		try(mr_moe_single(x, rf))
	})
}


mr_moe_single <- function(res, rf)
{
	require(dplyr)
	require(randomForest)
	metric <- res$info[1,] %>% dplyr::select(-c(id.exposure, id.outcome, steiger_filtered, outlier_filtered, nsnp_removed))

	methodlist <- names(rf)
	pred <- lapply(methodlist, function(m)
	{
		d <- tibble(
			method = m,
			MOE = predict(rf[[m]], metric, type="prob")[,2]
		)
		return(d)
	}) %>% bind_rows %>% arrange(desc(MOE))
	res$estimates$selection <- "DF + HF"
	res$estimates$selection[!res$estimates$outlier_filtered & res$estimates$steiger_filtered] <- "DF"
	res$estimates$selection[res$estimates$outlier_filtered & !res$estimates$steiger_filtered] <- "HF"
	res$estimates$selection[!res$estimates$outlier_filtered & !res$estimates$steiger_filtered] <- "Tophits"
	res$estimates$method <- paste(res$estimates$method, "-", res$estimates$selection)
	res$estimates <- inner_join(res$estimates, pred, by="method") %>% arrange(desc(MOE))
	return(res)
}
