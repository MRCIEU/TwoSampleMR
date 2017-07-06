
#' Wrapper for MR-PRESSO
#'
#' See https://github.com/rondolab/MR-PRESSO
#'
#' @param dat Output from harmonise_data
#' @param NbDistribution = 1000 Number of bootstraps
#' @param SignifThreshold = 0.05 Outlier significance threshold
#'
#' @export
#' @return List of results for every exposure/outcome combination
run_mr_presso <- function(dat, NbDistribution = 1000,  SignifThreshold = 0.05)
{
	require(MRPRESSO)
	dat <- subset(dat, mr_keep)
	d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)), select=c(exposure, outcome, id.exposure, id.outcome))
	res <- list()
	attributes(res)$id.exposure <- d$id.exposure
	attributes(res)$id.outcome <- d$id.outcome
	attributes(res)$exposure <- d$exposure
	attributes(res)$outcome <- d$id.exposure
	for(j in 1:nrow(d))
	{
		x <- subset(dat, exposure == d$exposure[j] & outcome == d$outcome[j])
		message(x$exposure[1], " - ", x$outcome[1])
		res[[j]] <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = x, NbDistribution = NbDistribution,  SignifThreshold = SignifThreshold)
	}
	return(res)
}
