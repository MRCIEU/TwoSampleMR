#' Convenient function to create a harmonised dataset
#'
#' <full description>
#' @param exposures default 2 and 301 (bmi and ldl)
#' @param outcomes default 7 and 1001 (edu and chd)
#'
#' @export
#' @return Harmonised data frame
make_dat <- function(exposures=c(2,301), outcomes=c(7,1001))
{
	a <- extract_instruments(exposures)
	b <- extract_outcome_data(a$SNP, outcomes)
	return(harmonise_data(a,b))
}
