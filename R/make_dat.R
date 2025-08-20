#' Convenient function to create a harmonised dataset
#'
#' Convenient function to create a harmonised dataset.
#'
#' @param exposures The default is `c("ieu-a-2", "ieu-a-301")` (BMI and LDL).
#' @param outcomes The default is `c("ieu-a-7", "ieu-a-1001")` (CHD and EDU).
#' @param proxies Look for proxies? Default = `TRUE`
#'
#' @export
#' @return Harmonised data frame
make_dat <- function(
  exposures = c("ieu-a-2", "ieu-a-301"),
  outcomes = c("ieu-a-7", "ieu-a-1001"),
  proxies = TRUE
) {
  a <- extract_instruments(exposures)
  b <- extract_outcome_data(a$SNP, outcomes, proxies = proxies)
  return(harmonise_data(a, b, action = 1))
}
