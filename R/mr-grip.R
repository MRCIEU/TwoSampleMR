#' MR-GRIP: Allele coding invariant version of MR-Egger
#'
#' @param b_exp Vector of genetic effects on exposure.
#' @param b_out Vector of genetic effects on outcome.
#' @param se_exp Standard errors of genetic effects on exposure.
#' @param se_out Standard errors of genetic effects on outcome.
#' @param parameters List of parameters.
#'
#' @export
#' @return List of with the following elements:
#' \describe{
#' \item{b}{MR estimate}
#' \item{se}{Standard error of MR estimate}
#' \item{pval}{p-value of MR estimate}
#' \item{Q, Q_df, Q_pval}{Heterogeneity stats}
#' \item{b.wi}{MR estimate adjusting for weak instruments}
#' \item{se.wi}{Standard error adjusting for weak instruments}
#' \item{pval.wi}{p-value adjusting for weak instruments}
#' \item{mod}{Summary of regression}
#' \item{dat}{Original data used for MR-GRIP}
#' }
mr_grip <- function(b_exp, b_out, se_exp, se_out, parameters) {
  stopifnot(length(b_exp) == length(b_out))
  stopifnot(length(se_exp) == length(se_out))
  stopifnot(length(b_exp) == length(se_out))

  nulllist <- list(
    b = NA,
    se = NA,
    pval = NA,
    nsnp = NA,
    Q = NA,
    Q_df = NA,
    Q_pval = NA,
    mod = NA,
    smod = NA,
    dat = NA
  )
  if (
    sum(!is.na(b_exp) && !is.na(b_out) && !is.na(se_exp) && !is.na(se_out)) < 3
  ) {
    return(nulllist)
  }

  dat <- data.frame(
    b_out = b_out,
    b_exp = b_exp,
    se_exp = se_exp,
    se_out = se_out
  )
  grip_out <- b_out * b_exp
  grip_exp <- b_exp^2
  # GRIP regression.  Includes intercept.  Weights designed to replicate IVW under no intercept.
  mod <- stats::lm(grip_out ~ grip_exp, weights = 1 / (grip_exp * se_out^2))
  smod <- summary(mod)
  b <- stats::coefficients(smod)[2, 1]
  se <- stats::coefficients(smod)[2, 2]
  b.adj <- NA
  se.adj <- NA
  pval.adj <- NA
  pval <- 2 * stats::pt(abs(b / se), length(b_exp) - 2L, lower.tail = FALSE)
  return(list(
    b = b,
    se = se,
    pval = pval,
    b.adj = b.adj,
    se.adj = se.adj,
    pval.adj = pval.adj,
    nsnp = length(b_exp),
    mod = smod,
    dat = dat
  ))
}
