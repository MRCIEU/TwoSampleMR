#' MR-GRIP: a modified MR-Egger model with the Genotype Recoding Invariance Property
#'
#' This implements the modified MR-Egger model with the Genotype Recoding Invariance Property (MR-GRIP) due to Dudbridge and Bowden et al. (2025).
#' It is well known that the results of MR-Egger are sensitive to which alleles are designated as the effect alleles.
#' A pragmatic convention is to orient all SNPs to have positive effects on the exposure, which has some advantages in interpretation but also brings some philosophical limitations.
#' The MR-GRIP model is a modification to the MR-Egger model in which each term is multiplied by the genotype-phenotype associations.
#' This makes each term in the model invariant to allele coding.
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
#' \item{b.wi}{MR estimate adjusting for weak instruments}
#' \item{se.wi}{Standard error adjusting for weak instruments}
#' \item{pval.wi}{p-value adjusting for weak instruments}
#' \item{mod}{Summary of regression}
#' \item{dat}{Original data used for MR-GRIP}
#' }
mr_grip <- function(b_exp, b_out, se_exp, se_out, parameters) {
  if (length(b_exp) != length(b_out)) stop("The lengths of b_exp and b_out are not equal.")
  if (length(se_exp) != length(se_out)) stop("The lengths of se_exp and se_out are not equal.")
  if (length(b_exp) != length(se_out)) stop("The lengths of b_exp and se_out are not equal.")

  nulllist <- list(
    b = NA,
    se = NA,
    pval = NA,
    nsnp = NA,
    mod = NA,
    smod = NA,
    dat = NA
  )
  if (
    sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3
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
  b.wi <- NA
  se.wi <- NA
  pval.wi <- NA
  pval <- 2 * stats::pt(abs(b / se), length(b_exp) - 2L, lower.tail = FALSE)
  return(list(
    b = b,
    se = se,
    pval = pval,
    b.wi = b.wi,
    se.wi = se.wi,
    pval.wi = pval.wi,
    nsnp = length(b_exp),
    mod = smod,
    dat = dat
  ))
}
