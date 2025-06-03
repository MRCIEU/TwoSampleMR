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
#' \item{b_i}{Intercept}
#' \item{se_i}{Standard error of intercept}
#' \item{pval_i}{p-value of intercept}
#' \item{b.adj}{MR estimate adjusting for weak instruments}
#' \item{se.adj}{Standard error adjusting for weak instruments}
#' \item{pval.adj}{p-value adjusting for weak instruments}
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
    b_i = NA,
    se_i = NA,
    pval_i = NA,
    b.adj = NA,
    se.adj = NA,
    pval.adj = NA,
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
  grip_weights <- 1 / (b_exp^2 * se_out^2)

  # GRIP regression.  Includes intercept.  Weights under NOME assumption.
  mod <- stats::lm(grip_out ~ grip_exp, weights = grip_weights)
  smod <- summary(mod)
  b <- stats::coefficients(smod)[2, 1]
  se <- stats::coefficients(smod)[2, 2]
  b_i <- stats::coefficients(smod)[1, 1]
  se_i <- stats::coefficients(smod)[1, 2]

  # Weak instrument adjustment
  grip_weights <- 1 / (se_out^2)
  numer <- sum(grip_weights) *
    sum(grip_weights * b_out * b_exp * (b_exp^2 - 3 * se_exp^2)) -
    sum(grip_weights * b_out * b_exp) * sum(grip_weights * (b_exp^2 - se_exp^2))
  denom <- sum(grip_weights) *
    sum(grip_weights * (b_exp^4 - 6 * b_exp^2 * se_exp^2 + 3 * se_exp^4)) -
    (sum(grip_weights * (b_exp^2 - se_exp^2)))^2
  b.adj <- numer / denom

  var_out <- mean((b_out - (b * grip_exp + b_i) / b_exp)^2) * b_exp^2
  numer <- sum(grip_weights)^2 *
    sum(grip_weights^2 * var_out * (b_exp^2 - 3 * se_exp^2)^2) +
    sum(grip_weights^2 * var_out) * sum(grip_weights * (b_exp^2 - se_exp^2))^2 -
    2 *
      sum(grip_weights) *
      sum(grip_weights * (b_exp^2 - se_exp^2)) *
      sum(grip_weights^2 * (b_exp^2 - se_exp^2) * var_out)
  se.adj <- sqrt(numer) / denom

  pval <- 2 * stats::pt(abs(b / se), length(b_exp) - 2, lower.tail = FALSE)
  pval_i <- 2 *
    stats::pt(abs(b_i / se_i), length(b_exp) - 2, lower.tail = FALSE)
  pval.adj <- 2 *
    stats::pt(abs(b.adj / se.adj), length(b_exp) - 2, lower.tail = FALSE)
  return(list(
    b = b,
    se = se,
    pval = pval,
    b_i = b_i,
    se_i = se_i,
    pval_i = pval_i,
    b.adj = b.adj,
    se.adj = se.adj,
    pval.adj = pval.adj,
    nsnp = length(b_exp),
    mod = mod,
    smod = smod,
    dat = dat
  ))
}
