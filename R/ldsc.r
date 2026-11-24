Ghuber <- function(u, k = 30, deriv = 0) {
  if (!deriv) {
    return(pmin(1, k / abs(u)))
  } else {
    return(abs(u) <= k)
  }
}

#' Univariate LDSC
#'
#' Imported here to help estimate sample overlap between studies
#'
#' @param Z summary Z-statistics for M variants
#' @param r2 average reference LD scores for M variants
#' @param N GWAS sample size for each variant (could be different across variants)
#' @param W variant weight
#'
#' @keywords internal
#' @return model fit
ldsc_h2_internal <- function(Z, r2, N, W = NULL) {
  if (is.null(W)) {
    W <- rep(1, length(Z))
  }
  tau <- (mean(Z^2) - 1) / mean(N * r2)
  Wv <- 1 / (1 + tau * N * r2)^2
  id <- which(Z^2 > 30)
  if (length(id) > 0) {
    Wv[id] <- sqrt(Wv[id])
  }
  mod <- MASS::rlm(I(Z^2) ~ I(N * r2), weight = W * Wv, psi = Ghuber, k = 30)
  return(summary(mod))
}


#' Bivariate LDSC
#'
#' Imported here to help estimate sample overlap between studies
#'
#' @param Zs Mx2 matrix of summary Z-statistics for M variants from two GWAS
#' @param r2 average reference LD scores for M variants
#' @param N1 sample size for the 1st GWAS
#' @param N2 sample size for the 2nd GWAS
#' @param Nc overlapped sample size between the two GWAS
#' @param W variant weight
#' @param h1 hsq for trait 1
#' @param h2 hsq for trait 2
#'
#' @return List of models
#' @references
#' Bulik-Sullivan,B.K. et al. (2015) An atlas of genetic correlations across human diseases and traits. Nat. Genet. 47, 1236–1241.
#'
#' Guo,B. and Wu,B. (2018) Principal component based adaptive association test of multiple traits using GWAS summary statistics. bioRxiv 269597; doi: 10.1101/269597
#'
#' Gua,B. and Wu,B. (2019) Integrate multiple traits to detect novel trait-gene association using GWAS summary data with an adaptive test approach. Bioinformatics. 2019 Jul 1;35(13):2251-2257. doi: 10.1093/bioinformatics/bty961.
#'
#' https://github.com/baolinwu/MTAR
#' @keywords internal
ldsc_rg_internal <- function(Zs, r2, h1, h2, N1, N2, Nc = 0, W = NULL) {
  if (is.null(W)) {
    W <- rep(1, length(r2))
  }

  Y <- Zs[, 1] * Zs[, 2]

  X <- (sqrt(N1) * sqrt(N2) + sqrt(Nc / N1 * N2)) * r2
  N1r2 <- N1 * r2
  N2r2 <- N2 * r2
  r0 <- 0

  ## 1st round
  if (any(Nc > 0)) {
    rcf <- as.vector(MASS::rlm(Y ~ X, psi = Ghuber)$coef)
    r0 <- rcf[1]
    gv <- rcf[-1]
  } else {
    gv <- as.vector(MASS::rlm(Y ~ X - 1, psi = Ghuber)$coef)
  }

  ## 2nd round
  Wv <- 1 / ((h1 * N1r2 + 1) * (h2 * N2r2 + 1) + (X * gv + r0)^2)
  id <- which(abs(Zs[, 1] * Zs[, 2]) > 30)
  if (length(id) > 0) {
    Wv[id] <- sqrt(Wv[id])
  }

  if (any(Nc > 0)) {
    rcf <- MASS::rlm(Y ~ X, weight = W * Wv, psi = Ghuber, k = 30)
  } else {
    rcf <- MASS::rlm(Y ~ X - 1, weight = W * Wv, psi = Ghuber, k = 30)
  }
  return(summary(rcf))
}


#' Univariate LDSC
#'
#' Imported here to help estimate sample overlap between studies.
#'
#' **Not available yet.** This function stops with an error of class `twosamplemr_not_yet_available`: it needs HapMap3 variants with their LD scores and allele frequencies, to infer ancestry, which opengwasr does not provide yet (opengwas/opengwasr#26).
#'
#' @param id ID to analyse
#' @param ancestry ancestry of traits 1 and 2 (AFR, AMR, EAS, EUR, SAS) or 'infer' (default) in which case it will try to guess based on allele frequencies
#' @param snpinfo Not used until the function is available again
#' @param splitsize How many SNPs to extract at one time. Default=`20000`
#'
#' @export
#' @return model fit
#' @references
#' Bulik-Sullivan,B.K. et al. (2015) An atlas of genetic correlations across human diseases and traits. Nat. Genet. 47, 1236–1241.
#'
#' Guo,B. and Wu,B. (2018) Principal component based adaptive association test of multiple traits using GWAS summary statistics. bioRxiv 269597; doi: 10.1101/269597
#'
#' Gua,B. and Wu,B. (2019) Integrate multiple traits to detect novel trait-gene association using GWAS summary data with an adaptive test approach. Bioinformatics. 2019 Jul 1;35(13):2251-2257. doi: 10.1093/bioinformatics/bty961.
#'
#' <https://github.com/baolinwu/MTAR>
ldsc_h2 <- function(id, ancestry = "infer", snpinfo = NULL, splitsize = 20000) {
  ldsc_not_yet_available("ldsc_h2")
}

#' Bivariate LDSC
#'
#' Imported here to help estimate sample overlap between studies.
#'
#' **Not available yet.** This function stops with an error of class `twosamplemr_not_yet_available`: it needs HapMap3 variants with their LD scores and allele frequencies, to infer ancestry, which opengwasr does not provide yet (opengwas/opengwasr#26).
#'
#' @param id1 ID 1 to analyse
#' @param id2 ID 2 to analyse
#' @param ancestry ancestry of traits 1 and 2 (AFR, AMR, EAS, EUR, SAS) or 'infer' (default) in which case it will try to guess based on allele frequencies
#' @param snpinfo Not used until the function is available again
#' @param splitsize How many SNPs to extract at one time. Default=`20000`
#'
#' @export
#' @return model fit
ldsc_rg <- function(id1, id2, ancestry = "infer", snpinfo = NULL, splitsize = 20000) {
  ldsc_not_yet_available("ldsc_rg")
}

ldsc_not_yet_available <- function(fn) {
  not_yet_available(
    fn,
    "HapMap3 variants with their LD scores and allele frequencies",
    "opengwas/opengwasr#26"
  )
}
