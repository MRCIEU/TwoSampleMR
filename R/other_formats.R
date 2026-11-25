#' Convert TwoSampleMR format to MendelianRandomization format
#'
#' The MendelianRandomization package offers MR methods that can
#' be used with the same data used in the TwoSampleMR package. This
#' function converts from the TwoSampleMR format to the MRInput class.
#'
#' @param dat Output from the [harmonise_data()] function.
#' @param get_correlations Default `FALSE`. If `TRUE` then extract the LD matrix for the SNPs from an OpenGWAS LD Store with [ld_matrix()]. From an LD Store, pairs below its Far Threshold are 0, and pairs it cannot measure (other chromosomes, beyond its far window) are set to 0 too, treated as independent, with a message giving how many.
#' @param pop If `get_correlations` is `TRUE` then the ancestry (or LD Store id) of the LD Store to use. The default, `NULL`, is the European LD Store on the API Source, and the Source's own LD resource on a Local Source (see [ld_matrix()]).
#' @param bfile If `get_correlations` is `TRUE`, the prefix of a PLINK fileset to use as the LD reference instead of an LD Store. Default = `NULL`
#' @param plink_bin Path to the plink binary for `bfile`. Default = `NULL`
#' @param src A Source from opengwasr, used with `get_correlations`: [opengwasr::og_api()], or a Local Source, [opengwasr::og_local()], whose own LD resource is then used. The default, [opengwasr::og_api()] (named as TwoSampleMR's), is made only when it is needed.
#'
#' @export
#' @return List of MRInput objects for each exposure/outcome combination
dat_to_MRInput <- function(
  dat,
  get_correlations = FALSE,
  pop = NULL,
  bfile = NULL,
  plink_bin = NULL,
  src = twosamplemr_source()
) {
  if (!requireNamespace("MendelianRandomization", quietly = TRUE)) {
    stop(
      "Package \"MendelianRandomization\" must be installed to use this function.",
      call. = FALSE
    )
  }
  combos <- unique(dat[, c("exposure", "outcome")])
  out <- lapply(seq_len(nrow(combos)), function(i) {
    x <- dat[dat$exposure == combos$exposure[i] & dat$outcome == combos$outcome[i], ]
    message("Converting:")
    message(" - exposure: ", x$exposure[1])
    message(" - outcome: ", x$outcome[1])
    if (get_correlations) {
      out <- ld_for_mr_input(x, pop, bfile, plink_bin, src)
      if (is.null(out)) {
        return(NULL)
      }
      x <- out$x
      ld <- out$ld

      MendelianRandomization::mr_input(
        bx = x$beta.exposure,
        bxse = x$se.exposure,
        by = x$beta.outcome,
        byse = x$se.outcome,
        exposure = x$exposure[1],
        outcome = x$outcome[1],
        snps = x$SNP,
        effect_allele = x$effect_allele.exposure,
        other_allele = x$other_allele.exposure,
        eaf = x$eaf.exposure,
        correlation = ld
      )
    } else {
      MendelianRandomization::mr_input(
        bx = x$beta.exposure,
        bxse = x$se.exposure,
        by = x$beta.outcome,
        byse = x$se.outcome,
        exposure = x$exposure[1],
        outcome = x$outcome[1],
        snps = x$SNP,
        effect_allele = x$effect_allele.exposure,
        other_allele = x$other_allele.exposure,
        eaf = x$eaf.exposure
      )
    }
  })
  names(out) <- paste(combos$exposure, combos$outcome, sep = ".")
  return(out)
}


#' Harmonise LD matrix against summary data
#'
#' LD matrix returns with rsid_ea_oa identifiers. Make sure that they are oriented to the same effect allele as the summary dataset. Summary dataset can be exposure dataset or harmonised dataset.
#'
#' @param x Exposure dataset or harmonised dataset
#' @param ld Output from [ld_matrix()]. A matrix named by ALID is renamed `<SNP>_<allele1>_<allele2>`, with `SNP` as in `x`.
#' @param src A Source from opengwasr, used only to find the ALIDs of rsIDs in `x` that have no position.
#'
#' @export
#' @return List of exposure dataset and harmonised LD matrix
harmonise_ld_dat <- function(x, ld, src = twosamplemr_source()) {
  ld <- name_ld_by_snp(x, ld, src)
  snpnames <- do.call(rbind, strsplit(rownames(ld), split = "_"))
  i1 <- snpnames[, 1] %in% x$SNP
  ld <- ld[i1, i1, drop = FALSE]
  snpnames <- snpnames[i1, , drop = FALSE]
  i2 <- x$SNP %in% snpnames[, 1]
  x <- x[i2, ]
  # Put the LD matrix in the same SNP order as x
  i3 <- match(x$SNP, snpnames[, 1])
  ld <- ld[i3, i3, drop = FALSE]
  snpnames <- snpnames[i3, , drop = FALSE]
  stopifnot(all(snpnames[, 1] == x$SNP))
  x$effect_allele.exposure <- as.character(x$effect_allele.exposure)
  x$other_allele.exposure <- as.character(x$other_allele.exposure)
  # Set1 x and ld alleles match
  snpnames <- data.frame(snpnames, stringsAsFactors = FALSE)
  snpnames <- merge(
    subset(x, select = c(SNP, effect_allele.exposure, other_allele.exposure)),
    snpnames,
    by.x = "SNP",
    by.y = "X1"
  )
  snpnames <- snpnames[match(x$SNP, snpnames$SNP), ]
  snpnames$keep <- (snpnames$X2 == snpnames$effect_allele.exposure &
    snpnames$X3 == snpnames$other_allele.exposure) |
    (snpnames$X3 == snpnames$effect_allele.exposure & snpnames$X2 == snpnames$other_allele.exposure)

  # What happens if everything is gone?
  if (nrow(x) == 0) {
    message(" - none of the SNPs could be aligned to the LD reference panel")
    return(NULL)
  }

  if (!all(snpnames$keep)) {
    message(
      " - the following SNPs could not be aligned to the LD reference panel: \n- ",
      paste(subset(snpnames, !keep)$SNP, collapse = "\n - ")
    )
  }

  snpnames$flip1 <- snpnames$X2 != snpnames$effect_allele.exposure
  x <- subset(x, SNP %in% snpnames$SNP)
  temp1 <- x$effect_allele.exposure[snpnames$flip1]
  temp2 <- x$other_allele.exposure[snpnames$flip1]
  x$beta.exposure[snpnames$flip1] <- x$beta.exposure[snpnames$flip1] * -1
  x$beta.outcome[snpnames$flip1] <- x$beta.outcome[snpnames$flip1] * -1
  x$effect_allele.exposure[snpnames$flip1] <- temp2
  x$other_allele.exposure[snpnames$flip1] <- temp1

  rownames(ld) <- snpnames$SNP
  colnames(ld) <- snpnames$SNP

  if (!all(snpnames$keep)) {
    message("Removing ", sum(!snpnames$keep), " variants due to harmonisation issues")
    ld <- ld[snpnames$keep, snpnames$keep, drop = FALSE]
    x <- x[snpnames$keep, ]
  }
  return(list(x = x, ld = ld))
}


#' Wrapper for MR-PRESSO
#'
#' See <https://github.com/rondolab/MR-PRESSO> for more details.
#'
#' @param dat Output from [harmonise_data()].
#' @param NbDistribution Number of bootstrap replications. The default is `1000`.
#' @param SignifThreshold Outlier significance threshold. The default is `0.05`.
#'
#' @export
#' @return List of results for every exposure/outcome combination
run_mr_presso <- function(dat, NbDistribution = 1000, SignifThreshold = 0.05) {
  dat <- subset(dat, mr_keep)
  d <- subset(
    dat,
    !duplicated(data.frame(id.exposure, id.outcome)),
    select = c(exposure, outcome, id.exposure, id.outcome)
  )
  res <- list()
  attributes(res)$id.exposure <- d$id.exposure
  attributes(res)$id.outcome <- d$id.outcome
  attributes(res)$exposure <- d$exposure
  attributes(res)$outcome <- d$outcome
  for (j in seq_len(nrow(d))) {
    x <- subset(dat, exposure == d$exposure[j] & outcome == d$outcome[j])
    message(x$exposure[1], " - ", x$outcome[1])
    res[[j]] <- MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = x,
      NbDistribution = NbDistribution,
      SignifThreshold = SignifThreshold
    )
  }
  return(res)
}

#' Convert dat to RadialMR format
#'
#' Creates a list of RadialMR format datasets for each exposure-outcome pair.
#'
#' @param dat Output from [harmonise_data()].
#'
#' @export
#' @return List of RadialMR format datasets
dat_to_RadialMR <- function(dat) {
  combos <- unique(dat[, c("exposure", "outcome")])
  out <- lapply(seq_len(nrow(combos)), function(i) {
    x <- dat[dat$exposure == combos$exposure[i] & dat$outcome == combos$outcome[i], ]
    message("Converting:")
    message(" - exposure: ", x$exposure[1])
    message(" - outcome: ", x$outcome[1])
    d <- subset(x, mr_keep)
    d <- RadialMR::format_radial(
      d$beta.exposure,
      d$beta.outcome,
      d$se.exposure,
      d$se.outcome,
      RSID = d$SNP
    )
    return(d)
  })
  names(out) <- paste(combos$exposure, combos$outcome, sep = ".")
  return(out)
}

#' Radial IVW analysis
#'
#' @param b_exp Vector of genetic effects on exposure.
#' @param b_out Vector of genetic effects on outcome.
#' @param se_exp Standard errors of genetic effects on exposure.
#' @param se_out Standard errors of genetic effects on outcome.
#' @param parameters List of parameters.
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{causal effect estimate}
#' \item{se}{standard error}
#' \item{pval}{p-value}
#' }
mr_ivw_radial <- function(b_exp, b_out, se_exp, se_out, parameters = default_parameters()) {
  if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2) {
    return(list(b = NA, se = NA, pval = NA, nsnp = NA))
  }
  d <- RadialMR::format_radial(
    BXG = b_exp,
    BYG = b_out,
    seBXG = se_exp,
    seBYG = se_out,
    RSID = seq_along(b_exp)
  )
  out <- RadialMR::ivw_radial(d, alpha = 0.05, weights = 3)
  b <- out$coef[1, 1]
  se <- out$coef[1, 2]
  pval <- 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
  Q_df <- out$df
  Q <- out$qstatistic
  Q_pval <- stats::pchisq(Q, Q_df, lower.tail = FALSE)
  return(list(
    b = b,
    se = se,
    pval = pval,
    nsnp = length(b_exp),
    Q = Q,
    Q_df = Q_df,
    Q_pval = Q_pval
  ))
}


#' Perform MRMix analysis on harmonised dat object
#'
#' See <https://github.com/gqi/MRMix> for more details.
#'
#' @param dat Output from [harmonise_data()]. Ensures that no eaf.exposure values are missing.
#' Missing eaf.exposure values are replaced by eaf.outcome values, and SNPs with an eaf.exposure that is still missing or is not strictly between 0 and 1 are removed with a warning.
#'
#' @export
#' @return List of results, with one list item for every exposure/outcome pair in dat object
run_mrmix <- function(dat) {
  combos <- unique(dat[, c("id.exposure", "id.outcome")])
  out <- lapply(seq_len(nrow(combos)), function(i) {
    x <- dat[dat$id.exposure == combos$id.exposure[i] & dat$id.outcome == combos$id.outcome[i], ]
    message("Analysing ", x$id.exposure[1], " against ", x$id.outcome[1])
    if (grepl("log odds", x$units.exposure[1])) {
      xunit <- "binary"
    } else if (grepl("^SD", x$units.exposure[1])) {
      xunit <- "n"
    } else {
      xunit <- "continuous"
    }
    if (grepl("log odds", x$units.outcome[1])) {
      yunit <- "binary"
    } else if (grepl("^SD", x$units.outcome[1])) {
      yunit <- "n"
    } else {
      yunit <- "continuous"
    }
    index <- is.na(x$eaf.exposure)
    if (any(index)) {
      warning(paste0(
        x$id.exposure[1],
        ".",
        x$id.outcome[1],
        " - Some eaf.exposure values are missing, using eaf.outcome in their place"
      ))
      x$eaf.exposure[index] <- x$eaf.outcome[index]
    }
    # standardize() scales by sqrt(2 * eaf * (1 - eaf)), so SNPs with an eaf of
    # 0 or 1 would give standard errors of 0
    index <- is.na(x$eaf.exposure) | x$eaf.exposure <= 0 | x$eaf.exposure >= 1
    if (any(index)) {
      warning(paste0(
        x$id.exposure[1],
        ".",
        x$id.outcome[1],
        " - Removing ",
        sum(index),
        " SNP(s) with eaf.exposure missing or not strictly between 0 and 1"
      ))
      x <- x[!index, ]
    }
    if (nrow(x) == 0) {
      warning(paste0(combos$id.exposure[i], ".", combos$id.outcome[i], " - No SNPs remaining"))
      return(NULL)
    }
    l <- MRMix::standardize(
      x$beta.exposure,
      x$beta.outcome,
      x$se.exposure,
      x$se.outcome,
      xunit,
      yunit,
      x$samplesize.exposure,
      x$samplesize.outcome,
      x$eaf.exposure
    )
    x$beta.exposure <- l$betahat_x_std
    x$beta.outcome <- l$betahat_y_std
    x$se.exposure <- l$sx_std
    x$se.outcome <- l$sy_std
    MRMix::MRMix(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome)
  })
  names(out) <- paste(combos$id.exposure, combos$id.outcome, sep = ".")
  out
}


# An LD Store names its matrix by ALID. This renames it `<SNP>_<allele1>_<allele2>`, with `SNP` as
# `x` names it, which is how a matrix made from a PLINK reference is named; variants of `x` that are
# not in the matrix are not renamed, and those of the matrix that are not in `x` are dropped.
name_ld_by_snp <- function(x, ld, src) {
  alid <- grepl(alid_pattern, rownames(ld))
  if (!all(alid)) {
    return(ld)
  }
  at <- match(rownames(ld), dat_alids(x, src))
  parts <- do.call(rbind, strsplit(rownames(ld), ":", fixed = TRUE))
  name <- paste(x$SNP[at], parts[, 3], parts[, 4], sep = "_")
  ld <- ld[!is.na(at), !is.na(at), drop = FALSE]
  rownames(ld) <- colnames(ld) <- name[!is.na(at)]
  ld
}


# The LD matrix of the variants of `x` for MendelianRandomization, harmonised to `x`, or NULL if none
# can be. From an LD Store the pairs below its Far Threshold are 0 (`fill_below_t_far`), and the
# remaining unmeasured pairs are set to 0 too; a PLINK Reference's matrix is complete already.
ld_for_mr_input <- function(x, pop, bfile, plink_bin, src) {
  message(" - obtaining LD matrix")
  ld <- ld_matrix(
    unique(x$SNP),
    pop = pop,
    bfile = bfile,
    plink_bin = plink_bin,
    fill_below_t_far = is.null(bfile),
    src = src
  )
  out <- harmonise_ld_dat(x, ld, src)
  if (is.null(out)) {
    return(NULL)
  }
  out$ld <- complete_ld(out$ld, plink = !is.null(bfile))
  out
}

# Pairs of variants that an LD Store does not hold (other chromosomes, or beyond its far window) are
# unmeasured: they are treated as independent, so they are set to 0, and the user is told how many.
complete_ld <- function(ld, plink = FALSE) {
  missing <- is.na(ld)
  if (plink || !any(missing)) {
    return(ld)
  }
  pairs <- sum(missing[upper.tri(missing)])
  message(
    " - ",
    pairs,
    " pair(s) of variants were set to 0 in the LD matrix: the LD Store does not hold them ",
    "because the variants are on different chromosomes or further apart than its far window, so ",
    "they are treated as independent."
  )
  ld[missing] <- 0
  ld
}
