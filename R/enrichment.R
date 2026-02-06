#' Fisher's combined test
#'
#' @param pval Vector of outcome p-values
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{MR estimate}
#' \item{se}{Standard error}
#' \item{pval}{p-value}
#' }
fishers_combined_test <- function(pval) {
  pval <- pval[is.finite(pval) & pval <= 1 & pval >= 0]
  index <- pval == 0
  if (any(index)) {
    warning("p-values of 0 are unreliable in Fisher's combined test.")
    pval[index] <- 1e-50
  }
  p <- stats::pchisq(-2 * sum(log(pval)), df = 2 * length(pval), lower.tail = FALSE)
  return(list(pval = p, nsnp = length(pval)))
}


#' Get list of available p-value enrichment methods
#'
#' @export
#' @return Data frame
enrichment_method_list <- function() {
  a <- list(
    list(
      obj = "fishers_combined_test",
      name = "Fisher's combined test",
      PubmedID = "",
      Description = ""
    )
  )
  a <- lapply(a, as.data.frame)
  a <- data.table::rbindlist(a, fill = TRUE, use.names = TRUE)
  data.table::setDF(a)
  a <- as.data.frame(lapply(a, as.character), stringsAsFactors = FALSE)
  return(a)
}


#' Perform enrichment analysis
#'
#'
#' @param dat Harmonised exposure and outcome data. Output from [harmonise_data()].
#' @param method_list List of methods to use in analysis. Default is \code{enrichment_method_list()$obj}.
#' See [enrichment_method_list()] for details.
#'
#' @export
#' @return data frame
enrichment <- function(dat, method_list = enrichment_method_list()$obj) {
  methl <- enrichment_method_list()
  combos <- unique(dat[, c("id.exposure", "id.outcome")])
  results <- lapply(seq_len(nrow(combos)), function(i) {
    x1 <- dat[dat$id.exposure == combos$id.exposure[i] & dat$id.outcome == combos$id.outcome[i], ]
    x <- subset(x1, !is.na(pval.outcome))
    if (nrow(x) == 0) {
      message(
        "No outcome p-values for analysis of '",
        x1$id.exposure[1],
        "' on '",
        x1$id.outcome[1],
        "'"
      )
      return(NULL)
    }
    res <- lapply(method_list, function(meth) {
      get(meth)(x$pval.outcome)
    })
    enrichment_tab <- data.frame(
      id.exposure = x1$id.exposure[1],
      id.outcome = x1$id.outcome[1],
      outcome = x$outcome[1],
      exposure = x$exposure[1],
      method = methl$name[methl$obj %in% method_list],
      nsnp = vapply(res, function(x) x$nsnp, numeric(1)),
      pval = vapply(res, function(x) x$pval, numeric(1)),
      stringsAsFactors = FALSE
    )
    enrichment_tab <- subset(enrichment_tab, !is.na(pval))
    return(enrichment_tab)
  })
  res <- data.table::rbindlist(results, fill = TRUE, use.names = TRUE)
  data.table::setDF(res)
  return(res)
}
