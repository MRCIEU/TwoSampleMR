# The one place where opengwasr's data frames become TwoSampleMR's.
#
# Everything that calls opengwasr, and everything that maps its columns to TwoSampleMR's formats,
# is here: the public functions (extract_instruments(), extract_outcome_data(), available_outcomes(),
# add_metadata(), mv_extract_exposures()) check their arguments and call these.

# --- Messages -------------------------------------------------------------------------------------

# A message of class `twosamplemr_deprecated_argument`, so that a script (or a test) can tell it from
# any other message.
deprecated_argument_message <- function(text) {
  message(structure(
    class = c("twosamplemr_deprecated_argument", "message", "condition"),
    list(message = paste0(text, "\n"), call = NULL)
  ))
}

# `opengwas_jwt` is read from the Source now.
deprecate_opengwas_jwt <- function(fn) {
  deprecated_argument_message(paste0(
    "`opengwas_jwt` is no longer used by `",
    fn,
    "()` and is ignored. ",
    "Credentials come from the Source: `opengwasr::og_api()` reads `OPENGWAS_JWT`, ",
    "or pass your own Source with `src`."
  ))
}

# The proxy arguments of extract_outcome_data() and mv_extract_exposures(), named by whichever ones
# the user passed.
deprecate_proxy_arguments <- function(passed, fn = "extract_outcome_data") {
  deprecated_argument_message(paste0(
    "The ",
    paste0("`", passed, "`", collapse = ", "),
    " argument",
    if (length(passed) > 1) "s" else "",
    " of `",
    fn,
    "()` no longer ",
    if (length(passed) > 1) "have" else "has",
    " any effect. LD proxies are not looked up: the Reference-Completed Store Releases that ",
    "opengwasr reads already hold Imputed Associations for most variants a study did not report. ",
    "They come back as the variant itself, with `association_status.outcome` set to \"imputed\" ",
    "(see `harmonise_data()` for the count)."
  ))
}

deprecate_splitsize <- function() {
  deprecated_argument_message(paste0(
    "`splitsize` is no longer needed by `extract_outcome_data()` and is ignored: ",
    "opengwasr splits large queries itself."
  ))
}

# --- Analyses -------------------------------------------------------------------------------------

# Identifiers as ieugwasr::legacy_ids() wrote them: a bare number is an `ieu-a-` id (or `met-a-`,
# `met-b-`, `met-c-` in the ranges those batches took over), and `UKB-a:` style prefixes are
# `ukb-a-`. Anything else is passed to opengwasr as it is.
legacy_analysis_ids <- function(x) {
  x <- as.character(x)
  bare <- !is.na(x) & !grepl("\\D", x)
  x[bare] <- paste0("ieu-a-", x[bare])
  for (prefix in c("UKB-a:", "UKB-b:", "UKB-c:", "IEU-a:")) {
    index <- startsWith(x, prefix)
    new <- paste0(tolower(sub(":", "", prefix)), "-")
    x[index] <- sub(prefix, new, x[index], fixed = TRUE)
  }
  to_met <- function(x, range, batch) {
    index <- x %in% paste0("ieu-a-", range)
    x[index] <- sub("ieu-a-", batch, x[index], fixed = TRUE)
    x
  }
  x <- to_met(x, 303:754, "met-a-")
  x <- to_met(x, 119:269, "met-b-")
  to_met(x, 838:960, "met-c-")
}

# The canonical ids of the Analyses `ids` name. An id that names no Analysis is warned about and
# left out (as TwoSampleMR always did), unless none is known, which is the opengwasr error.
resolve_analyses <- function(src, ids) {
  ids <- unique(legacy_analysis_ids(ids))
  info <- tryCatch(
    opengwasr::og_analyses(src, ids),
    opengwasr_analysis_not_found = function(e) {
      known <- setdiff(ids, e$ids)
      if (length(known) == 0) {
        stop(e)
      }
      warning(structure(
        class = c("twosamplemr_unknown_analyses", "warning", "condition"),
        list(
          message = paste0(
            "Left out because they name no Analysis: ",
            paste(e$ids, collapse = ", "),
            "."
          ),
          call = NULL,
          ids = e$ids
        )
      ))
      opengwasr::og_analyses(src, known)
    }
  )
  info$id
}

# The catalogue, as `available_outcomes()` returns it.
catalogue <- function(src) {
  a <- opengwasr::og_analyses(src)
  a$trait <- a$analysis_label
  a
}

# The catalogue for the functions that merge study characteristics into MR results. The old
# catalogue's `subcategory` and `nsnp` are not in the new one: results are not grouped by
# subcategory, and the number of SNPs on the array is not reported.
catalogue_for_results <- function() {
  ao <- available_outcomes()
  if (!"subcategory" %in% names(ao)) {
    message(
      "The OpenGWAS catalogue no longer has `subcategory` or `nsnp`, so outcomes are not ",
      "grouped by subcategory and the number of SNPs on the array is not reported."
    )
    ao$subcategory <- ""
  }
  ao
}

# --- Exposure data --------------------------------------------------------------------------------

# Top hits of `analyses` (canonical ids) in TwoSampleMR's exposure format.
opengwasr_instruments <- function(src, analyses, p1, clump, r2, kb) {
  d <- opengwasr::og_tophits(
    src,
    analyses,
    threshold = p1,
    limit = 10000,
    clump = as.logical(clump),
    r2 = r2,
    kb = kb
  )
  if (nrow(d) == 0) {
    return(NULL)
  }
  d$SNP <- variant_name(d)
  # Case and control counts of a continuous trait are all NA, and not worth a column.
  for (column in c("ncase", "ncontrol")) {
    if (all(is.na(d[[column]]))) d[[column]] <- NULL
  }
  d$phenotype <- paste0(d$analysis_label, " || id:", d$id)
  key <- paste(d$id, tolower(d$SNP))
  res <- format_data(
    d,
    type = "exposure",
    snps = NULL,
    phenotype_col = "phenotype",
    snp_col = "SNP",
    chr_col = "chromosome",
    pos_col = "position",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col = "n",
    ncase_col = "ncase",
    ncontrol_col = "ncontrol",
    min_pval = 1e-200,
    id_col = "id"
  )
  index <- match(paste(res$id.exposure, res$SNP), key)
  # format_data() lower-cases SNP, which would break an ALID.
  alid <- is.na(d$rsid[index])
  res$SNP[alid] <- d$alid[index][alid]
  res <- carry_columns(
    res,
    d,
    index,
    c("effect_scale", "association_status", "clump_ld", "clump_ld_choice", "clump_r2", "clump_kb"),
    "exposure"
  )
  res$data_source.exposure <- "igd"
  res
}

# --- Outcome data ---------------------------------------------------------------------------------

# Associations of `snps` in `analyses` (canonical ids) in TwoSampleMR's outcome format.
opengwasr_outcomes <- function(src, snps, analyses) {
  d <- opengwasr::og_associations(src, snps, analyses)
  if (is.null(nrow(d)) || nrow(d) == 0) {
    return(NULL)
  }
  d <- format_d(d, snps)
  if (nrow(d) > 0) {
    d$data_source.outcome <- "igd"
    return(d)
  }
  NULL
}

#' Format the table that opengwasr returns
#'
#' @param d Data frame from [opengwasr::og_associations()]
#' @param snps The SNPs that were asked for, so that a variant asked for by ALID keeps its ALID.
#' @return Data frame
#' @keywords internal
format_d <- function(d, snps = NULL) {
  d1 <- data.frame(
    SNP = variant_name(d, snps),
    chr = as.character(d$chromosome),
    pos = as.character(d$position),
    beta.outcome = as.numeric(d$beta),
    se.outcome = as.numeric(d$se),
    samplesize.outcome = as.numeric(d$n),
    pval.outcome = as.numeric(d$pval),
    eaf.outcome = as.numeric(d$eaf),
    effect_allele.outcome = as.character(d$effect_allele),
    other_allele.outcome = as.character(d$other_allele),
    outcome = as.character(d$analysis_label),
    id.outcome = as.character(d$id),
    stringsAsFactors = FALSE
  )
  if (!all(is.na(d$ncase))) d1$ncase.outcome <- as.numeric(d$ncase)
  if (!all(is.na(d$ncontrol))) d1$ncontrol.outcome <- as.numeric(d$ncontrol)
  d1 <- carry_columns(d1, d, seq_len(nrow(d)), c("effect_scale", "association_status"), "outcome")
  d <- d1

  if (nrow(d) == 0) {
    message("No matches")
    return(d)
  }

  d$originalname.outcome <- d$outcome
  d$outcome.deprecated <- paste0(d$outcome, " || ", d$consortium.outcome, " || ", d$year.outcome)
  d$outcome <- paste0(d$outcome, " || id:", d$id.outcome)

  rem <- is.na(d$beta.outcome) & is.na(d$pval.outcome)
  d <- subset(d, !rem)

  # For any that have missing SE but available beta and pval, infer SE
  index <- (is.na(d$se.outcome) | d$se.outcome == 0) &
    (!is.na(d$beta.outcome) & !is.na(d$pval.outcome))
  if (any(index)) {
    d$se.outcome[index] <- get_se(d$beta.outcome[index], d$pval.outcome[index])
  }

  d <- cleanup_outcome_data(d)

  mrcols <- c("beta.outcome", "se.outcome", "effect_allele.outcome")
  d$mr_keep.outcome <- complete.cases(d[, mrcols])
  if (!all(d$mr_keep.outcome)) {
    missinginfosnps <- paste(subset(d, !mr_keep.outcome)$SNP, collapse = " ")
    warning(
      "The following SNP(s) are missing required information for the MR tests and will be excluded: ",
      missinginfosnps
    )
  }
  if (!any(d$mr_keep.outcome)) {
    warning(
      "None of the provided SNPs can be used for MR analysis, they are missing required information."
    )
  }

  return(d)
}

# --- Metadata -------------------------------------------------------------------------------------

# What each metadata field of add_metadata() is called in `opengwasr::og_analyses()`.
analyses_names <- c(
  sample_size = "sample_size",
  ncase = "n_cases",
  ncontrol = "n_controls",
  effect_scale = "stored_effect_scale",
  sd = "original_sd"
)

# The `units` label that effect scales had before: `units.exposure` and `units.outcome` are read by
# add_rsq(), standardise_units(), get_rsq() and run_mrmix().
units_from_effect_scale <- function(effect_scale) {
  units <- c(sd = "SD", log_or = "log odds")[effect_scale]
  names(units) <- NULL
  unknown <- unique(effect_scale[is.na(units) & !is.na(effect_scale)])
  if (length(unknown) > 0) {
    message(
      "There is no `units` label for the effect scale ",
      paste0("\"", unknown, "\"", collapse = ", "),
      ", so it is NA."
    )
  }
  units
}

# Metadata of the Analyses `id` (as the data name them, which may be legacy ids), with columns named
# `<field>.<what>` as add_metadata() adds them. NULL where none of the ids is an Analysis.
opengwasr_metadata <- function(src, id, what, cols) {
  info <- tryCatch(
    opengwasr::og_analyses(src, id),
    opengwasr_analysis_not_found = function(e) {
      known <- setdiff(id, e$ids)
      if (length(known) == 0) data.frame() else opengwasr::og_analyses(src, known)
    }
  )
  if (nrow(info) == 0) {
    message(what, ": none of the IDs found in database")
    return(NULL)
  }
  # Match on the id the data carry, which may be a legacy id.
  info$id <- info$requested_id
  for (col in cols) {
    from <- if (col %in% names(analyses_names)) analyses_names[[col]] else col
    info[[col]] <- if (from %in% names(info)) info[[from]] else NA
  }
  if ("effect_scale" %in% cols) {
    info$units <- units_from_effect_scale(info$effect_scale)
    cols <- c(cols, "units")
  }
  info <- subset(info, select = c("id", cols))
  names(info) <- paste0(names(info), ".", what)
  names(info)[names(info) == paste0("sample_size.", what)] <- paste0("samplesize.", what)
  info
}

# --- Helpers --------------------------------------------------------------------------------------

# The columns of an opengwasr result that TwoSampleMR's formats have no slot for, copied across by
# Analysis and variant after format_data() has dropped them. `from` is the opengwasr data frame,
# `to` the formatted one, and `index` says which row of `from` each row of `to` came from.
carry_columns <- function(to, from, index, columns, type) {
  for (column in columns) {
    if (column %in% names(from)) {
      to[[paste0(column, ".", type)]] <- from[[column]][index]
    }
  }
  to
}

# opengwasr names variants `rsid`, which is NA for a variant without an rsID; `alid` always exists.
variant_name <- function(d, requested = NULL) {
  name <- ifelse(is.na(d$rsid), d$alid, d$rsid)
  if (!is.null(requested)) {
    # Someone who asked for an ALID gets that ALID back, even where the Analysis also knows an rsID.
    name[d$alid %in% requested] <- d$alid[d$alid %in% requested]
  }
  as.character(name)
}
