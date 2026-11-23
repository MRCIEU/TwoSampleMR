# Helpers for the functions that read OpenGWAS through opengwasr.

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

# The proxy arguments of extract_outcome_data() and mv_extract_exposures(), named by whichever ones the user passed.
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
