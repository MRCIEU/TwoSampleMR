#' Add meta data to extracted data
#'
#' Previously the meta data was returned alongside association information.
#' This is mostly unnecessary as it is needlessly repeating information.
#' This is a convenience function that reinstates that information.
#' Can be applied to either exposure data, outcome data, or harmonised data
#'
#' The meta data come from opengwasr (see [opengwasr::og_analyses()]). A field is added as
#' `<field>.exposure` or `<field>.outcome`, only where the data hold no value for it already.
#' The fields are named as they always were, and each one is read from the Analysis metadata
#' column that now holds it:
#'
#' * `sample_size` is `sample_size`, added as `samplesize`;
#' * `ncase` and `ncontrol` are `n_cases` and `n_controls`;
#' * `effect_scale` is `stored_effect_scale`: `"sd"` for effects per SD of the trait, `"log_or"`
#'   for log odds ratios;
#' * `sd` is `original_sd`, the SD of the trait in its original units where OpenGWAS recorded one;
#' * `unit` no longer exists, because effects are per SD, not in a study's own units. Asking for it
#'   gives a message, and `effect_scale` instead.
#'
#' Any other column of the Analysis metadata, such as `assigned_ancestry`, can be named too.
#' An id that names no Analysis is left without meta data, with a message.
#'
#' @param dat Either exposure data, outcome data or harmonised data
#' @param cols Which metadata fields to add. Default = `c("sample_size", "ncase", "ncontrol", "effect_scale", "sd")`
#' @param src A Source from opengwasr. The default, [opengwasr::og_api()], is made only when it is needed.
#'
#' @export
#' @return Data frame
add_metadata <- function(
  dat,
  cols = c("sample_size", "ncase", "ncontrol", "effect_scale", "sd"),
  src = opengwasr::og_api()
) {
  stopifnot(is.data.frame(dat))
  stopifnot("id.exposure" %in% names(dat) | "id.outcome" %in% names(dat))
  if ("unit" %in% cols) {
    deprecated_argument_message(paste0(
      "`unit` no longer exists in the OpenGWAS metadata: effects are per SD or log odds. ",
      "Adding `effect_scale` instead."
    ))
    cols <- unique(ifelse(cols == "unit", "effect_scale", cols))
  }
  # What each field is called in `opengwasr::og_analyses()`.
  analyses_names <- c(
    sample_size = "sample_size",
    ncase = "n_cases",
    ncontrol = "n_controls",
    effect_scale = "stored_effect_scale",
    sd = "original_sd"
  )
  get_info <- function(id, what = "exposure", cols) {
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
    info <- subset(info, select = c("id", cols))
    names(info) <- paste0(names(info), ".", what)
    names(info)[names(info) == paste0("sample_size.", what)] <- paste0("samplesize.", what)
    return(info)
  }

  order_col <- random_string()
  dat[[order_col]] <- seq_len(nrow(dat))
  if ("id.exposure" %in% names(dat)) {
    exposure_id <- unique(dat[["id.exposure"]])
    info <- get_info(id = exposure_id, what = "exposure", cols = cols)
    if (!is.null(info)) {
      for (x in names(info)) {
        if (!x %in% names(dat)) {
          dat[[x]] <- NA
        }

        for (id in unique(info[["id.exposure"]])) {
          dat[[x]][is.na(dat[[x]]) & dat[["id.exposure"]] == id] <- info[[x]][
            info[["id.exposure"]] == id
          ]
        }
      }
    }
  }

  if ("id.outcome" %in% names(dat)) {
    outcome_id <- unique(dat[["id.outcome"]])
    info <- get_info(id = outcome_id, what = "outcome", cols = cols)
    if (!is.null(info)) {
      for (x in names(info)) {
        if (!x %in% names(dat)) {
          dat[[x]] <- NA
        }

        for (id in unique(info[["id.outcome"]])) {
          dat[[x]][is.na(dat[[x]]) & dat[["id.outcome"]] == id] <- info[[x]][
            info[["id.outcome"]] == id
          ]
        }
      }
    }
  }

  dat <- dat[order(dat[[order_col]]), ]
  dat <- dat[, !names(dat) %in% order_col]
  # dat <- fix_ukb_d(dat)
  return(dat)
}
