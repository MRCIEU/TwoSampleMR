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
#' When `effect_scale` is added, `units.exposure` and `units.outcome` are added too, derived from it
#' (`"sd"` is `"SD"`, `"log_or"` is `"log odds"`; any other scale gives `NA`, with a message), as the
#' functions that read them, such as [add_rsq()] and [standardise_units()], expect.
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
  get_info <- function(id, what, cols) opengwasr_metadata(src, id, what, cols)

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
