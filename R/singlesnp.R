#' Perform 2 sample MR on each SNP individually
#'
#' @param dat Output from [harmonise_data()].
#' @param parameters List of parameters. The default is `default_parameters()`.
#' @param single_method Function to use for MR analysis. The default is `"mr_wald_ratio"`.
#' @param all_method Functions to use for MR analysis. The default is `c("mr_ivw", "mr_egger_regression")`.
#'
#' @export
#' @return List of data frames
mr_singlesnp <- function(
  dat,
  parameters = default_parameters(),
  single_method = "mr_wald_ratio",
  all_method = c("mr_ivw", "mr_egger_regression")
) {
  if (!"samplesize.outcome" %in% names(dat)) {
    dat$samplesize.outcome <- NA
  }

  stopifnot("outcome" %in% names(dat))
  stopifnot("exposure" %in% names(dat))
  stopifnot("beta.exposure" %in% names(dat))
  stopifnot("beta.outcome" %in% names(dat))
  stopifnot("se.exposure" %in% names(dat))
  stopifnot("se.outcome" %in% names(dat))

  dat_dt <- data.table::as.data.table(dat)
  combos <- unique(dat_dt[, .(id.exposure, id.outcome)])

  # Pre-compute method names outside the loop
  method_list <- mr_method_list()
  method_names <- vapply(
    all_method,
    function(m) {
      paste0("All - ", method_list$name[method_list$obj == m])
    },
    character(1)
  )

  results <- lapply(seq_len(nrow(combos)), function(i) {
    exp_id <- combos$id.exposure[i]
    out_id <- combos$id.outcome[i]
    X <- dat_dt[id.exposure == exp_id & id.outcome == out_id]
    x <- X[mr_keep == TRUE]
    nsnp <- nrow(x)
    if (nsnp == 0) {
      x <- X[1, ]
      d <- data.frame(
        id.exposure = exp_id,
        id.outcome = out_id,
        SNP = "No available data",
        b = NA,
        se = NA,
        p = NA,
        samplesize = NA,
        outcome = x$outcome[1],
        exposure = x$exposure[1],
        stringsAsFactors = FALSE
      )
      return(d)
    }
    l <- lapply(1:nsnp, function(j) {
      with(
        x,
        get(single_method)(
          beta.exposure[j],
          beta.outcome[j],
          se.exposure[j],
          se.outcome[j],
          parameters
        )
      )
    })
    for (j in seq_along(all_method)) {
      l[[nsnp + j]] <- with(
        x,
        get(all_method[j])(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters)
      )
    }

    d <- data.frame(
      id.exposure = exp_id,
      id.outcome = out_id,
      SNP = c(as.character(x$SNP), method_names),
      b = vapply(l, function(y) y$b, numeric(1)),
      se = vapply(l, function(y) y$se, numeric(1)),
      p = vapply(l, function(y) y$pval, numeric(1)),
      samplesize = x$samplesize.outcome[1],
      stringsAsFactors = FALSE
    )
    d$outcome <- x$outcome[1]
    d$exposure <- x$exposure[1]
    return(d)
  })
  res <- data.table::rbindlist(results, fill = TRUE, use.names = TRUE)
  data.table::setDF(res)
  res <- subset(
    res,
    select = c(exposure, outcome, id.exposure, id.outcome, samplesize, SNP, b, se, p)
  )
  return(res)
}


#' Forest plot
#'
#' If the data frame contains a `category` column, SNPs will be coloured and
#' grouped by category in the forest plot (e.g., cluster assignments from
#' MR-Clust). See Figure 4 of Vabistsevits et al. (2024) for examples.
#'
#' @param singlesnp_results from [mr_singlesnp()].
#' @param exponentiate Plot on exponential scale. The default is `FALSE`.
#' @param category_colours Named character vector of colours for categories.
#'   Names should match the values in the `category` column. If `NULL`
#'   (the default), the Okabe-Ito colourblind-friendly palette is used.
#'
#' @references
#' Vabistsevits, M., Davey Smith, G., Richardson, T.G. et al.
#' Mammographic density mediates the protective effect of early-life body size
#' on breast cancer risk.
#' *Nature Communications*, **15**, 4021 (2024).
#' \doi{10.1038/s41467-024-48105-7}
#'
#' @export
#' @return List of plots
#' @examples
#' \dontrun{
#' # Basic forest plot
#' bmi_exp_dat <- extract_instruments(outcomes = "ieu-a-2")
#' chd_out_dat <- extract_outcome_data(
#'   snps = bmi_exp_dat$SNP, outcomes = "ieu-a-7"
#' )
#' dat <- harmonise_data(bmi_exp_dat, chd_out_dat)
#' res <- mr_singlesnp(dat)
#' mr_forest_plot(res)
#'
#' # Forest plot with RadialMR outliers
#' radial_dat <- dat_to_RadialMR(dat)
#' radial_res <- RadialMR::ivw_radial(radial_dat[[1]], alpha = 0.05, weights = 3)
#' outlier_snps <- radial_res$outliers$SNP
#' snp_rows <- !grepl("^All", res$SNP)
#' res$category <- NA_character_
#' res$category[snp_rows] <- ifelse(
#'   res$SNP[snp_rows] %in% outlier_snps, "Outlier", "Main"
#' )
#' mr_forest_plot(res)
#'
#' # With custom colours
#' mr_forest_plot(res, category_colours = c(Main = "grey50", Outlier = "red"))
#' }
mr_forest_plot <- function(singlesnp_results, exponentiate = FALSE, category_colours = NULL) {
  dat_dt <- data.table::as.data.table(singlesnp_results)
  combos <- unique(dat_dt[, .(id.exposure, id.outcome)])

  res <- lapply(seq_len(nrow(combos)), function(i) {
    exp_id <- combos$id.exposure[i]
    out_id <- combos$id.outcome[i]
    d <- as.data.frame(dat_dt[id.exposure == exp_id & id.outcome == out_id])
    if (sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
    levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
    am <- grep("All", d$SNP, value = TRUE)
    d$up <- d$b + 1.96 * d$se
    d$lo <- d$b - 1.96 * d$se
    d$tot <- 0.01
    d$tot[d$SNP %in% am] <- 1
    d$SNP <- as.character(d$SNP)

    has_category <- "category" %in% names(d)

    if (has_category) {
      d$category <- as.character(d$category)
      d$category[is.na(d$category)] <- "Uncategorized"

      # Order SNPs by overall effect size
      snp_rows <- d[!d$SNP %in% am, ]
      snp_rows <- snp_rows[order(snp_rows$b), ]
      nom <- snp_rows$SNP

      # Assign summary category to "All" rows
      d$category[d$SNP %in% am] <- "Summary"

      # Insert separator row
      sep_row <- d[nrow(d), ]
      sep_row$SNP <- ""
      sep_row$b <- NA
      sep_row$se <- NA
      sep_row$up <- NA
      sep_row$lo <- NA
      sep_row$category <- "Summary"
      d <- rbind(d, sep_row)

      d$SNP <- ordered(d$SNP, levels = c(am, "", nom))

      # Build colour palette for categories
      cats <- sort(unique(snp_rows$category))
      n_cats <- length(cats)
      if (!is.null(category_colours)) {
        cat_colours <- category_colours
      } else {
        # Okabe-Ito colourblind-friendly palette
        palette <- c(
          "#000000", "#0072B2", "#D55E00", "#009E73",
          "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#999999"
        )[seq_len(min(n_cats, 9))]
        if (n_cats > 9) {
          hues <- seq(15, 375, length.out = n_cats - 8)[seq_len(n_cats - 9)]
          palette <- c(palette, grDevices::hcl(h = hues, c = 100, l = 65))
        }
        cat_colours <- palette
        names(cat_colours) <- cats
      }
      cat_colours["Summary"] <- "black"

      # Black-coloured categories and Uncategorized get thin lines;
      # coloured categories and Summary get thick lines
      thin_cats <- c(
        names(cat_colours)[cat_colours %in% c("black", "#000000")],
        "Uncategorized"
      )
      thin_cats <- setdiff(thin_cats, "Summary")
      d$lw_cat <- ifelse(d$category %in% thin_cats, "thin", "thick")

    } else {
      # Original ordering: non-All SNPs sorted by effect size
      snp_rows <- d[!d$SNP %in% am, ]
      nom <- snp_rows$SNP[order(snp_rows$b)]

      # Insert separator row
      d <- rbind(d, d[nrow(d), ])
      d$SNP[nrow(d) - 1] <- ""
      d$b[nrow(d) - 1] <- NA
      d$up[nrow(d) - 1] <- NA
      d$lo[nrow(d) - 1] <- NA

      d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
    }

    xint <- 0
    if (exponentiate) {
      d$b <- exp(d$b)
      d$up <- exp(d$up)
      d$lo <- exp(d$lo)
      xint <- 1
    }

    # Build colour aesthetic and scale based on category mode
    if (has_category) {
      colour_aes <- ggplot2::aes(colour = category)
      colour_scale <- ggplot2::scale_colour_manual(values = cat_colours)
    } else {
      colour_aes <- ggplot2::aes(colour = as.factor(tot))
      colour_scale <- ggplot2::scale_colour_manual(values = c("black", "red"))
    }

    if (has_category) {
      lw_aes <- ggplot2::aes(linewidth = lw_cat)
      lw_scale <- ggplot2::scale_linewidth_manual(values = c(thin = 0.3, thick = 0.5))
    } else {
      lw_aes <- ggplot2::aes(linewidth = as.factor(tot))
      lw_scale <- ggplot2::scale_linewidth_manual(values = c(0.3, 1))
    }

    bar_aes <- utils::modifyList(
      utils::modifyList(ggplot2::aes(xmin = lo, xmax = up), lw_aes),
      colour_aes
    )

    if (utils::packageVersion("ggplot2") <= "3.5.2") {
      p <- ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) +
        ggplot2::geom_vline(xintercept = xint, linetype = "dotted") +
        ggplot2::geom_errorbarh(bar_aes, height = 0) +
        ggplot2::geom_point(colour_aes) +
        ggplot2::geom_hline(
          ggplot2::aes(yintercept = which(levels(SNP) %in% "")),
          colour = "grey"
        ) +
        colour_scale +
        lw_scale
    } else {
      p <- ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) +
        ggplot2::geom_vline(xintercept = xint, linetype = "dotted") +
        ggplot2::geom_errorbar(bar_aes, width = 0, orientation = "y") +
        ggplot2::geom_point(colour_aes) +
        ggplot2::geom_hline(
          ggplot2::aes(yintercept = which(levels(SNP) %in% "")),
          colour = "grey"
        ) +
        colour_scale +
        lw_scale
    }

    p +
      ggplot2::theme(
        legend.position = "none",
        axis.text.y = ggplot2::element_text(size = 8),
        axis.ticks.y = ggplot2::element_line(linewidth = 0),
        axis.title.x = ggplot2::element_text(size = 8)
      ) +
      ggplot2::labs(
        y = "",
        x = paste0("MR effect size for\n'", d$exposure[1], "' on '", d$outcome[1], "'")
      )
  })
  res
}


#' Density plot
#'
#' @param singlesnp_results from [mr_singlesnp()].
#' @param mr_results Results from [mr()].
#' @param exponentiate Plot on exponentiated scale. The default is `FALSE`.
#' @param bandwidth Density bandwidth parameter.
#'
#' @export
#' @return List of plots
mr_density_plot <- function(
  singlesnp_results,
  mr_results,
  exponentiate = FALSE,
  bandwidth = "nrd0"
) {
  dat_dt <- data.table::as.data.table(singlesnp_results)
  combos <- unique(dat_dt[, .(id.exposure, id.outcome)])

  res <- lapply(seq_len(nrow(combos)), function(i) {
    exp_id <- combos$id.exposure[i]
    out_id <- combos$id.outcome[i]
    d <- as.data.frame(dat_dt[id.exposure == exp_id & id.outcome == out_id])
    if (sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    d$SNP <- as.character(d$SNP)

    d2 <- subset(d, !grepl("All - ", SNP))
    d1 <- subset(mr_results, id.exposure == d2$id.exposure[1] & id.outcome == d2$id.outcome[1])

    xint <- 0
    if (exponentiate) {
      d$b <- exp(d$b)
      d$up <- exp(d$up)
      d$lo <- exp(d$lo)
      xint <- 1
    }

    ggplot2::ggplot(d2, ggplot2::aes(x = b)) +
      ggplot2::geom_vline(xintercept = xint, linetype = "dotted") +
      ggplot2::geom_density(ggplot2::aes(weight = 1 / se), bw = bandwidth) +
      ggplot2::geom_point(y = 0, colour = "red", ggplot2::aes(size = 1 / se)) +
      ggplot2::geom_vline(data = mr_results, ggplot2::aes(xintercept = b, colour = method)) +
      ggplot2::scale_colour_brewer(type = "qual") +
      ggplot2::labs(x = "Per SNP MR estimate")
  })
  res
}

#' Funnel plot
#'
#' Create funnel plot from single SNP analyses.
#'
#' @param singlesnp_results from [mr_singlesnp()].
#'
#' @export
#' @return List of plots
mr_funnel_plot <- function(singlesnp_results) {
  dat_dt <- data.table::as.data.table(singlesnp_results)
  combos <- unique(dat_dt[, .(id.exposure, id.outcome)])

  res <- lapply(seq_len(nrow(combos)), function(i) {
    exp_id <- combos$id.exposure[i]
    out_id <- combos$id.outcome[i]
    d <- as.data.frame(dat_dt[id.exposure == exp_id & id.outcome == out_id])
    if (sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    am <- grep("All", d$SNP, value = TRUE)
    d$SNP <- gsub("All - ", "", d$SNP)
    am <- gsub("All - ", "", am)
    ggplot2::ggplot(subset(d, !SNP %in% am), ggplot2::aes(y = 1 / se, x = b)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(
        data = subset(d, SNP %in% am),
        ggplot2::aes(xintercept = b, colour = SNP)
      ) +
      # ggplot2::scale_colour_brewer(type="qual") +
      ggplot2::scale_colour_manual(
        values = c(
          "#a6cee3",
          "#1f78b4",
          "#b2df8a",
          "#33a02c",
          "#fb9a99",
          "#e31a1c",
          "#fdbf6f",
          "#ff7f00",
          "#cab2d6",
          "#6a3d9a",
          "#ffff99",
          "#b15928"
        )
      ) +
      ggplot2::labs(y = expression(1 / SE[IV]), x = expression(beta[IV]), colour = "MR Method") +
      ggplot2::theme(legend.position = "top", legend.direction = "vertical")
  })
  res
}
