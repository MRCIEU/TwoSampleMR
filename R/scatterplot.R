#' Create scatter plot with fitted lines showing the causal effect estimate for different MR estimators
#'
#' Create scatter plot with fitted lines showing the causal effect estimate for different MR estimators.
#'
#' @param mr_results Output from [mr()].
#' @param dat Output from [harmonise_data()].
#' @export
#' @return List of plots
mr_scatter_plot <- function(mr_results, dat) {
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  combos <- unique(dat[, c("id.exposure", "id.outcome")])
  mrres <- lapply(seq_len(nrow(combos)), function(i) {
    d <- dat[dat$id.exposure == combos$id.exposure[i] & dat$id.outcome == combos$id.outcome[i], ]
    if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(
      mr_results,
      id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1]
    )
    mrres$a <- 0
    if ("MR Egger" %in% mrres$method) {
      temp <- mr_egger_regression(
        d$beta.exposure,
        d$beta.outcome,
        d$se.exposure,
        d$se.outcome,
        default_parameters()
      )
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }

    if ("MR Egger (bootstrap)" %in% mrres$method) {
      temp <- mr_egger_regression_bootstrap(
        d$beta.exposure,
        d$beta.outcome,
        d$se.exposure,
        d$se.outcome,
        default_parameters()
      )
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }

    if ("MR GRIP" %in% mrres$method) {
      temp <- mr_grip(
        d$beta.exposure,
        d$beta.outcome,
        d$se.exposure,
        d$se.outcome,
        default_parameters()
      )
      mrres$a[mrres$method == "MR GRIP"] <- stats::lm(
        d$beta.outcome - temp$b * d$beta.exposure ~ 1,
        weights = 1 / d$se.outcome^2
      )$coef
    }

    if (utils::packageVersion("ggplot2") <= "3.5.2") {
      ggplot2::ggplot(
        data = d,
        ggplot2::aes(x = beta.exposure, y = beta.outcome)
      ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            ymin = beta.outcome - se.outcome,
            ymax = beta.outcome + se.outcome
          ),
          colour = "grey",
          width = 0
        ) +
        ggplot2::geom_errorbarh(
          ggplot2::aes(
            xmin = beta.exposure - se.exposure,
            xmax = beta.exposure + se.exposure
          ),
          colour = "grey",
          height = 0
        ) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(
          data = mrres,
          ggplot2::aes(intercept = a, slope = b, colour = method),
          show.legend = TRUE
        ) +
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
        ggplot2::labs(
          colour = "MR Estimate",
          x = paste("SNP effect on", d$exposure[1]),
          y = paste("SNP effect on", d$outcome[1])
        ) +
        ggplot2::theme(legend.position = "top", legend.direction = "vertical") +
        ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
    } else {
      ggplot2::ggplot(
        data = d,
        ggplot2::aes(x = beta.exposure, y = beta.outcome)
      ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            ymin = beta.outcome - se.outcome,
            ymax = beta.outcome + se.outcome
          ),
          colour = "grey",
          width = 0
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            xmin = beta.exposure - se.exposure,
            xmax = beta.exposure + se.exposure
          ),
          colour = "grey",
          width = 0,
          orientation = "y"
        ) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(
          data = mrres,
          ggplot2::aes(intercept = a, slope = b, colour = method),
          show.legend = TRUE
        ) +
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
        ggplot2::labs(
          colour = "MR Estimate",
          x = paste("SNP effect on", d$exposure[1]),
          y = paste("SNP effect on", d$outcome[1])
        ) +
        ggplot2::theme(legend.position = "top", legend.direction = "vertical") +
        ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
    }
  })
  names(mrres) <- paste(combos$id.exposure, combos$id.outcome, sep = ".")
  mrres
}


blank_plot <- function(message) {
  ggplot2::ggplot(data.frame(a = 0, b = 0, n = message)) +
    ggplot2::geom_text(ggplot2::aes(x = a, y = b, label = n)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}
