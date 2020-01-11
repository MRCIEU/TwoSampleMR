load("~/repo/mr_base_paper/results/lpa_ldl_trigs.RData")
disease <- subset(res, id.outcome %in% subset(outcomes, category == "Disease")$id)
risk_factor <- subset(res, id.outcome %in% subset(outcomes, category == "Risk factor")$id)


# Load example data

load(system.file("data/forest_plot_example.rdata", package="TwoSampleMR"))

# Create different forest plots

forest_plot(mr_result, in_columns=TRUE, by_category=TRUE, xlab="Effect of")
forest_plot(mr_result, in_columns=TRUE, by_category=TRUE, exponentiate=TRUE, trans="log2")

forest_plot(mr_result, in_columns=TRUE, by_category=FALSE)
forest_plot(mr_result, in_columns=TRUE, by_category=FALSE, exponentiate=TRUE, trans="log2")

forest_plot(mr_result, in_columns=FALSE, by_category=TRUE)
forest_plot(mr_result, in_columns=FALSE, by_category=FALSE)

