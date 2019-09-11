context("Associations")
library(TwoSampleMR)

test_that("variants",
	o <- variants_gene("ENSG00000123374")
	expect_gt(nrow(o), 0)
)
