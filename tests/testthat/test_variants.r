context("Variants")
library(TwoSampleMR)

o1 <- variants_gene("ENSG00000123374")
o2 <- variants_gene("ENSG00000123374", 100000)

test_that("variants",
{
	expect_gt(nrow(o1), 0)
	expect_gt(nrow(o2), nrow(o1))
})


