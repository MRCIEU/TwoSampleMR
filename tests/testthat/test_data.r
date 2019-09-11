context("Associations")
library(TwoSampleMR)

a <- extract_instruments(c("IEU-a-300", "IEU-a-301"), force_server=TRUE)
b <- extract_outcome_data(a$SNP, c("IEU-a-2","IEU-a-7"))

test_that("extract_outcome_data is working", {
	expect_gt(nrow(b), 0)
})

ab <- harmonise_data(a,b)

test_that("harmonisation", {
	expect_equal(
		unique(paste(ab$id.exposure, ab$id.outcome)) %>% length,
		4
	)
})

