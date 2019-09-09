context("API")
library(TwoSampleMR)
library(dplyr)
library(magrittr)
library(httr)

test_that("status", 
{
	stat <- api_status()
	expect_true(is.list(stat))
	expect_gte(length(stat), 2)
})



test_that("gwasinfo", 
{
	expect_true(
		nrow(api_query('gwasinfo/2',access_token=NULL)) == 1
	)
	expect_equal(
		nrow(api_query('gwasinfo', query=list(id=c(2,1001)))), 
		2
	)
	expect_gt(
		nrow(gwasinfo()),
		100
	)
})




test_that("variants",

o <- variants_gene("ENSG00000123374")


)
