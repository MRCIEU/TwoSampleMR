context("singlesnp")

dat <- make_dat("ieu-a-2", "ieu-a-7")

test_that("singlesnp", {
	w <- mr_singlesnp(dat)
	expect_true(nrow(w) == sum(dat$mr_keep) + 2)
})


test_that("singlesnp_plot", {
	w <- mr_singlesnp(dat)
	p <- mr_forest_plot(w)
	expect_true(length(p) == 1)
})

test_that("density plot", {
	w <- mr_singlesnp(dat)
	m <- mr(dat, method_list="mr_ivw")
	p <- mr_density_plot(w, m)
	expect_true(length(p) == 1)
})

test_that("singlesnp_plot", {
	w <- mr_singlesnp(dat)
	p <- mr_funnel_plot(w)
	expect_true(length(p) == 1)
})
