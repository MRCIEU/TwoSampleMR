test_that("Create local data", {
  skip_if(Sys.getenv('TWOSAMPLEMR_ENABLE_OPENGWAS_TESTS') != TRUE, "Tests requiring OpenGWAS")
skip("Skip local data creation unless you have good access to the API.")

exp_dat <- extract_instruments("ieu-a-2")
out_dat <- extract_outcome_data(exp_dat$SNP, "ieu-a-7")
dat <- make_dat("ieu-a-2", "ieu-a-7") %>% add_metadata()
dat2 <- make_dat()

save(exp_dat, out_dat, dat, dat2, file=file.path("inst", "extdata", "test_commondata.RData"), compress = "xz")
})
