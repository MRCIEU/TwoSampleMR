skip()

exp_dat <- extract_instruments("ieu-a-2")
out_dat <- extract_outcome_data(exp_dat$SNP, "ieu-a-7")
dat <- make_dat("ieu-a-2", "ieu-a-7") %>% add_metadata()
dat2 <- make_dat()

save(exp_dat, out_dat, dat, dat2, file=file.path("inst", "extdata", "test_commondata.RData"))
