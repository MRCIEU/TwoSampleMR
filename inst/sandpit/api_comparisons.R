exposure_dat <- read_exposure_data("inst/data/telomere_length.txt")

options(mrbaseapi="http://api.mrbase.org/")
exp1 <- clump_data(exposure_dat)

options(mrbaseapi="http://scmv-webapps.epi.bris.ac.uk:5000/")
exp2 <- clump_data(exposure_dat)

all(exp1==exp2)


options(mrbaseapi="http://scmv-webapps.epi.bris.ac.uk:5000/")
inst1 <- extract_instruments(89, clump=TRUE)

options(mrbaseapi="http://api.mrbase.org/")
inst2 <- extract_instruments(89, clump=TRUE)

table(inst1==inst2)
for (i in seq_len(ncol(inst1))) {
	print(sum(inst1[,i] == inst2[,i]))
}

options(mrbaseapi="http://scmv-webapps.epi.bris.ac.uk:5000/")
ao1 <- available_outcomes()

options(mrbaseapi="http://api.mrbase.org/")
ao2 <- available_outcomes()

table(ao1 == ao2)


options(mrbaseapi="http://scmv-webapps.epi.bris.ac.uk:5000/")
out1 <- extract_outcome_data(inst1$SNP, c(1,2), proxies=TRUE)

options(mrbaseapi="http://api.mrbase.org/")
out2 <- extract_outcome_data(c(inst2$SNP, exposure_dat$SNP), c(1,2), proxies=TRUE)

table(out1 == out2)



dat1 <- harmonise_data(exposure_dat, out2)
dat2 <- harmonise_data(inst2, out2)
