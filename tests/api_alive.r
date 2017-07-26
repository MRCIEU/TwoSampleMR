library(TwoSampleMR)

a <- extract_instruments(2)
b <- clump_data(a)
message(nrow(b))
c <- ld_matrix(a$SNP)
message(nrow(c))
d <- extract_outcome_data(a$SNP, 7)
e <- available_outcomes()
message(nrow(e))

