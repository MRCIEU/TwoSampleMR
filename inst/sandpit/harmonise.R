# Assumptions - effect allele is present in exposure and outcome
# Possibilities

# Both alleles in E and O and freq for both
# - 

# 2E + 2O + Ef + Of
# 2E + 2O + Ef
# 2E + 2O + Of
# 2E + 2O

# 1E + 2O + Ef + Of
# 1E + 2O + Ef
# 1E + 2O + Of
# 1E + 2O

# 2E + 1O + Ef + Of
# 2E + 1O + Ef
# 2E + 1O + Of
# 2E + 1O

# 1E + 1O + Ef + Of
# 1E + 1O + Ef
# 1E + 1O + Of
# 1E + 1O




a <- read.table("inst/extdata/alleles.txt", he=T, stringsAsFactors=FALSE)
SNP <- a$SNP
A1 <- a$A1
A2 <- a$A2
B1 <- a$B1
B2 <- a$B2
betaA <- a$betaA
betaB <- a$betaB
fA <- a$fA
fB <- a$fB


exposure_dat <- data.frame(
	SNP=a$SNP,
	effect_allele.exposure=a$A1,
	other_allele.exposure=a$A2,
	beta.exposure=a$betaA,
	eaf.exposure=a$fA,
	exposure="A",
	stringsAsFactors=FALSE
)


outcome_dat <- data.frame(
	SNP=a$SNP,
	effect_allele.outcome=a$B1,
	other_allele.outcome=a$B2,
	beta.outcome=a$betaB,
	eaf.outcome=NA,
	id.outcome="B",
	stringsAsFactors=FALSE
)

dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat)

a <- read.table("inst/data/alleles.txt", he=T, stringsAsFactors=FALSE)
SNP <- a$SNP
A1 <- a$A1
A2 <- a$A2
B1 <- a$B1
B2 <- a$B2
betaA <- a$betaA
betaB <- a$betaB
fA <- a$fA
fB <- a$fB

fA[4] <- NA
harmonise_21(SNP, A1, A2, B1, betaA, betaB, fA, rep(NA, length(SNP)), 0.08, 2)
harmonise_22(SNP, A1, A2, B1, B2, betaA, betaB, fA, fB, 0.08, 3)
