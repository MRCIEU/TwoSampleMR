library(dplyr)

al <- matrix(c(
"A", "AT", "AT", "A",
"A", "AT", "A", "AT",
"A", "AT", "D", "I",
"A", "AT", "I", "D",
"AT", "A", "D", "I",
"AT", "A", "I", "D",
"D", "I", "AT", "A",
"D", "I", "A", "AT",
"I", "D", "AT", "A",
"I", "D", "A", "AT",
"I", "D", "AT", "AT",
"I", "D", "D", "I"
), nrow=4) %>% t

TwoSampleMR:::recode_indels_22(al[,1], al[,2], al[,3], al[,4]) %>% cbind(al, .)
TwoSampleMR:::recode_indels_21(al[,1], al[,2], al[,3]) %>% cbind(al, .)
TwoSampleMR:::recode_indels_12(al[,1], al[,3], al[,4]) %>% cbind(al, .)




a <- extract_instruments(2)
b <- extract_outcome_data(a$SNP, 7)
ab <- harmonise_data(a,b)
mr(ab, method_list="mr_ivw")


a1 <- a
b1 <- b

a1$effect_allele.exposure[1] <- "GAAAA"
b1$other_allele.outcome[23] <- "GAAAA"

ab1 <- harmonise_data(a1,b1)
mr(ab1, method_list="mr_ivw")



a1 <- a
b1 <- b

a1$effect_allele.exposure[1] <- "GAAAA"
b1$other_allele.outcome[23] <- "I"
b1$effect_allele.outcome[23] <- "D"

ab1 <- harmonise_data(a1,b1)
mr(ab1, method_list="mr_ivw")



a1 <- a
b1 <- b

b1$other_allele.outcome <- NA

ab1 <- harmonise_data(a1,b1)
mr(ab1, method_list="mr_ivw")



a1 <- a
b1 <- b

a1$effect_allele.exposure[1] <- "GAAAA"
b1$other_allele.outcome[23] <- "I"
b1$effect_allele.outcome[23] <- "D"
b1$other_allele.outcome <- NA

ab1 <- harmonise_data(a1,b1)
mr(ab1, method_list="mr_ivw")



a1 <- a
b1 <- b

a1$effect_allele.exposure[1] <- "GAAAA"
b1$other_allele.outcome[23] <- "I"
b1$effect_allele.outcome[23] <- "D"
b1$other_allele.outcome <- NA

ab1 <- harmonise_data(a1,b1)
mr(ab1, method_list="mr_ivw")



x <- data.frame(
	SNP=c("a","b", "c", "d", "e", "f"),
	beta=1,
	se=1,
	effect_allele=c("D", "ACTG", "A", "G", "C", "TD"),
	other_allele=c("I", "F", "T", "C", "G", "A"),
	stringsAsFactors=FALSE
)

format_data(x)


