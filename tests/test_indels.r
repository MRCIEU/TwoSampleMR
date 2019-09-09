
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

test_that("recode indels", {

	expect_true(all(
		recode_indels_22(al[,1], al[,2], al[,3], al[,4])$keep == 
		c(T,T,T,T,T,T,T,T,T,T,F,T)
	))

	expect_true(all(
		recode_indels_21(al[,1], al[,2], al[,3])$keep == 
		c(T,T,T,T,T,T,F,F,F,F,F,F)
	))

	expect_true(all(
		recode_indels_12(al[,1], al[,3], al[,4])$keep == 
		c(T,T,F,F,F,F,T,T,T,T,F,F)
	))

})




a <- extract_instruments(2)
b <- extract_outcome_data(a$SNP, 7)
ab <- harmonise_data(a,b)
b <- mr(ab, method_list="mr_ivw")$b


test_that("harmonise indels1", {
	a1 <- a
	b1 <- b

	a1$effect_allele.exposure[1] <- "GAAAA"
	b1$other_allele.outcome[23] <- "GAAAA"
	ab1 <- harmonise_data(a1,b1)
	expect_equal(mr(ab1, method_list="mr_ivw")$b, b)

	a1 <- a
	b1 <- b

	a1$effect_allele.exposure[1] <- "GAAAA"
	b1$other_allele.outcome[23] <- "I"
	b1$effect_allele.outcome[23] <- "D"

	ab2 <- harmonise_data(a1,b1)
	expect_equal(mr(ab2, method_list="mr_ivw")$b, b)
})









test_that("harmonise indels1", {
	a1 <- a
	b1 <- b

	b1$other_allele.outcome <- NA

	ab3 <- harmonise_data(a1,b1)
mr(ab3, method_list="mr_ivw")



a1 <- a
b1 <- b

a1$effect_allele.exposure[1] <- "GAAAA"
b1$other_allele.outcome[23] <- "I"
b1$effect_allele.outcome[23] <- "D"
b1$other_allele.outcome <- NA

ab4 <- harmonise_data(a1,b1)
mr(ab4, method_list="mr_ivw")




x <- data.frame(
	SNP=c("a","b", "c", "d", "e", "f"),
	beta=1,
	se=1,
	effect_allele=c("D", "ACTG", "A", "G", "C", "TD"),
	other_allele=c("I", "F", "T", "C", "G", "A"),
	stringsAsFactors=FALSE
)

format_data(x)


