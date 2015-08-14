a <- read.table("~/Desktop/analysis1_allstudies.txt", he=T, sep="\t")
b <- subset(a, outcome %in% c("cardiogram", "bladdercancer"))
tel <- subset(b, select=c(snp.exposure, allele1, allele2, effect.exposure, se.exposure, chr_name, chr_pos, outcome, effect_allele, other_allele, beta.outcome, se.outcome, eaf.outcome, eaf.trait))
names(tel) <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "beta.exposure", "se.exposure", "chr_name", "chrom_start", "outcome", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", "se.outcome", "eaf.outcome", "eaf.exposure")
tel$exposure <- "Telomere length"



card <- subset(tel, outcome=="cardiogram", select=c(SNP, effect_allele.outcome, other_allele.outcome, beta.outcome, se.outcome, eaf.outcome))
names(card) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "eaf")


telo <- subset(tel, outcome=="cardiogram", select=c(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, eaf.exposure))
names(telo) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "eaf")


blad <- subset(tel, outcome=="bladdercancer", select=c(SNP, effect_allele.outcome, other_allele.outcome, beta.outcome, se.outcome, eaf.outcome))
names(blad) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "eaf")


write.table(card, file="inst/data/cardiogram.txt", row=F, col=T, qu=F)
write.table(blad, file="inst/data/bladdercancer.txt", row=F, col=T, qu=F)
write.table(telo, file="inst/data/telomere_length.txt", row=F, col=T, qu=F)

save(tel, file="inst/data/tel.RData")
