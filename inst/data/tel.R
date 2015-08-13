a <- read.table("~/Desktop/analysis1_allstudies.txt", he=T, sep="\t")
b <- subset(a, outcome %in% c("cardiogram", "bladdercancer"))
tel <- subset(b, select=c(snp.exposure, allele1, allele2, effect.exposure, se.exposure, chr_name, chr_pos, outcome, effect_allele, other_allele, effect.outcome, se.outcome, eaf.outcome, eaf.trait))
names(tel) <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "effect.exposure", "se.exposure", "chr_name", "chrom_start", "outcome", "effect_allele.outcome", "other_allele.outcome", "effect.outcome", "se.outcome", "eaf.outcome", "eaf.exposure")
tel$exposure <- "Telomere length"
save(tel, file="tel.RData")
