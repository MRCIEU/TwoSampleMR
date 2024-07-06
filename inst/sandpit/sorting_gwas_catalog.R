
# In EXCEL substitute missing cells for NA

gwascat.file<-paste("~/gwascatalog_",Sys.Date(),".txt",sep="") #for newest version
download.file("https://www.ebi.ac.uk/gwas/api/search/downloads/alternative", gwascat.file,  method="curl",quiet = FALSE,cacheOK = FALSE)

# a <- read.table("~/Downloads/gwas_catalog_v1.0-downloaded_2015-09-21_2.txt", he=T, sep="\t", quote='"', comment="", stringsAsFactors=FALSE)

a<-read.table(gwascat.file,header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE,stringsAsFactors=F)
# a<-read.table("~/Downloads/gwascatalog.txt",header=TRUE,sep='\t',quote="",comment.char="",as.is=TRUE,stringsAsFactors=F)
# a<-read.table(gwascat.file,he=T, sep="\t", quote='"', comment="", stringsAsFactors=FALSE)
# a<-read.table("~/Downloads/gwascatalog.txt",he=T, sep="\t", quote='"', comment="", stringsAsFactors=FALSE)

b <- subset(a, select=c(DISEASE.TRAIT, PUBMEDID, FIRST.AUTHOR, DATE, SNPS,STRONGEST.SNP.RISK.ALLELE, REGION, REPORTED.GENE.S., OR.or.BETA, X95..CI..TEXT., P.VALUE, P.VALUE..TEXT., RISK.ALLELE.FREQUENCY, INITIAL.SAMPLE.DESCRIPTION,REPLICATION.SAMPLE.DESCRIPTION,MAPPED_TRAIT,MAPPED_TRAIT_URI,DATE.ADDED.TO.CATALOG))

# Make risk allele numeric
b$RISK.ALLELE.FREQUENCY <- as.numeric(b$RISK.ALLELE.FREQUENCY)

#exclude SNPs with missing rsids
b<-b[grep("rs",b$SNPS,ignore.case=TRUE,),] #exclude SNPs without an rsid 
b[grep("-",b$STRONGEST.SNP.RISK.ALLELE,ignore.case=TRUE,invert=TRUE),]

# Get effect allele
pos.allele<-gregexpr("[ATGC]",b$STRONGEST.SNP.RISK.ALLELE)
b$effect_allele<-lapply(seq_along(b$STRONGEST.SNP.RISK.ALLELE),FUN=function(x) substr(b$STRONGEST.SNP.RISK.ALLELE[x],unlist(pos.allele[x]),unlist(pos.allele[x])))

# Get year
b$DATE <- as.Date(b$DATE, format="%d-%b-%y")
b$year <- format(b$DATE, "%Y")

# Try to get the units 
Start<-unlist(lapply(b$X95..CI..TEXT.,FUN=function(x) unlist(gregexpr("] ",x))+2))
Stop<-nchar(b$X95..CI..TEXT.)
b$units<-unlist(lapply(seq_along(Start) ,FUN=function(x) substr(b$X95..CI..TEXT.[x],Start[x],Stop[x])))
b$units[which(unlist(regexpr("[:A-Za-z:]",b$units))==-1)]<-NA
b$units[which(b$units=="[NR]")]<-NA
b$units[which(b$units=="NR")]<-NA
b$units[which(b$units=="")]<-NA
b$units<-gsub("\\(","",b$units)
b$units<-gsub(")","",b$units)
b$units<-gsub("\\[","",b$units)
b$units<-gsub("]","",b$units)
b$units<-gsub("s\\.d\\.","SD",b$units) #without escapes . is interpretated as any character

b$type<-NA
b$type[unlist(lapply(c("older","higher","taller","increase","better","more", # higher
		"younger","lower","shorter","decrease","dcrease","decrea","fewer","worse", #lower
		"SD","unit","kg/m2","cm","msec","variance explained","% variance"), #other
		FUN=function(x) grep(x,b$units,ignore.case=TRUE)))]<-"continuous"

# Assume that anything with OR.or.BETA < 0.5 is not an odds ratio / is a continuous phenotype
b$type[which(b$OR.or.BETA<0.5 & is.na(b$units))]<-"continuous?"
# Guess that anything with OR.or.BETA > 2.0 is not an odds ratio / is a continuous phenotype
b$type[which(b$OR.or.BETA>2.0 & is.na(b$units))]<-"continuous?"
# Guess that OR.or.BETA 1.0 to 2.0 is an odds ratio / is a binary phenotype
b$type[which((b$OR.or.BETA>1.0 | b$OR.or.BETA<2.0) & is.na(b$units))] <-"binary?"

b$direction<-NA
b$direction[unlist(lapply(c("older","higher","taller","increase","better","more"), # higher
		FUN=function(x) grep(x,b$units,ignore.case=TRUE)))] <- "higher"
b$direction[unlist(lapply(c("younger","lower","shorter","decrease","dcrease","decrea","fewer","worse"), #lower
		FUN=function(x) grep(x,b$units,ignore.case=TRUE)))] <- "lower"

# Try to get standard errors and units from confidence intervals
# calculate two sets of standard errors, assuming OR.or.BETA is and isn't an odds ratio

#experimental new script#
#########################
pos<-regexpr("-",c("1.174-1,457" , "1.46,2.33", "0.49,0.098", "1.28,2.202","1.1-1.2"))
test1<-substr(c("1.174-1,457" , "1.46,2.33", "0.49,0.098", "1.28,2.202","1.1-1.2"),1,pos-1)
nums<-gregexpr("[:0-9:]",c("1.174-1,457" , "1.46,2.33", "0.49,0.098", "1.28,2.202","1.1-1.2"))
end<-lapply(seq_along(nums),FUN=function(x) unlist(nums[x])[length(unlist(nums[x]))]) # this finds the position of the last number in the sequence
test2<-substr(c("1.174-1,457" , "1.46,2.33", "0.49,0.098", "1.28,2.202","1.1-1.2"),pos+1,end)

test1
test2
as.numeric(test2)


c<-b
b<-b[grep(",",b$ci95,invert=T),]


pos.start<-regexpr("\\[",b$X95..CI..TEXT.)+1
pos.end<-regexpr("]",b$X95..CI..TEXT.)-1
ci95<-substr(b$X95..CI..TEXT.,pos.start,pos.end)
ci95
ci95[which(ci95=="")]<-NA
pos<-regexpr("-",ci95)
ci95[!is.na(ci95)][grep("-",ci95[!is.na(ci95)],invert=TRUE)]
nums<-gregexpr("[:0-9:]",ci95)
end<-lapply(seq_along(nums),FUN=function(x) unlist(nums[x])[length(unlist(nums[x]))]) # this finds the position of the last number in the sequence
num1<-substr(ci95,1,pos-1)
num2<-substr(ci95,pos+1,end)
unique(num2)

b[grep("mg/dl decrease",b$X95..CI..TEXT.),]
pos2<-regexpr(",",ci95)
pos[pos2!=-1]<-pos2[pos2!=-1] #for the scenario where instead of '-' there is a ','

se<-(as.numeric(num2)-as.numeric(num1))/(1.96*2)
lnse<-(log(as.numeric(num2))-log(as.numeric(num1)))/(1.96*2)
snps.table$se<-se
snps.table$lnse<-lnse
snps.table$ci95<-ci95



###
b$X95..CI..TEXT.[grep("NR", b$X95..CI..TEXT.)]<- NA
index1 <- grep("\\] \\(", b$X95..CI..TEXT.)
index1l <- grepl("\\] \\(", b$X95..CI..TEXT.)
temp1 <- do.call(rbind, strsplit(b$X95..CI..TEXT., split="] \\("))
temp1[,1] <- gsub("\\(|\\)|\\[|\\]| ", "", temp1[,1])
temp1[,1] <- gsub("^-", "@", temp1[,1])
temp1[,1] <- gsub("--", "-@", temp1[,1])
temp1 <- cbind(temp1, do.call(rbind, strsplit(temp1[,1], split="-"))[,1:2])
temp1 <- data.frame(temp1, stringsAsFactors=FALSE)
temp1[,3] <- gsub("^\\.", "0.", temp1[,3])
temp1[,4] <- gsub("^\\.", "0.", temp1[,4])
temp1[,3] <- gsub("^@\\.", "@0.", temp1[,3])
temp1[,4] <- gsub("^@\\.", "@0.", temp1[,4])
temp1[,3] <- gsub("@", "-", temp1[,3])
temp1[,4] <- gsub("@", "-", temp1[,4])
temp1[,3] <- as.numeric(temp1[,3])
temp1[,4] <- as.numeric(temp1[,4])
temp1[,2] <- gsub("\\)", "", temp1[,2])
temp1[,2] <- gsub("\\(", "", temp1[,2])
temp1[,2] <- gsub("^\\[", "", temp1[,2])
temp1[,2] <- gsub("\\]$", "", temp1[,2])
temp1[is.na(temp1[,4]),3] <- NA
temp1[is.na(temp1[,3]),4] <- NA
i1 <- which(temp1[,3] > temp1[,4])
x <- temp1[i1,3]
temp1[i1,3] <- temp1[i1,4]
temp1[i1,4] <- x


b$UNITS <- temp1$X2
b$UNITS[!grepl("[[:alpha:]]", b$UNITS) & !is.na(b$UNITS) & b$OR.or.BETA > 0.5 & !index1l] <- "Odds ratio"
b$UNITS[grepl("^\\[", b$UNITS) & grepl("\\]$", b$UNITS) & !is.na(b$UNITS) & b$OR.or.BETA > 0.5 & !index1l] <- "Odds ratio"
b$UNITS[grepl("^[0-9]", b$UNITS) & grepl("[0-9]$", b$UNITS) & !is.na(b$UNITS) & b$OR.or.BETA > 0.5 & !index1l] <- "Odds ratio"
b$UNITS[grep("[[:alpha:]]", b$UNITS, invert=TRUE)] <- NA
b$ci_lower <- temp1$X3
b$ci_upper <- temp1$X4


## modify decrease / increase units so that they are all the same
# older younger
# decrease increase
# taller shorter
# higher lower
# convert effect sizes and OR appropriately

b$OR.or.BETA.fixed <- b$OR.or.BETA

i2 <- grepl("decrease", b$UNITS)
b$OR.or.BETA.fixed[i2] <- b$OR.or.BETA[i2] * -1

i2 <- grepl("shorter", b$UNITS)
b$OR.or.BETA.fixed[i2] <- b$OR.or.BETA[i2] * -1

i2 <- grepl("lower", b$UNITS)
b$OR.or.BETA.fixed[i2] <- b$OR.or.BETA[i2] * -1

i2 <- grepl("younger", b$UNITS)
b$OR.or.BETA.fixed[i2] <- b$OR.or.BETA[i2] * -1

b$OR.or.BETA.fixed[which(b$UNITS=="Odds ratio")] <- log(b$OR.or.BETA[which(b$UNITS=="Odds ratio")])

b$SE <- abs(b$ci_upper - b$ci_lower) / 3.92
b$SE[which(b$UNITS=="Odds ratio")] <- abs(log(b$ci_upper[which(b$UNITS=="Odds ratio")]) - log(b$ci_lower[which(b$UNITS=="Odds ratio")])) / 3.92

b$UNITS <- gsub("decrease", "increase", b$UNITS)
b$UNITS <- gsub("lower", "higher", b$UNITS)
b$UNITS <- gsub("shorter", "taller", b$UNITS)
b$UNITS <- gsub("younger", "older", b$UNITS)


b$REPORTED.GENE.S.2 <- sapply(strsplit(b$REPORTED.GENE.S., split=","), function(x) x[1])

b1 <- subset(b, select=c(DISEASE.TRAIT, P.VALUE..TEXT., PUBMEDID, year, SNP, REGION, REPORTED.GENE.S.2, effect_allele, OR.or.BETA.fixed, SE, P.VALUE, UNITS, RISK.ALLELE.FREQUENCY))

names(b1) <- c("Phenotype", "Phenotype info", "PubmedID", "Year", "SNP", "Region", "Gene", "Allele", "Effect", "SE", "P-value", "Units", "eaf")

b1$Phenotype <- as.factor(b1$Phenotype)
b1$Year <- as.factor(b1$Year)
i1 <- which(b1["P-value"] > 1)
b1[i1,"P-value"] <- as.numeric(gsub("\\+", "-", as.character(b1[i1,"P-value"])))

library(biomaRt)
b1$Allele[b1$Allele == "?"] <- NA
i3 <- !is.na(b1$Allele)

Mart <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
Attr<-listAttributes(Mart)
ensembl<-getBM(attributes=c("refsnp_id","chr_name","chrom_start","allele", "minor_allele", "minor_allele_freq"),filters="snp_filter",values=b1$SNP[i3],mart=Mart)
ensembl <- subset(ensembl, !duplicated(refsnp_id))
temp <- subset(b1, select=c(SNP, Allele))
temp$index <- seq_len(nrow(temp))
temp <- merge(temp, ensembl, by.x="SNP", by.y="refsnp_id", all.x=TRUE)
temp <- temp[order(temp$index),]
alleles <- data.frame(t(sapply(strsplit(temp$allele, split="/"), function(x) x[1:2])), stringsAsFactors=FALSE)
alleles$effect_allele <- temp$Allele

i4 <- sapply(seq_len(nrow(alleles)), function(i) alleles[i, which(alleles[i, 1:2] != alleles[i, 3])[1]])
i4 <- sapply(i4, function(x) if(is.null(x)) NA else x)
b1$other_allele <- i4
b1$eaf[b1$eaf >= 1 | b1$eaf <= 0] <- NA

eafindex1 <- which(temp$Allele == temp$minor_allele & is.na(b1$eaf))
eafindex2 <- which(temp$Allele != temp$minor_allele & is.na(b1$eaf))

b1$eaf[eafindex1] <- temp$minor_allele_freq[eafindex1]
b1$eaf[eafindex2] <- 1 - temp$minor_allele_freq[eafindex2]


gwas_catalog <- b1
save(gwas_catalog, file="gwas_catalog.RData")
