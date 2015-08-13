#source("/home/ph14916/MR/R/SNP_lookups_pipeline.R")
assign("last.warning", NULL, envir = baseenv())
rm(list = ls())
#setwd("/projects/MRC-IEU/users/ph14916/MR/telomeres/data") 
#library(biomaRt)
#setwd("/projects/MRC-IEU/users/ph14916/MR/vitE/") 
#setwd("/projects/MRC-IEU/users/ph14916/MR/selenium/")
#setwd("/projects/MRC-IEU/users/ph14916/MR/fetal_Hb") 
setwd("/projects/MRC-IEU/users/ph14916/MR/telomeres_GWASuniform") 

#source("/home/ph14916/MR/R/find_r2_allele_freq_1k_function.R")
source("/home/ph14916/MR/R/SNPlookups_GWASuniform_functions.R")

trait="telomeres"
instruments.file="data/telomere_instruments_GWASuniform.txt"
snp.exposure<-"SNPs"
effect="beta"
se="se"
eaf.trait="eaf.exclFamhs"
allele1="effect.allele" #effect allele
allele2="other.allele"
trait.col=NA #set to NA if trait.col is same as trait

###################################################################
#get position info from ENSEMBL if not present in instruments file#
###################################################################
snp.tab<-read.table(instruments.file,sep="\t",head=T,colClasses="character")
snp.tab<-snp.tab[,names(snp.tab)!="chr_name"]
snp.tab<-snp.tab[,names(snp.tab)!="chrom_start"]
snp.tab[,snp.exposure]<-gsub(" ","",snp.tab[,snp.exposure])
ensembl<-ensembl_get_position(snp=snp.tab[,snp.exposure])
snp.tab<-merge(snp.tab,ensembl,by.x=snp.exposure,by.y="refsnp_id")
snp.tab<-snp.tab[!is.na(snp.tab$chrom_start),]
if(is.na(trait.col)) { 
  snp.tab$exposure<-trait
  trait.col<-"exposure"
}
snp.tab[,trait.col]<-gsub("/","_",snp.tab[,trait.col])
write.table(snp.tab,instruments.file,sep="\t",col.names=T,row.names=F,quote=F)


#download r2 and allele frequency info from 1000 genomes
find_r2_1k_function(trait=trait,instruments.file=instruments.file,
	chr_id="chr_name",chr_pos="chrom_start",snp.exposure=snp.exposure,genome.build.hg="hg38",
	liftover.dir= "/projects/MRC-IEU/programs/twosampleMR/liftover/",
	tabix.dir="/projects/MRC-IEU/programs/twosampleMR/tabix-0.2.6/",
	vcftools.dir<-"/projects/MRC-IEU/programs/twosampleMR/vcftools_v0.1.13/",
	sample.1k="/projects/MRC-IEU/programs/twosampleMR/integrated_call_samples_v3.20130502.ALL.panel",
	super.pop="EUR")
	
#gwasplusSNPs[gwasplusSNPs$snp %in% c("rs10936599","rs10936601","rs12696304","rs1317082","rs4387287","rs9419958","rs9420907"),]

#############
#SNP lookups#
#############
databases.excl<-c("ImmuneCellScience","metabolome") # exclude warninthese databases to speed things up
outcome.database<-outcome.database[!outcome.database %in% databases.excl]
#do the SNP lookups in each database
outcome.database# to see a list of available databases

snp.lookups(outcome.database="metabolome",trait=trait,snp.exposure=snp.exposure,trait.col=trait.col,effect=effect,se=se,chr_id="chr_name",chr_pos="chrom_start",eaf.trait=eaf.trait,
	allele1=allele1,allele2=allele2,
	instruments.file=instruments.file,
	liftover.dir="/projects/MRC-IEU/programs/twosampleMR/liftover/",
	gwas.dir="/projects/MRC-IEU/publicdata/gwas_uniform/",build.plus.strand.same=F)

#database="metabolome"
#chr_id="chr_name"
#chr_pos="chrom_start"
#liftover.dir="/projects/MRC-IEU/programs/twosampleMR/liftover/"
#gwas.dir="/projects/MRC-IEU/publicdata/gwas_uniform/"
#build.plus.strand.same=F
	
warnings()
############################################
#collate all SNP lookups into a single file#
############################################
unlink("data/analysis1_allstudies.txt")
files.list<-dir("data")[grep("analysis1",dir("data"))]
all.results<-NULL
for(file.name in files.list){
	print(file.name)
	all.results[[file.name]]<-read.table(paste("data/",file.name,sep=""),sep="\t",colClasses="character",head=T,fill=T,quote="")
}
all.results<-do.call(rbind,all.results)
write.table(all.results,"data/analysis1_allstudies.txt",sep="\t",col.names=T,row.names=F,quote=F)

#######################################
#create a summary table of SNP lookups#
#######################################

all.results<-read.table("data/analysis1_allstudies.txt",sep="\t",head=T,colClasses="character",quote="")
res.tab<-NULL
for(database in unique(all.results$database)){
	print("database")
	print(database)
	results.subset.database<-all.results[all.results$database==database,]
	outcomes<-unique(results.subset.database$outcome)
	for(outcome in outcomes){
		print("outcome")
		print(outcome)
		results.subset.outcome<-results.subset.database[results.subset.database$outcome==outcome,]
		mean.N_case.outcome<-round(mean(as.numeric(	results.subset.outcome$N_case.outcome)),1)
		mean.N_control.outcome<-round(mean(as.numeric(	results.subset.outcome$N_control.outcome)),1)
		mean.N_total.outcome<-round(mean(as.numeric(results.subset.outcome$N_total.outcome)),1)
		n.snps<-length(	results.subset.outcome$snp.exposure)
		database<-unique(results.subset.outcome$database.outcome)
		strand<-unique(results.subset.outcome$strand.outcome)
		genome.build<-unique(results.subset.outcome$genome.build.outcome)
		data.tab<-data.frame(matrix(c(database,outcome,mean.N_case.outcome,mean.N_control.outcome,mean.N_total.outcome,n.snps,strand,genome.build),nrow=length(mean.N_case.outcome),ncol=8))
		names(data.tab)<-c("database","outcome","N_cases","N_controls","N_total","N_snps","strand","genome.build")
		res.tab[[outcome]]<-data.tab
	}
}
res.tab<-do.call(rbind,res.tab)

nice.names<-read.table("/projects/MRC-IEU/users/ph14916/MR/results_nice_names_ncases_ncontrols_v2.txt",sep="\t",colClasses="character",head=T,fill=T,quote="")
#res.nice<-merge(all.results,nice.names,by.x=c("database.outcome","outcome"),by.y=c("database","outcome"))
res.tab2<-merge(res.tab,nice.names,by=c("database","outcome"),all.x=T)
names(res.tab2)[names(res.tab2)=="N_cases"]<-"N_cases_from_database"
names(res.tab2)[names(res.tab2)=="N_controls"]<-"N_controls_from_database"
names(res.tab2)[names(res.tab2)=="ncases"]<-"N_cases_from_paper"
names(res.tab2)[names(res.tab2)=="ncontrols"]<-"N_controls_from_paper"
write.table(res.tab2,"data/summary_of_studies.txt",sep="\t",col.names=T,row.names=F,quote=F)
