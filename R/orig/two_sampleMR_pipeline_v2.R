###############
#two sample MR#
###############
# setwd("/projects/MRC-IEU/users/ph14916/MR/lipids/")
assign("last.warning", NULL, envir = baseenv())
rm(list = ls())

setwd("/projects/MRC-IEU/users/ph14916/MR/telomeres/")
trait<-"telomeres"
source("/home/ph14916/MR/R/twosampleMR_function_v7.R")

library(plyr)
library(meta)
library(biomaRt)

res.tab<-read.table("data/analysis1_allstudies.txt",sep="\t",head=T,colClasses="character",quote="")
#source("/home/ph14916/MR/R/two_sample_MR_quality_control_v2.R")
source("/home/ph14916/MR/R/two_sample_MR_quality_control.R")

dim(res.tab)
databases<-unique(res.tab$database)
databases<-databases[!databases %in% c("ibd_ichip","metabolome","ImmuneCellScience")]
databases
#databases="IBS_dbGAP"
#databases=  "alcohol_dependence_dbGAP"
#databases=  "ibd_ucolitis"

twosampleMR(trait=trait,res.tab=res.tab,databases=databases,genome.build="hg38",liftover.dir="/projects/MRC-IEU/programs/twosampleMR/liftover/",out.file="v1")

database=  ; outcome=    
exposure= M18467 
warnings()
exposure<- "telomeres"
outcome<-  "diastolic_blood_pressure"
database<-  "diastolic_blood_pressure"
trait<-trait
genome.build="hg38"
liftover.dir="/projects/MRC-IEU/programs/twosampleMR/liftover/"
out.file="GWASsig"

###################################################
#Collate the MR results and create some nice names#
###################################################

rm(list = ls())
library(plyr)
wk.dir<-"/projects/MRC-IEU/users/ph14916/MR/lipids/"
setwd(wk.dir) 

files.list<-dir("results/MRresults")
files.list<-files.list[grep("v1",files.list)]
#grep("ssgac",files.list)
#files.list[699:700]
#files.list<-files.list[4]
all.results<-NULL
for(file.name in files.list){
	print(file.name)
	all.results[[file.name]]<-read.table(paste("results/MRresults/",file.name,sep=""),sep="\t",colClasses="character",head=T)
}
all.results<-do.call(rbind.fill,all.results)
dim(all.results)
nice.names<-read.table("/projects/MRC-IEU/users/ph14916/MR/results_nice_names_ncases_ncontrols.txt",sep="\t",colClasses="character",head=T,fill=T,quote="")
res.nice<-merge(all.results,nice.names,by=c("database","outcome"),all.x=T)
write.table(res.nice,"results/all_results.txt",sep="\t",col.names=T,row.names=F,quote=F)

	
#create nice names, merge all.results with previous nice_names file and also metabolic traits info file, then open in excel and edit missing fields; then merge all.results with that file
#nice.names<-read.table("results/results_nice_names.txt",sep="\t",head=T,colClasses="character",quote = "")
#meta.tab<-read.table("/projects/MRC-IEU/users/ph14916/GWAS_summary_data/metabolome/metabolic_traits_metabolome.txt",sep="\t",head=T,colClasses="character",fill=T,quote="")
#meta.tab<-meta.tab[,c("Metabolite.ID","Metabolite","Pathway","Super.pathway")]
#all.results<-merge(all.results,nice.names,by=c("database","outcome"),all.x=T)
#all.results<-merge(all.results,meta.tab,by.x="outcome",by.y="Metabolite.ID",all.x=T)

#all.results$nice_names[!is.na(all.results$Metabolite)]<-all.results$Metabolite[!is.na(all.results$Metabolite)]
#all.results$group[!is.na(all.results$Metabolite)]<-"metabolite"
#all.results$subgroup[!is.na(all.results$Metabolite)]<-all.results$Super.pathway[!is.na(all.results$Metabolite)]
#all.results<-all.results[,names(all.results)!="Metabolite"]
#all.results<-all.results[,names(all.results)!= "Super.pathway"]

#head(all.results[all.results$group=="metabolite" & !is.na(all.results$group),c("outcome", "nice_names" , "group","subgroup","Pathway")])

#all.results$Pathway[!is.na(all.results$Pathway)]

#write.table(unique(all.results[,c("database","outcome","nice_names","group","subgroup","Pathway")]),"results/summary.txt",sep="\t",col.names=T,row.names=F,quote=F)


