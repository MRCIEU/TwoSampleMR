#source("/home/ph14916/MR/R/snplookups_function_v6.R")
#library(plyr)
##rm(list = ls())
#setwd("/home/ph14916/MR/results") 

#
snp.lookups <- function(outcome.database,eaf.trait,allele1,allele2,liftover.dir,trait,
	snp.exposure,trait.col,effect,se,chr_id,chr_pos,instruments.file,gwas.dir,build.plus.strand.same) { #database corresponding to disease GWAS, e.g. diagram, diagramMetabochip, cardiogram, cardiogramplusc4d,ilcco
	for(database in outcome.database){
		if(!"data" %in% dir()) dir.create("data")
		if(!"log" %in% dir()) dir.create("log")
		in.file<-paste("",instruments.file,sep="")
		trait.assoc<-read.table(in.file,sep="\t",head=T,colClasses="character",quote="",fill=T)
		pos<-unlist(lapply(c(snp.exposure,eaf.trait, allele1, allele2,effect,se,chr_id,chr_pos,trait.col),FUN=function(x) which(names(trait.assoc)==x)))
		trait.assoc<-trait.assoc[,pos]
		names(trait.assoc)<-c("snp.exposure","eaf.trait","allele1","allele2","effect.exposure","se.exposure","chr_name","chr_pos","exposure")
		write.table(unique(trait.assoc$snp.exposure),"data/snplist.txt",col.names=F,row.names=F,quote=F)
		trait.assoc<-trait.assoc[!duplicated(trait.assoc[,c("snp.exposure","exposure")]),]
		print(database)
		pos<-which(names.database==database)
		gwas.file<-paste(gwas.dir,gwas.files[pos],sep="") #file path for GWAS database
			
		#SNP lookups for ILCCO lung cancer database, ilcco summary associations are distributed across multiple files corresponding to different chromosomes and SNP id is in chr#:123 format, not rsid format
		if(sum(grep("ilcco",database))==1){ #ilcco SNP ids are in bed format not rsid format, coordinates are in build 37 / hg19 format
			#use lift over to convert SNP positions to hg19/GRCh37
			#note that gwas catalogue SNPs chr_name & chrom_start are hg38/CRCh38 format from ENSEMBL
			#SNPs in metabGWAS Shin nature genetics 2014 are build 36
			#bed format looks like this: chr1    743267  743268  rs3115860
			chr<-as.numeric(trait.assoc$chr_name)
			pos<-as.numeric(trait.assoc$chr_pos)
			snps<-trait.assoc$snp
			bed.matrix<-data.frame(matrix(c(paste("chr",chr,sep=""),pos,pos+1,snps),nrow=length(snps),ncol=4))
			write.table(bed.matrix,"data/bed.txt",sep=" ",col.names=F,row.names=F,quote=F)
			lift.CMD<-paste(liftover.dir,"liftOver data/bed.txt ",liftover.dir,"hg38ToHg19.over.chain.gz data/output.bed data/unlifted.bed",sep="")
			system(lift.CMD)
			Bed<-read.table("data/output.bed",sep="\t",head=F,stringsAsFactors=F)
			names(Bed)[names(Bed)=="V2"]<-"Chr_pos_hg19"
			names(Bed)[names(Bed)=="V1"]<-"chr_name.bed"
			names(Bed)[names(Bed)=="V4"]<-"snp"
			snps<-paste(Bed$chr_name.bed,":",Bed$Chr_pos_hg19,sep="")
			snps<-gsub("chr","",snps)
			Bed$bed<-snps
			write.table(unique(snps),"data/snplist.txt",col.names=F,row.names=F,quote=F)			
		}
	
		fgrep.cmd<-paste("fgrep -hwf data/snplist.txt ",
			gwas.file,
			" > data/gwasplusSNPs.txt",
			sep="") #-h supresses file prefix
			
		mulifile.database<-c("metabolome","GTEx","ImmuneCellScience")
		if(database %in% mulifile.database ){
			fgrep.cmd<-paste("fgrep -wf data/snplist.txt ",
			gwas.file,
			" > data/gwasplusSNPs.txt",
			sep="") #don't supress file prefix
		}
	
		system(fgrep.cmd)
		
		head.cmd<-paste("head -1 ",
			gwas.file,
			" > data/filehead.txt",
		sep="")	
		
		
		if(sum(grep("icbp",database))==1){
			head.cmd<-paste("sed -n 22,22p ", 
				gwas.file,
				" > data/filehead.txt",
			sep="")	
		}
		
		if(sum(grep("ilcco",database))==1){
			head.cmd<-paste("head -1 ",gwas.dir,"ILCCO/ALL/TRICL_ICR_MDACC_IARC_NIH_GWAMA_overall_chr1_dartmouth.out > data/filehead.txt",sep="") 
		}		
		
		if(sum(grep("metabolome",database))==1){
			head.cmd<-paste("head -1 ",gwas.dir,"metabolome/metabolites_meta/M00053.metal.pos.txt > data/filehead.txt",sep="") 
		}		
		
#		if(sum(grep("GTEx",database))==1){
#			head.cmd<-paste("head -1 ",gwas.dir,"GTEx/sigonly/Adipose_Subcutaneous.portal.eqtl  > data/filehead.txt",sep="") 
#		}	
		
		if(sum(grep("GTEx",database))==1){
			head.cmd<-paste("head -1 ",gwas.dir,"GTEx/alldata/Cells_Transformed_fibroblasts.cis.eqtl  > data/filehead.txt",sep="") 
		}	
		
	
		if(sum(grep("ImmuneCellScience",database))==1){
			head.cmd<-paste("head -1 ",gwas.dir,"ImmuneCellScience/2-GWASResults/GWA_Lin_1.txt  > data/filehead.txt",sep="") 
		}		
		
		if(sum(grep("telomere_outcomes",database))==1){
			head.cmd<-paste("head -1 ",gwas.dir,"telomere_outcomes/harmonized/AMD.txt  > data/filehead.txt",sep="") 
		}			
				
			
		system(head.cmd)
		
		app.cmd<-"cat data/filehead.txt  data/gwasplusSNPs.txt > data/temp ; mv data/temp  data/gwasplusSNPs.txt"
		system(app.cmd)
		
		gwasplusSNPs<-read.table("data/gwasplusSNPs.txt",head=T,colClasses="character",fill=T)
		if(sum(grep("icbp",database))==1 | sum(grep("dbGAP",database))==1 ){ 
			gwasplusSNPs<-read.table("data/gwasplusSNPs.txt",sep="\t",head=T,colClasses="character")
		}
		print(nrow(gwasplusSNPs))

		#nrows in gwasplusSNPs
		if(nrow(gwasplusSNPs)==0){		
			write.table(paste("the selected trait SNPs are not present in ",database," GWAS",sep=""),
				paste("log/",trait,"_",database,"_SNPlookups_log.txt",sep=""),col.names=F,row.names=F,quote=F)
				cat("the selected exposure SNPs are not present in the",database," database \n")
			next
		}
				                                     
		#generic script for renaming columns.  
		hwe.ids<-c("pHWE","pHWE..control.","hwe_p","HWE.P.value..All.sample..N.698.","Pexact..HWE.p.value.", "ControlHWE_P","control_hwe","unaff_P_HWE","hweP","HWE.p.values","hwe_min","HWE..Pearsosn..ctrls","hweP","HWE..ctrls.","HWE.P","pexhwe")
		hetp.ids<-c("HetPVal","HetPVa","Heterogeneity.study","q_p.value","het_pvalue","Q_Pval")
		r2.info.id<-c("RSQR_imp","rsq_imp","Rsq_imp","info_score" ,"info","Info","INFO" ,"BCAC.r2.imp","r2_info_median","OCAC.r2.imp","r2_info",
		"Imputation.Info.score..Impute2.","ImputationRsq","impute_info", "impute.v2.info.score","Phenotype_frequentist_add_Gender_C1_C2_C3_C4_C5_C6_C7_C8_C9_C10_C11_C12_C13_C14_C15_C16_C17_C18_C19_C20_score_info",
		"frequentist_add_info","qual")
		r2.proxy.id<-c("proxyR2")
		call.rate.ids<-c("Call.rate..control.","call.rate..95..confidence.","CallRate_con" ,"Call.rate","Call.Rate","Call.rate.con","CallRate_con","call")
		ncase.ids<-c("N_cases","ncases","number.of.cases","N_case","N_CASES","Ncases")
		ncon.ids<-c("N_controls","ncontrols","number.of.controls","N_cont","N_CONTROLS","Ncontrols","N_control","ncontrol","ncont")
		ntotal.ids<-c("N","TotalSampleSize","Sample_Size","sample.size","Sample.size","NMISS","n") #NMISS number of non missing genotypes in PLINK output
		p.ids<-c("pval_add","Pval","P","PVAL","P_Value","P.value","p","GC.Pvalue","P.value","pvalue","P.value","pval","SCAN.P","meta_2tP","Pvalue","p_sanger","p.Value","dbGC_P","p")
		chr.ids<-c("Chr","CHR","Chr.ID","Chr.ID","Chromosome","chromosome","chr_name","chr","X.CHROM","hg19chrc")
		bp.ids<-c("BP","POS","Chr.Position","Pos","Chr.Position","bp","Position","position", "chrom_start","pos")
		strand.find.ids<-c("forward","FORWARD","Forward","+","-","reverse","REVERSE","Reverse","top","bottom","strand")
		genome.build.find.ids<-c("Hg19","Hg18","Hg38","Hg39","hg19","hg18","hg38","hg39","b37","b36","b.36","b.37","b38","b39","hg18chr","hg19chrc","build")
		nstudies.ids<-c("n_studies")
		direction.ids<-c("direction","Direction","SCAN")
		lambda.ids<-c("lambda.estimate")
		pgc.ids<-c("pgc")
	
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs)) %in% tolower(hwe.ids)]<-"P.hwe.outcome"
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(hetp.ids)]<-"het.p.value.outcome"
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(r2.info.id)]<-"r2.info.imp.outcome"
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(call.rate.ids)]<-"genotype.call.rate.outcome"
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(ncase.ids)]<-"N_case.outcome"
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(ncon.ids)]<-"N_control.outcome"	
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs)) %in% tolower(ntotal.ids)]<-"N_total.outcome"	
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(nstudies.ids)]<-"N_studies.outcome" 
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(chr.ids)]<-"Chr_id.outcome" 
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(bp.ids)]<-"Chr_pos.outcome" 
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(p.ids)]<-"p.value.outcome" 	
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(strand.find.ids)]<-"possible_strand.outcome" 	
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(genome.build.find.ids)]<-"possible_build.outcome" 	
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(direction.ids)]<-"direction.outcome" 	
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(lambda.ids)]<-"lambda.outcome" 	
		names(gwasplusSNPs)[tolower(names(gwasplusSNPs))  %in% tolower(pgc.ids)]<-"pgc.outcome" 	
	
	
		if(database=="diagramplusmetabochip"){ #alleles are aligned to the forward strand of NCBI build 36; note that in metabolome GWAS (shin, Nature Genetics 2014), alleles were also aligned to the forward strange of build 36
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNPID"]<-"snp.outcome" #SNP rsid in diesase GWAS
			gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$OR))
			gwasplusSNPs$se.outcome<-(log(as.numeric(gwasplusSNPs$OR_95U))-log(as.numeric(gwasplusSNPs$OR_95L)))/(1.96*2)
			names(gwasplusSNPs)[names(gwasplusSNPs)=="EFFECT_ALLELE"]<-"effect_allele" #effect allele in disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="OTHER_ALLELE"]<-"other_allele" #effect allele in disease GWAS
			gwasplusSNPs$eaf.outcome<-NA 
			gwasplusSNPs$strand.outcome <- "forward"  
			gwasplusSNPs$genome.build.outcome<-"hg18"  #build 36 is hg18					
		}
		
		if(database=="diagram"){ #alleles are aligned to the forward strand of NCBI build 36; note that in metabolome GWAS (shin, Nature Genetics 2014), alleles were also aligned to the forward strange of build 36
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNP"]<-"snp.outcome" #SNP rsid in diesase GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="RISK_ALLELE"]<-"effect_allele" #effect allele in disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="OTHER_ALLELE"]<-"other_allele" #effect allele in disease GWAS
			gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$OR))
			gwasplusSNPs$se.outcome<-(log(as.numeric(gwasplusSNPs$OR_95U))-log(as.numeric(gwasplusSNPs$OR_95L)))/(1.96*2)
			gwasplusSNPs$eaf.outcome<-NA 
			gwasplusSNPs$strand.outcome<- "forward"  
			gwasplusSNPs$genome.build.outcome<-"hg18"  #build 36 is hg18
		}
		if(database=="cardiogramplusc4d"){
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNP_ID"]<-"snp.outcome" #SNP rsid in diesase GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="log_odds"]<-"effect.outcome" #log odds ratio for disease in GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="log_odds_se"]<-"se.outcome" #log odds ratio for disease in GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="reference_allele"]<-"effect_allele" #effect allele in disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="other_allele"]<-"other_allele" #effect allele in disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="eaf"]<-"eaf.outcome" #eaf in disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="ref_allele_frequency"]<-"eaf.outcome" #eaf in disease GWAS
		}
		
		if(database=="cardiogram"){
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNP"]<-"snp.outcome" #SNP rsid in diesase GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="log_odds"]<-"effect.outcome" #log odds ratio for disease in GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="log_odds_se"]<-"se.outcome" #log odds ratio for disease in GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="reference_allele"]<-"effect_allele" #effect allele in disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="other_allele"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="ref_allele_frequency"]<-"eaf.outcome" #eaf in disease GWAS
			gwasplusSNPs$strand.outcome <- "+" #confirmed through correspondence with Stavroula Kanoni 
			gwasplusSNPs$genome.build.outcome<-"hg18"  #build 36 us hg18
		}
		
		if(sum(grep("ilcco",database))==1){#coordinates in ILCCO are in build 37 / hg19 format
			gwasplusSNPs<-merge(gwasplusSNPs,unique(Bed),by.x="rs_number",by.y="bed")
			names(gwasplusSNPs)[names(gwasplusSNPs)=="snp"]<-"snp.outcome" #SNP rsid in disease GWAS	
			info.files<-c("EAF_all.csv","EAF_adeno.csv","EAF_squam.csv")
			database.files<-c("ilcco","ilcco_adeno","ilcco_squam")
			info.file<-info.files[database.files %in% database]
			system(paste("fgrep -hwf data/snplist.txt ",gwas.dir,"ILCCO/",info.file," > data/snp_info.ilcco.txt",sep="")) #quality metrics for ilcco SNPs
			system(paste("head -1 ",gwas.dir,"ILCCO/",info.file," > data/filehead.txt",sep=""))
			system("cat data/filehead.txt  data/snp_info.ilcco.txt  > data/temp ; mv data/temp  data/snp_info.ilcco.txt")
			info.ilcco<-read.table("data/snp_info.ilcco.txt",sep=",",head=T,colClasses="character")		
			names(info.ilcco)[names(info.ilcco)=="rs_number"]<-"snp.outcome"
			gwasplusSNPs<-merge(gwasplusSNPs,info.ilcco,by="snp.outcome")
			names(gwasplusSNPs)[names(gwasplusSNPs)=="EAF_s1"]<-"eaf.outcome" #eaf in disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="reference_allele"]<-"effect_allele" #effect allele in disease GWAS
			gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$OR))
			names(gwasplusSNPs)[names(gwasplusSNPs)=="OR_se"]<-"se.outcome" #effect.outcome se called OR se in ilcco		
			Vars<-unlist(strsplit(gwasplusSNPs$rs_number,split=":"))
			Chr_id<-Vars[seq(1,length(Vars),by=2)]
			Chr_id<-as.numeric(gsub("chr","",Chr_id))
			Chr_pos<-as.numeric(Vars[seq(2,length(Vars),by=2)])
			gwasplusSNPs$Chr_id.outcome<-Chr_id
			gwasplusSNPs$Chr_pos.outcome<-Chr_pos	
			gwasplusSNPs$genome.build.outcome<-"hg19"
			names(gwasplusSNPs)[names(gwasplusSNPs)=="N_control"]<-"N_control.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="N_case"]<-"N_case.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="r2.info.imp.outcome"]<-"r2.info.imp.outcome"
			#columns in info.ilcco: name  rs_number chr position info_s1  RSQ_s2 RSQ_s3 info_s4  EAF_s1 EAF_s2 EAF_s3 EAF_s4 N_case N_control
		}
		
		if(sum(grep("smoking",database))==1){#coordinates in hg18 format. INFO=1 for all (unavailable); OR is not an OR !!! Is the linear regression beta for the continuous variables and logistic regression beta for the discrete variables. 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNP"]<-"snp.outcome" #SNP rsid in diesase GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="OR"]<-"effect.outcome" #OR is not an OR !!! Is the linear regression beta for the continuous variables and logistic regression beta for the discrete variables. 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SE"]<-"se.outcome"
			names(gwasplusSNPs)[names(gwasplusSNPs)=="A1"]<-"effect_allele" #effect allele in disease GWAS, assuming A1 is effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="A2"]<-"other_allele" #assuming A2 is the other allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="FRQ_U"]<-"eaf.outcome" #eaf in disease GWAS
			
		}
			
		#GIANT		
		if(sum(grep("giant",database))==1){
			ea.ids<-c("A1","Allele1") #effect allele confirmed on GIANT website
			oa.ids<-c("A2","Allele2")
			snp.ids<-c("SNP","MarkerName") 
			eaf.ids<-c("Freq1.Hapmap","Freq.Allele1.HapMapCEU","FreqAllele1HapMapCEU")
			se.ids<-c("se","SE")
			beta.ids<-c("b") 
			head(gwasplusSNPs)
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% snp.ids]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% beta.ids]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% se.ids]<-"se.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% ea.ids]<-"effect_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% oa.ids]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% eaf.ids]<-"eaf.outcome" #SNP rsid in trait/disease GWAS
			gwasplusSNPs$genome.build.outcome<-"hg19"
			gwasplusSNPs$strand.outcome<-"+"
		}
		
		if(sum(grep("AMD",database))==1){
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Marker"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Overall"]<-"direction.outcome"
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Allele1"]<-"effect_allele" #assuming allele1 is the effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Allele2"]<-"other_allele" 
			gwasplusSNPs$effect.outcome<-NA
			gwasplusSNPs$se.outcome<-NA
			gwasplusSNPs$eaf.outcome<-NA
		}
		
		#ALS
		if(sum(grep("ALS_disease",database))==1){
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNP"]<-"snp.outcome" #SNP rsid in trait/disease GWAS		
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Freq_B"]<-"eaf.outcome" #assuming  Freq_B is the effect allele			
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Allele_B"]<-"effect_allele" #assuming Allele_B is the effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Allele_A"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Effect"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="StdErr"]<-"se.outcome" #se		
		}
		
		if(sum(grep("magic",database))==1){
			#minor allele is not the effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="snp"]<-"snp.outcome" #SNP rsid in trait/disease GWAS			
			names(gwasplusSNPs)[names(gwasplusSNPs)=="effect_allele"]<-"effect_allele" #effect allele in disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="other_allele"]<-"other_allele" #effect allele in disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="effect"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="stderr"]<-"se.outcome" #se
				
			#identify effect allele frequency using 1000 genomes european data
			freq.tab<-read.table(paste("data/allele_freq_1000genomes_european_",trait,".txt",sep=""),sep="\t",head=T,colClasses="character")
			gwasplusSNPs<-merge(gwasplusSNPs,freq.tab,by.x="snp.outcome",by.y="snp")
			gwasplusSNPs$eaf.outcome<-gwasplusSNPs$maf
			gwasplusSNPs$effect_allele<-tolower(gwasplusSNPs$effect_allele)
			gwasplusSNPs$other_allele<-tolower(gwasplusSNPs$other_allele)
			gwasplusSNPs$minor_allele<-tolower(gwasplusSNPs$minor_allele)
			gwasplusSNPs$major_allele<-tolower(gwasplusSNPs$major_allele)
			gwasplusSNPs$eaf.outcome[gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele]<-1-as.numeric(gwasplusSNPs$eaf.outcome[gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele])
						
			#find effect allele frequency for effect alleles coded using different strand from 1000 genomes
			pos.diff<-gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele & gwasplusSNPs$effect_allele!=gwasplusSNPs$major_allele
			gwasplusSNPs$strand.diff<-paste(pos.diff)
			strand1<-c("g","c","t","a")
			strand2<-c("c","g","a","t")
			minor_allele<-gwasplusSNPs$minor_allele[pos.diff]
			major_allele<-gwasplusSNPs$major_allele[pos.diff]
			eaf.outcome<-gwasplusSNPs$maf_EUR1k[pos.diff]
			effect_allele<-gwasplusSNPs$effect_allele[pos.diff]
			other_allele<-gwasplusSNPs$other_allele[pos.diff]
			minor_allele_other_strand<-unlist(lapply(1:length(minor_allele),FUN=function(i) strand2[strand1 %in% minor_allele[i]]))
			major_allele_other_strand<-unlist(lapply(1:length(major_allele),FUN=function(i) strand2[strand1 %in% major_allele[i]]))
			eaf.outcome[effect_allele != minor_allele_other_strand]<-1-as.numeric(eaf.outcome[effect_allele != minor_allele_other_strand])
			gwasplusSNPs$minor_allele[pos.diff]<-minor_allele_other_strand
			gwasplusSNPs$major_allele[pos.diff]<-major_allele_other_strand
			gwasplusSNPs$eaf.outcome[pos.diff]<-eaf.outcome
	#			gwasplusSNPs[,c("effect_allele","other_allele","minor_allele","major_allele","eaf.outcome","maf_EUR1k","strand.diff")]
			if(any(gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele & gwasplusSNPs$effect_allele!=gwasplusSNPs$major_allele)) stop(paste("effect alleles are on different strands in ",database," database and 1000genomes",sep=""))
			
			#set effect allele frequency to NA for palindromic SNPs
			pos.pal<-unique(unlist(lapply(1:length(strand1),FUN=function(i) which(gwasplusSNPs$effect_allele == strand1[i] & gwasplusSNPs$other_allele == strand2[i]))))
			gwasplusSNPs$eaf.outcome[pos.pal]<-NA
		}
		
		if(sum(grep("icbp",database))==1){
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNP.ID"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			#gwasplusSNPs$eaf.outcome<-0.1 #artificaially set to 0.1. Use maf of dbSNP to decide whether to drop palindromic SNPs
			#gwasplusSNPs$eaf.outcome[gwasplusSNPs$Coded.Allele != gwasplusSNPs$Minor.allele]<-0.9
			#gwasplusSNPs$eaf.outcome<-as.numeric(gwasplusSNPs$eaf.outcome)
			gwasplusSNPs$other_allele[gwasplusSNPs$Coded.Allele!=gwasplusSNPs$Allele1] <- gwasplusSNPs$Allele1[gwasplusSNPs$Coded.Allele != gwasplusSNPs$Allele1]
			gwasplusSNPs$other_allele[gwasplusSNPs$Coded.Allele!=gwasplusSNPs$Allele2] <- gwasplusSNPs$Allele2[gwasplusSNPs$Coded.Allele != gwasplusSNPs$Allele2]
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Coded.Allele"]<-"effect_allele" 			
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Effect"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SE"]<-"se.outcome" #se
			gwasplusSNPs$eaf.outcome<-NA
			databases<-c("icbp003588","icbp003589" ,"icbp003590","icbp003591")
			outcomes<-c("systolic_blood_pressure","diastolic_blood_pressure","pulse_pressure","mean_arterial_pressure")
			databases[databases==database]
	#		SDs<-c(19,11,14,13) #SDs for "systolic_blood_pressure","diastolic_blood_pressure","pulse_pressure","mean_arterial_pressure", respectively
	#		gwasplusSNPs$effect.outcome<-as.numeric(gwasplusSNPs$effect.outcome)/SDs[databases==database]
	#		gwasplusSNPs$se.outcome<-as.numeric(gwasplusSNPs$se.outcome)/SDs[databases==database]
			database<-outcomes[databases==database]
			#"icbp003588"/* Name: Genome-wide association analysis on systolic blood pressure Accession: pha003588.1*/
			#"icbp003589" /*Name: Genome-wide association analysis on diastolic blood pressure Accession: pha003589.1*/
			#"icbp003590" /*Name: Genome-wide association analysis on pulse pressure Accession: pha003590.1 */
			#"icbp003591" /*Name: Genome-wide association analysis on mean arterial pressure:pha003591.1 */
		}		
	
		
		
		if(sum(grep("glgc",database))==1){
			names(gwasplusSNPs)[names(gwasplusSNPs)=="rsid"]<-"snp.outcome" #SNP rsid in trait/disease GWAS	
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Freq.A1.1000G.EUR"]<-"eaf.outcome" #assuming A1 is the effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="A1"]<-"effect_allele" #confirmed on GLGC website
			names(gwasplusSNPs)[names(gwasplusSNPs)=="A2"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="beta"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="se"]<-"se.outcome" #se
			Vars<-unlist(strsplit(gwasplusSNPs$SNP_hg18,split=":"))
			Chr_id<-Vars[seq(1,length(Vars),by=2)]
			Chr_id<-as.numeric(gsub("chr","",Chr_id))
			Chr_pos<-as.numeric(Vars[seq(2,length(Vars),by=2)])
			gwasplusSNPs$Chr_id.disease<-Chr_id
			gwasplusSNPs$Chr_pos.disease<-Chr_pos
		}
		
		#ADHD
		if("pgc_adhd"==database){
			names(gwasplusSNPs)[names(gwasplusSNPs)=="snpid"]<-"snp.outcome" #SNP rsid in trait/disease GWAS		
			names(gwasplusSNPs)[names(gwasplusSNPs)=="CEUmaf"]<-"eaf.outcome" #assuming effect allele is the minor allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="a1"]<-"effect_allele" #assuming a1 is the effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="a2"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="zscore"]<-"zscore" 
			gwasplusSNPs$eaf.outcome[gwasplusSNPs$eaf.outcome=="."]<-NA
			gwasplusSNPs$effect.outcome<-NA
			gwasplusSNPs$se.outcome<-NA		
		}
		
		#schizophrenia
		if(sum(grep("pgc_scz",database))==1){
			names(gwasplusSNPs)[names(gwasplusSNPs)=="snpid"]<-"snp.outcome" #SNP rsid in trait/disease GWAS		
			names(gwasplusSNPs)[names(gwasplusSNPs)=="a1"]<-"effect_allele" #assuming a1 is the effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="a2"]<-"other_allele" 
			gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$or))
			names(gwasplusSNPs)[names(gwasplusSNPs)=="se"]<-"se.outcome" #se
			gwasplusSNPs$eaf.outcome<-NA
		}
		
		#Major depression disorder, dipolar disorder, autism
		if(any(c("pgc_bip","pgc_mdd","pgc_aut") %in% database)){ 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="snpid"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="CEUaf"]<-"eaf.outcome" #assuming CEUaf is effect allele frequency			
			names(gwasplusSNPs)[names(gwasplusSNPs)=="a1"]<-"effect_allele" #assuming a1 is the effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="a2"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="se"]<-"se.outcome" #se
	  gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$or))	
			gwasplusSNPs$eaf.outcome[gwasplusSNPs$eaf.outcome=="."]<-NA
		}
		
		#International Genomics of Alzheimer’s Project
		if(any(c("igap_stage1") %in% database)){ 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="MarkerName"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			gwasplusSNPs$eaf.outcome<-NA			
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Effect_allele"]<-"effect_allele"  
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Non_Effect_allele"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Beta"]<-"effect.outcome"
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SE"]<-"se.outcome" #se
		}
		
		#GCAN anorexia nervosa
		if(any(c("gcan_anorexia") %in% database)){ 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNP"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			gwasplusSNPs$eaf.outcome<-NA	
			names(gwasplusSNPs)[names(gwasplusSNPs)=="reference_allele"]<-"effect_allele"  
			names(gwasplusSNPs)[names(gwasplusSNPs)=="other_allele"]<-"other_allele" 
			gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$OR))
			names(gwasplusSNPs)[names(gwasplusSNPs)=="OR_se"]<-"se.outcome" #se
		}
		
		#International Inflammatory Bowel Disease Genetics Consortium crohns
		if("ibd_crohns" == database){ 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNP"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="ALLELE"] <- "effect_allele"
			gwasplusSNPs$eaf.outcome<-NA	
			gwasplusSNPs$other_allele<-NA 
	#			names(gwasplusSNPs)[names(gwasplusSNPs)=="SCAN"]<-"SCAN_Z" #SCAN is probably the pooled Z but not sure		
			gwasplusSNPs$effect.outcome<-NA
			gwasplusSNPs$se.outcome<-NA
		}
		
		if("ibd_ichip" %in% database){ 
			#there are duplicate SNPs in this database, and the effects aren't always consistent; something to do with the study design and use of proxy SNPs
			gwasplusSNPs$ICHIP_SNP[gwasplusSNPs$ICHIP_SNP=="same"]<-gwasplusSNPs$GWAS_SNP[gwasplusSNPs$ICHIP_SNP=="same"]
			gwasplusSNPs$snp.outcome<-gwasplusSNPs$GWAS_SNP
			gwasplusSNPs$snp.outcome[!gwasplusSNPs$GWAS_SNP %in% unique(trait.assoc$snp.exposure)]<-gwasplusSNPs$ICHIP_SNP[!gwasplusSNPs$GWAS_SNP %in% unique(trait.assoc$snp.exposure) ]
			if(any(!gwasplusSNPs$snp.outcome %in% unique(trait.assoc$snp.exposure))) stop("1 or more SNPs in lookup results for ibd_ichip not present in instruments file")
	
			allele.snps<-gwasplusSNPs$GWAS_A1A2
			allele.snps[!gwasplusSNPs$GWAS_SNP %in% unique(trait.assoc$snp.exposure)]<-gwasplusSNPs$ICHIP_A1A2[!gwasplusSNPs$GWAS_SNP %in% unique(trait.assoc$snp.exposure)]
			gwasplusSNPs$allele.snp<-allele.snps
			gwasplusSNPs$effect_allele<-substr(allele.snps,1,1)
			gwasplusSNPs$other_allele<-substr(allele.snps,2,2)
			
			#find effect allele frequency
			freq.tab<-read.table(paste("data/allele_freq_1000genomes_european_",trait,".txt",sep=""),sep="\t",head=T,colClasses="character")
			gwasplusSNPs<-merge(gwasplusSNPs,freq.tab,by.x="snp.outcome",by.y="snp")
			gwasplusSNPs$eaf.outcome<-gwasplusSNPs$maf_EUR1k
			gwasplusSNPs$effect_allele<-tolower(gwasplusSNPs$effect_allele)
			gwasplusSNPs$other_allele<-tolower(gwasplusSNPs$other_allele)
			gwasplusSNPs$minor_allele<-tolower(gwasplusSNPs$minor_allele)
			gwasplusSNPs$major_allele<-tolower(gwasplusSNPs$major_allele)
			gwasplusSNPs$eaf.outcome[gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele]<-1-as.numeric(gwasplusSNPs$eaf.outcome[gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele])
		
			#find effect allele frequency for effect alleles coded using different strand from 1000 genomes
			pos.diff<-gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele & gwasplusSNPs$effect_allele!=gwasplusSNPs$major_allele
			gwasplusSNPs$strand.diff<-paste(pos.diff)
			strand1<-c("g","c","t","a")
			strand2<-c("c","g","a","t")
			minor_allele<-gwasplusSNPs$minor_allele[pos.diff]
			major_allele<-gwasplusSNPs$major_allele[pos.diff]
			eaf.outcome<-gwasplusSNPs$maf_EUR1k[pos.diff]
			effect_allele<-gwasplusSNPs$effect_allele[pos.diff]
			other_allele<-gwasplusSNPs$other_allele[pos.diff]
			minor_allele_other_strand<-unlist(lapply(1:length(minor_allele),FUN=function(i) strand2[strand1 %in% minor_allele[i]]))
			major_allele_other_strand<-unlist(lapply(1:length(major_allele),FUN=function(i) strand2[strand1 %in% major_allele[i]]))
			eaf.outcome[effect_allele != minor_allele_other_strand]<-1-as.numeric(eaf.outcome[effect_allele != minor_allele_other_strand])
			gwasplusSNPs$minor_allele[pos.diff]<-minor_allele_other_strand
			gwasplusSNPs$major_allele[pos.diff]<-major_allele_other_strand
			gwasplusSNPs$eaf.outcome[pos.diff]<-eaf.outcome
	#			gwasplusSNPs[,c("effect_allele","other_allele","minor_allele","major_allele","eaf.outcome","maf_EUR1k","strand.diff")]
			if(any(gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele & gwasplusSNPs$effect_allele!=gwasplusSNPs$major_allele)) stop(paste("effect alleles are on different strands in ",database," database and 1000genomes",sep=""))
			
			#set effect allele frequency to NA for palindromic SNPs
			pos.pal<-unique(unlist(lapply(1:length(strand1),FUN=function(i) which(gwasplusSNPs$effect_allele == strand1[i] & gwasplusSNPs$other_allele == strand2[i]))))
			gwasplusSNPs$eaf.outcome[pos.pal]<-NA
			
			or.ids<-c("CD_OR","UC_OR","IBD_OR") 
			ci.ids<-c("CD_OR_CI","UC_OR_CI","IBD_OR_CI") 
			p.ids<-c("CD_META_P","UC_META_P","IBD_META_P")
	
			gwas.diseases<-NULL
			for(i in 1:3) { 
				OR<-or.ids[i]
				CI<-ci.ids[i]
				P<-p.ids[i]
				names.gwas<-names(gwasplusSNPs)[!names(gwasplusSNPs) %in% c(or.ids,ci.ids,p.ids)]
				names.gwas<-c(names.gwas,OR,CI,P)
				gwas.disease<-gwasplusSNPs[,names.gwas]
				pos.end<-regexpr("_OR",OR)-1
				gwas.disease$outcome<-substr(OR,1,pos.end)
				names(gwas.disease)[names(gwas.disease)==OR]<-"OR"
				names(gwas.disease)[names(gwas.disease)==CI]<-"CI"
				names(gwas.disease)[names(gwas.disease)==P]<-"p.value.outcome"
				gwas.diseases[[i]]<-gwas.disease
			}
			
			gwasplusSNPs<-do.call("rbind",gwas.diseases)
			gwasplusSNPs$OR[gwasplusSNPs$OR=="-"]<-NA
			gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$OR))
			gwasplusSNPs$CI[gwasplusSNPs$CI=="-"]<-NA
			gwasplusSNPs<-gwasplusSNPs[!is.na(gwasplusSNPs$CI),]
			CIs<-unlist(strsplit(gwasplusSNPs$CI,split="-"))
			lci<-CIs[seq(1,length(CIs),2)]
			uci<-CIs[seq(2,length(CIs),2)]
			gwasplusSNPs$se.outcome<-(log(as.numeric(uci))-log(as.numeric(lci)))/(1.96*2)    
		}
			
		#International Inflammatory Bowel Disease Genetics Consortium colitis
		if(any(c("ibd_ucolitis") %in% database)){ 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNP"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs)=="RISK"]<-"effect_allele" 
			gwasplusSNPs$eaf.outcome<-NA #COF could be EAF
			gwasplusSNPs$other_allele<-NA 
			Vars<-unlist(strsplit(gwasplusSNPs$OR.95..,split="\\("))
			OR<-Vars[seq(1,length(Vars),by=2)]
			gwasplusSNPs$effect.outcome<-log(as.numeric(OR))
			CI<-Vars[seq(2,length(Vars),by=2)]
			CIs<-unlist(strsplit(CI,split="-"))
			LCI<-as.numeric(CIs[seq(1,length(CIs),by=2)])
			UCI<-CIs[seq(2,length(CIs),by=2)]
			UCI<-as.numeric(gsub(")","",UCI))
			gwasplusSNPs$se.outcome<-(log(UCI)-log(LCI))/(1.96*2)
		}
		
		#rheumatoid_arthritis
		if(any(c("rheumatoid_arthritis") %in% database)){ 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SNP"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			gwasplusSNPs$eaf.outcome<-NA			
			names(gwasplusSNPs)[names(gwasplusSNPs)=="minor_al"]<-"effect_allele"  
			names(gwasplusSNPs)[names(gwasplusSNPs)=="major_al"]<-"other_allele" 
			gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$meta_OR ))
			gwasplusSNPs$se.outcome<-(log(as.numeric(gwasplusSNPs$OR_95.CI_up ))-log(as.numeric(gwasplusSNPs$OR_95.CI_lo )))/(1.96*2)
			gwasplusSNPs$N_control<-as.numeric(gwasplusSNPs$controls_MM)+as.numeric(gwasplusSNPs$controls_Mm) +as.numeric(gwasplusSNPs$controls_mm)
			gwasplusSNPs$N_case<-as.numeric(gwasplusSNPs$cases_MM)+as.numeric(gwasplusSNPs$cases_Mm) +as.numeric(gwasplusSNPs$cases_mm)
			gwasplusSNPs$N_total<-gwasplusSNPs$N_control+gwasplusSNPs$N_case
			gwasplusSNPs<-gwasplusSNPs[gwasplusSNPs$se.outcome!=0,]
		}
		
		#GEnetic Factors for OSteoporosis Consortium
		if(any(c("gefos_LSbmdpooled","gefos_FNbmdpooled") %in% database)){ 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="MarkerName"]<-"snp.outcome" #SNP rsid in trait/disease GWAS			
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Freq.Allele1.HapMapCEU"]<-"eaf.outcome" #assuming Allele1 is the effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Allele1"]<-"effect_allele" #assuming allele1 is the effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Allele2"]<-"other_allele" 
			gwasplusSNPs$effect.outcome<-NA
			gwasplusSNPs$se.outcome<-NA
		}
		
		#SSGAC MA_EA_1st_stage, Social Science Genetic Association Consortium
		if(any(c("ssgac_educ_attain") %in% database)){ 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="MarkerName"]<-"snp.outcome" #SNP rsid in trait/disease GWAS		
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Allele1"]<-"effect_allele" #assuming allele1 is the effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Allele2"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="BetaZScale"]<-"zscore" 
			gwasplusSNPs$eaf.outcome<-NA
			gwasplusSNPs$effect.outcome<-NA
			gwasplusSNPs$se.outcome<-NA					
		}
		
		#SSGAC MA_EA_1st_stage, Social Science Genetic Association Consortium
		if(any(c("ssgac_eduyrs","ssgac_college") %in% database)){ 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="MarkerName"]<-"snp.outcome" #SNP rsid in trait/disease GWAS			
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Effect_Allele"]<-"effect_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Other_Allele"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="EAF"]<-"eaf.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="Beta"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="OR"]<-"OR" 
			names(gwasplusSNPs)[names(gwasplusSNPs)=="SE"]<-"se.outcome" 
			if(all(names(gwasplusSNPs)!="effect.outcome")) {
				gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs[,names(gwasplusSNPs) =="OR"]))
			}
		}
		
		#EGG (Early Growth Genetics) Consortium
		if(any(c("egg_bw2","egg_bl","egg_hc","egg_obesity") %in% database)){ 
			snp.ids<-c("RSID","SNP")
			eaf.ids<-c("EAF","Freq","EFFECT_ALLELE_FREQ")		
			ea.ids<-c("EA","EFFECT_ALLELE")
			oa.ids<-c("NEA","OTHER_ALLELE","NON_EFFECT_ALLELE")
			beta.ids<-c("BETA","EFFECT")
			se.ids<-c("SE","se")		
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% snp.ids]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% eaf.ids]<-"eaf.outcome" #SNP rsid in trait/disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% ea.ids]<-"effect_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% oa.ids]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% beta.ids]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) %in% se.ids]<-"se.outcome" 					
			if(all(names(gwasplusSNPs)!="eaf.outcome")) gwasplusSNPs$eaf.outcome<-NA #no eaf for "EGG/EGG_BW2_DISCOVERY.txt.gz", "EGG/EGG-GWAS-BL.txt.gz", "EGG/EGG_Obesity_Meta_Analysis_1.txt.gz"
			gwasplusSNPs$eaf.outcome<-gsub(",",".",gwasplusSNPs$eaf.outcome) #EGG_HC_DISCOVERY.txt.gz has commas instead of full stops 
			gwasplusSNPs$effect.outcome<- gsub(",",".",gwasplusSNPs$effect.outcome) 		#EGG_HC_DISCOVERY.txt.gz has commas instead of full stops 
			gwasplusSNPs$se.outcome<- gsub(",",".",gwasplusSNPs$se.outcome) #EGG_HC_DISCOVERY.txt.gz has commas instead of full stops 
			gwasplusSNPs$p.value.outcome<- gsub(",",".",gwasplusSNPs$p.value.outcome) #EGG_HC_DISCOVERY.txt.gz has commas instead of full stops 
			if(database=="egg_hc"){
				gwasplusSNPs$N_total.outcome<-gsub(",",".",gwasplusSNPs$N_total.outcome) #EGG_HC_DISCOVERY.txt.gz has commas instead of full stops 
			}
		}
		     
		#GWAScatalogplusMetab assume all outcomes are diseases
		if("gwasCatplusMetab_disease" %in% database){ 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "SNPs"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Risk.Allele.Frequency"]<-"eaf.outcome" #SNP rsid in trait/disease GWAS			
			names(gwasplusSNPs)[names(gwasplusSNPs) == "effect.allele.GwCat"]<-"effect_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "other.allele.ens.GwCat"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "lnbeta"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "lnse"]<-"se.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "PUBMEDID"]<-"outcome.pubmedid" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "ci95"]<-"outcome.ci95" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "First.Author"]<-"outcome.author" 	
			gwasplusSNPs$outcome<-gwasplusSNPs$Disease.Trait
	
		}     
		
		#GWAScatalogplusMetab assume all outcomes are diseases
		if(any(c("gwasCatplusMetab_notdisease") %in% database)){ 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "SNPs"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Risk.Allele.Frequency"]<-"eaf.outcome" #SNP rsid in trait/disease GWAS			
			names(gwasplusSNPs)[names(gwasplusSNPs) == "effect.allele.GwCat"]<-"effect_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "other.allele.ens.GwCat"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "beta"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "se"]<-"se.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "PUBMEDID"]<-"outcome.pubmedid" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "ci95"]<-"outcome.ci95" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "First.Author"]<-"outcome.author" 
			gwasplusSNPs$outcome<-gwasplusSNPs$Disease.Trait
		}     
		
		if(sum(grep("dbGAP",database))==1){ 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "SNP.ID"]<-"snp.outcome"
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Odds.ratio"]<-"OR"
			names(gwasplusSNPs)[names(gwasplusSNPs) == "CI.low"]<-"OR.lci"
			names(gwasplusSNPs)[names(gwasplusSNPs) == "CI.high"]<-"OR.uci"	
			gwasplusSNPs$effect_allele<-gwasplusSNPs$Allele1
			gwasplusSNPs$other_allele<-gwasplusSNPs$Allele2
			if(any(names(gwasplusSNPs)=="Code.Allele")){
				gwasplusSNPs$other_allele<-NA
				gwasplusSNPs$other_allele[gwasplusSNPs$Code.Allele!=gwasplusSNPs$Allele1] <- gwasplusSNPs$Allele1[gwasplusSNPs$Code.Allele != gwasplusSNPs$Allele1]
				gwasplusSNPs$other_allele[gwasplusSNPs$Code.Allele!=gwasplusSNPs$Allele2] <- gwasplusSNPs$Allele2[gwasplusSNPs$Code.Allele != gwasplusSNPs$Allele2]
				gwasplusSNPs$effect_allele<-gwasplusSNPs$Code.Allele
			}
			gwasplusSNPs$effect.outcome<-NA
			gwasplusSNPs$se.outcome<-NA
			if("OR" %in% names(gwasplusSNPs)){
				gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$OR))
				if("OR.lci" %in% names(gwasplusSNPs)){
					gwasplusSNPs$se.outcome<-(log(as.numeric(gwasplusSNPs$OR.uci)) - log(as.numeric(gwasplusSNPs$OR.lci)))/(1.96*2)
					
				}
				if(!"OR.lci" %in% names(gwasplusSNPs) | all(gwasplusSNPs$OR.uci=="")){
					z<-qnorm(as.numeric(gwasplusSNPs$p.value.outcome)/2,lower.tail=F)
					gwasplusSNPs$se.outcome<-abs(gwasplusSNPs$effect.outcome/z)
				}
				
			}
	
			if("effect_allele" %in% names(gwasplusSNPs)){ 
				freq.tab<-read.table(paste("data/allele_freq_1000genomes_european_",trait,".txt",sep=""),sep="\t",head=T,colClasses="character")
				gwasplusSNPs<-merge(gwasplusSNPs,freq.tab,by.x="snp.outcome",by.y="snp")
				gwasplusSNPs$eaf.outcome<-gwasplusSNPs$maf
				gwasplusSNPs$effect_allele<-tolower(gwasplusSNPs$effect_allele)
				gwasplusSNPs$other_allele<-tolower(gwasplusSNPs$other_allele)
				gwasplusSNPs$minor_allele<-tolower(gwasplusSNPs$minor_allele)
				gwasplusSNPs$major_allele<-tolower(gwasplusSNPs$major_allele)
				gwasplusSNPs$eaf.outcome[gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele]<-1-as.numeric(gwasplusSNPs$eaf.outcome[gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele])
				
				#find allele frequency for effect alleles coded using different strand from 1000 genomes
				pos.diff<-gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele & gwasplusSNPs$effect_allele!=gwasplusSNPs$major_allele
				gwasplusSNPs$strand.diff<-paste(pos.diff)
				strand1<-c("g","c","t","a")
				strand2<-c("c","g","a","t")
				minor_allele<-gwasplusSNPs$minor_allele[pos.diff]
				major_allele<-gwasplusSNPs$major_allele[pos.diff]
				eaf.outcome<-gwasplusSNPs$maf_EUR1k[pos.diff]
				effect_allele<-gwasplusSNPs$effect_allele[pos.diff]
				other_allele<-gwasplusSNPs$other_allele[pos.diff]
				minor_allele_other_strand<-unlist(lapply(1:length(minor_allele),FUN=function(i) strand2[strand1 %in% minor_allele[i]]))
				major_allele_other_strand<-unlist(lapply(1:length(major_allele),FUN=function(i) strand2[strand1 %in% major_allele[i]]))
				eaf.outcome[effect_allele != minor_allele_other_strand]<-1-as.numeric(eaf.outcome[effect_allele != minor_allele_other_strand])
				gwasplusSNPs$minor_allele[pos.diff]<-minor_allele_other_strand
				gwasplusSNPs$major_allele[pos.diff]<-major_allele_other_strand
				gwasplusSNPs$eaf.outcome[pos.diff]<-eaf.outcome
	#			gwasplusSNPs[,c("effect_allele","other_allele","minor_allele","major_allele","eaf.outcome","maf_EUR1k","strand.diff")]
				if(any(gwasplusSNPs$effect_allele!=gwasplusSNPs$minor_allele & gwasplusSNPs$effect_allele!=gwasplusSNPs$major_allele)) stop(paste("effect alleles are on different strands in ",database," database and 1000genomes",sep=""))
				
				#set effect allele frequency to NA for palindromic SNPs
				pos.pal<-unique(unlist(lapply(1:length(strand1),FUN=function(i) which(gwasplusSNPs$effect_allele == strand1[i] & gwasplusSNPs$other_allele == strand2[i]))))
				gwasplusSNPs$eaf.outcome[pos.pal]<-NA
			}
		
			
			if(!"effect_allele" %in% names(gwasplusSNPs) & !"other_allele" %in%  names(gwasplusSNPs)){ # if effect allele and other allele are missing assume effect allele was the minor allele and calculate missing allele info using 1000 genomes european data
				freq.tab<-read.table(paste("data/allele_freq_1000genomes_european_",trait,".txt",sep=""),sep="\t",head=T,colClasses="character")
				gwasplusSNPs<-merge(gwasplusSNPs,freq.tab,by.x="snp.outcome",by.y="snp")
				names(gwasplusSNPs)[names(gwasplusSNPs)=="minor_allele"]<-"effect_allele"
				names(gwasplusSNPs)[names(gwasplusSNPs)=="major_allele"]<-"other_allele"
				names(gwasplusSNPs)[names(gwasplusSNPs)=="maf_EUR1k"]<-"eaf.outcome"
			}
			
		}
		
		if(sum(grep("Haem",database))==1){ 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "SNP"]<-"snp.outcome" #SNP rsid in trait/disease GWAS
			gwasplusSNPs$eaf.outcome<-NA
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Effect_Allele"]<-"effect_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Non_Effect_Allele"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Beta"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) ==  "SE"]<-"se.outcome" 	
		}
		
		if(sum(unlist(lapply(c("immunobase","t1dbase"),FUN=function(x) grep(x,database))))==1){ 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Marker"]<-"snp.outcome" #SNP rsid in trait/disease GWAS		
			#Immunobase has no effect allele information. It seems that the minor allele was probably the effect allele http://hmg.oxfordjournals.org/content/21/23/5202.long
			#infer the effect allele information using 1000geomes european data
			freq.tab<-read.table(paste("data/allele_freq_1000genomes_european_",trait,".txt",sep=""),sep="\t",head=T,colClasses="character")
			gwasplusSNPs<-merge(gwasplusSNPs,freq.tab,by.x="snp.outcome",by.y="snp")
			names(gwasplusSNPs)[names(gwasplusSNPs)=="minor_allele"]<-"effect_allele"
			names(gwasplusSNPs)[names(gwasplusSNPs)=="major_allele"]<-"other_allele"
			names(gwasplusSNPs)[names(gwasplusSNPs)=="maf_EUR1k"]<-"eaf.outcome"
			
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Odds Ratio"]<-"OR"
			if(sum(log(as.numeric(gwasplusSNPs$OR)))!=0){
				gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$OR))
				z<-qnorm(as.numeric(gwasplusSNPs$PValue)/2,lower.tail=F)
				gwasplusSNPs$se.outcome<-gwasplusSNPs$effect.outcome/z
			}
			if(!"effect.outcome" %in% names(gwasplusSNPs)){
				gwasplusSNPs$effect.outcome<-NA
				gwasplusSNPs$se.outcome<-NA
			}
		}       
	
		if(sum(grep("melanoma",database))==1){ 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "SNP"]<-"snp.outcome" #SNP rsid in trait/disease GWAS		
			names(gwasplusSNPs)[names(gwasplusSNPs) == "A1"]<-"effect_allele"  #assuming A1 effect allele
			gwasplusSNPs$other_allele<-NA 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "OR"]<-"OR"
			gwasplusSNPs$effect.outcome<-log(as.numeric(gwasplusSNPs$OR))
			names(gwasplusSNPs)[names(gwasplusSNPs) == "SE"]<-"se.outcome"	
			nums<-gregexpr("[:0-9:]",gwasplusSNPs$P.par)
			End<-lapply(1:length(nums),FUN=function(x) unlist(nums[x])[length(unlist(nums[x]))])
			gwasplusSNPs$p.value.outcome<-substr(gwasplusSNPs$P.par,1,End)
			gwasplusSNPs$eaf.outcome<-NA
		}
	
		if(sum(grep("CKDGen",database))==1){ 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "rsID"]<-"snp.outcome" #SNP rsid in trait/disease GWAS		
			names(gwasplusSNPs)[names(gwasplusSNPs) == "allele1"]<-"effect_allele"  #assuming allele1 effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs) == "allele2"]<-"other_allele"  #assuming allele1 effect allele
			names(gwasplusSNPs)[names(gwasplusSNPs) == "freqA1"]<-"eaf.outcome"  #assuming allele1 effect allele
			gwasplusSNPs$effect.outcome<-NA
			gwasplusSNPs$se.outcome<-NA
		}
		
		if(sum(grep("metabolome",database))==1){ 			
			gwasplusSNPs$genome.build.outcome<-"hg18" # hg and b36 equivalent 
			gwasplusSNPs$strand.outcome<-"forward"
			tempvar<-unlist(strsplit(gwasplusSNPs$MarkerName,split=":"))
			gwasplusSNPs$metabolonID<-tempvar[seq(1,length(tempvar),2)]
			gwasplusSNPs$snp.outcome<-tempvar[seq(2,length(tempvar),2)]
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Allele1"]<-"effect_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Allele2"]<-"other_allele" 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "Effect"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) ==  "StdErr"]<-"se.outcome" 	
			names(gwasplusSNPs)[names(gwasplusSNPs) ==  "Freq1"]<-"eaf.outcome" 
	#			map.tab<-read.table(paste(gwas.dir,"metabolome/metaboliteMap.txt",sep=""),sep="\t",head=T,colClasses="character",quote = "")
			metab.info<-read.table(paste(gwas.dir,"metabolome/metabolic_traits_metabolome.txt",sep=""),sep="\t",head=T,stringsAsFactors=F,quote = "")
	#			head(metab.info)
			sd1<-metab.info$SD
			sd2<-metab.info$SD.1
			n1<-gsub("\"","",metab.info$N)
			n1<-gsub(",","",n1)
			n1<-gsub("-","",n1)
			n1<-as.numeric(n1)
			n2<-gsub("\"","",metab.info$N.1)
			n2<-gsub(",","",n2)
			n2<-gsub("-","",n2)
			n2<-as.numeric(n2)
			metab.info$SDp<-sqrt(((sd1^2*(n1-1))+(sd2^2*(n2-1)))/(n2+n1-2))
			pos.trait.start<-regexpr("metabolites_meta/",gwasplusSNPs$metabolonID)+17
			pos.trait.end<-regexpr(".metal",gwasplusSNPs$metabolonID)-1
			gwasplusSNPs$Metabolite.ID<-substr(gwasplusSNPs$metabolonID,pos.trait.start,pos.trait.end)
			gwasplusSNPs<-merge(gwasplusSNPs,metab.info,by="Metabolite.ID") #M01642 dropped because missing from metab.info (31 SNP associations); also missing from "metabolome/metaboliteMap.txt"
			gwasplusSNPs$effect.outcome<-as.numeric(gwasplusSNPs$effect.outcome)/gwasplusSNPs$SDp #standardized effect into SD units
			gwasplusSNPs$se.outcome<-as.numeric(gwasplusSNPs$se.outcome)/gwasplusSNPs$SDp #standardized effect into SD units
			gwasplusSNPs$outcome<-gwasplusSNPs$Metabolite.ID
	
	#			names(gwasplusSNPs)[names(gwasplusSNPs) ==  "metabolonDescription"]<-"outcome" 
	#			gwasplusSNPs$Metabolite.ID
	#			gwasplusSNPs$Metabolite
	#			gwasplusSNPs$Metabolite<-gsub("\"","",gwasplusSNPs$Metabolite)
	#			gwasplusSNPs$Pathway<-gsub("\"","",gwasplusSNPs$Pathway)
	#			gwasplusSNPs$Super.pathway
	#			gwasplusSNPs$outcome<-paste(gwasplusSNPs$Metabolite.ID,"_",gwasplusSNPs$Metabolite,"_",gwasplusSNPs$Pathway,"_",gwasplusSNPs$Super.pathway,sep="")
		}
		
		if(sum(grep("GTEx",database))==1){ 		
			tissue.snp<-unlist(strsplit(gwasplusSNPs$SNP,split="eqtl:"))
			tissue<-tissue.snp[seq(1,length(tissue.snp),2)]
			start.pos<-regexpr("alldata/",	tissue)+8
			stop.pos<-regexpr(".cis",	tissue)-1
			tissue<-substr(tissue,start.pos,stop.pos)
			gene.id<-gwasplusSNPs$gene
#			gene.name<-gwasplusSNPs$Gene_Name
			gwasplusSNPs$outcome<-paste(gene.id,tissue,sep="_")
			gwasplusSNPs$snp.outcome<-tissue.snp[seq(2,length(tissue.snp),2)]
			names(gwasplusSNPs)[names(gwasplusSNPs) == "beta"]<-"effect.outcome" 	
			gwasplusSNPs$se.outcome<-as.numeric(gwasplusSNPs$effect.outcome)/as.numeric(gwasplusSNPs$t.stat)
			system(paste("fgrep -hwf data/snplist.txt ",gwas.dir,"GTEx/GTEx_genot_imputed_variants_info4_maf05_CR95_CHR_POSb37_ID_REF_ALT.txt  > tempfile",sep=""))
			system(paste("head -1 ",gwas.dir,"GTEx/GTEx_genot_imputed_variants_info4_maf05_CR95_CHR_POSb37_ID_REF_ALT.txt > filehead.txt",sep=""))
			system("cat filehead.txt  tempfile > temp ; mv temp  tempfile")
			map.tab<-read.table("tempfile", comment.char="",sep=" ",head=T,colClasses="character")
			gwasplusSNPs<-merge(gwasplusSNPs,map.tab,by.x="snp.outcome",by.y="ID")
			names(gwasplusSNPs)[names(gwasplusSNPs) == "REF"]<-"other_allele" 	#confirmed by correspondence 
			names(gwasplusSNPs)[names(gwasplusSNPs) == "ALT"]<-"effect_allele"  #confirmed by correspondence
			gwasplusSNPs$eaf.outcome<-NA #effect allele is not the minor allele (confirmed by correspondence)
			gwasplusSNPs$genome.build.outcome<-"hg19" #2.  Our definitions of REF and ALT are always with respect to the hg19/GRCh37 genome reference.  The effect size we report is with respect to an increase or decrease in expression among individuals with ALT genotypes.
			gwasplusSNPs$strand.outcome <- "+" # The reference strand is always the forward or positive strand -- the same strand used by the hg19 genome reference.
		}
		
		if("ImmuneCellScience" %in% database){	
			trait.snp<-unlist(strsplit(gwasplusSNPs$name,split=".txt:"))
			trait.names<-trait.snp[seq(1,length(trait.snp),2)]
			pos.start<-regexpr("GWASResults/",gwasplusSNPs$name)+13
			pos.end<-nchar(trait.names)
			gwasplusSNPs$outcome<-substr(trait.names,pos.start,pos.end)
			gwasplusSNPs$snp.outcome<-trait.snp[seq(2,length(trait.snp),2)]
			#info.tab<-read.table(paste(gwas.dir,"ImmuneCellScience/4_Trait_Analysis/4_Trait_Analysis.txt",sep=""),sep="\t",colClasses="character",head=T)
			names(gwasplusSNPs)[names(gwasplusSNPs) == "allele1"]<-"other_allele" 	
			names(gwasplusSNPs)[names(gwasplusSNPs) == "effallele"]<-"effect_allele"  #effallele seems to be same as allele2
			names(gwasplusSNPs)[names(gwasplusSNPs) == "beta"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) ==  "sebeta"]<-"se.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) ==  "effallelefreq"]<-"eaf.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) ==  "build"]<-"genome.build.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) ==  "strand"]<-"strand.outcome" 
		}		
	
		if("Uric_acid" %in% database){		
			names(gwasplusSNPs)[names(gwasplusSNPs) == "variant"]<-"snp.outcome" #SNP rsid in trait/disease GWAS		
			names(gwasplusSNPs)[names(gwasplusSNPs) == "A2"]<-"other_allele" 	#confirmed in readme file, A1, allele for which effect (beta) is reported, A2, alternate allele (hg18, plus strand).
			names(gwasplusSNPs)[names(gwasplusSNPs) == "A1"]<-"effect_allele"  #confirmed in readme file, A1, allele for which effect (beta) is reported, A2, alternate allele (hg18, plus strand).
			names(gwasplusSNPs)[names(gwasplusSNPs) == "beta"]<-"effect.outcome" 
			names(gwasplusSNPs)[names(gwasplusSNPs) ==  "se"]<-"se.outcome" 
			gwasplusSNPs$eaf.outcome<-NA
			gwasplusSNPs$genome.build.outcome<-"hg18"
			gwasplusSNPs$strand.outcome<-"+"
			#Meta-analysis mean effect size (beta) is the inverse-variance weighted estimate derived from individual discovery study; se its standard error, meta-analysis P-value is P. 
			#Meta-analysis estimates are corrected for inflation of test statistics using genomic control at the individual study level. 
		}
	
		if(!database %in% c("ibd_ichip","telomere_outcomes","metabolome","gwasCatplusMetab_notdisease","gwasCatplusMetab_disease","GTEx","ImmuneCellScience")){ # these are databases that have multiple outcome traits
			gwasplusSNPs$outcome<-database
		}	     	
		
		gwasplusSNPs$database.outcome<-database
	
		if(all(!names(gwasplusSNPs) %in% "database.outcome")) stop("database.outcome column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "outcome")) stop("outcome column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "snp.outcome")) stop("snp.outcome column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "eaf.outcome")) stop("eaf.outcome column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "effect_allele")) stop("effect_allele column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "other_allele")) stop("other_allele column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "effect.outcome")) stop("effect.outcome column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "se.outcome")) stop("se.outcome column missing from gwasplusSNPs table")
		
		if(all(names(gwasplusSNPs)!="P.hwe.outcome")) gwasplusSNPs$P.hwe.outcome<-NA
		if(all(names(gwasplusSNPs)!="het.p.value.outcome")) gwasplusSNPs$het.p.value.outcome<-NA
		if(all(names(gwasplusSNPs)!="r2.info.imp.outcome")) gwasplusSNPs$r2.info.imp.outcome<-NA
		if(all(names(gwasplusSNPs)!="genotype.call.rate.outcome")) gwasplusSNPs$genotype.call.rate.outcome<-NA
		if(all(names(gwasplusSNPs)!="N_total.outcome")) gwasplusSNPs$N_total.outcome<-NA
		if(all(names(gwasplusSNPs)!="N_studies.outcome")) gwasplusSNPs$N_studies.outcome <-NA
		if(all(names(gwasplusSNPs)!="N_case.outcome")) gwasplusSNPs$N_case.outcome<-NA
		if(all(names(gwasplusSNPs)!="N_control.outcome")) gwasplusSNPs$N_control.outcome<-NA
		if(all(names(gwasplusSNPs)!="Chr_id.outcome")) gwasplusSNPs$Chr_id.outcome<-NA
		if(all(names(gwasplusSNPs)!="Chr_pos.outcome")) gwasplusSNPs$Chr_pos.outcome <-NA
		if(all(names(gwasplusSNPs)!="p.value.outcome")) gwasplusSNPs$p.value.outcome <-NA
		if(all(names(gwasplusSNPs)!="genome.build.outcome")) gwasplusSNPs$genome.build.outcome<-NA
		if(all(names(gwasplusSNPs)!="strand.outcome")) gwasplusSNPs$strand.outcome <-NA
		if(all(names(gwasplusSNPs)!="possible_strand.outcome")) gwasplusSNPs$possible_strand.outcome<-NA
		if(all(names(gwasplusSNPs)!="possible_build.outcome")) gwasplusSNPs$possible_build.outcome <-NA
		if(all(names(gwasplusSNPs)!="direction.outcome")) gwasplusSNPs$direction.outcome <-NA
		if(all(names(gwasplusSNPs)!="r2.proxy.outcome")) gwasplusSNPs$r2.proxy.outcome <-NA
		if(all(names(gwasplusSNPs)!="lambda.outcome")) gwasplusSNPs$lambda.outcome <-NA
		if(all(names(gwasplusSNPs)!="pgc.outcome")) gwasplusSNPs$pgc.outcome <-NA
		
		gwasplusSNPs<-gwasplusSNPs[,c("database.outcome","outcome","snp.outcome","effect_allele","other_allele","eaf.outcome","effect.outcome","se.outcome","p.value.outcome","Chr_id.outcome","Chr_pos.outcome","het.p.value.outcome","P.hwe.outcome","r2.info.imp.outcome","genotype.call.rate.outcome","N_total.outcome","N_studies.outcome","N_case.outcome","N_control.outcome" ,"genome.build.outcome","strand.outcome","possible_strand.outcome","possible_build.outcome","direction.outcome","r2.proxy.outcome","lambda.outcome","pgc.outcome")]
		
		eaf.outcome.notmiss.study<-T
		if(all(is.na(gwasplusSNPs$eaf.outcome))){
			eaf.outcome.notmiss.study<-F #don't have eaf for some studies, e.g these studies c("diagram","icbp","pgc_scz","gcan_anorexia","ibd_ucolitis","rheumatoid_arthritis") 
		}
		
		res.tab<-merge(trait.assoc,gwasplusSNPs,by.x="snp.exposure",by.y="snp.outcome")
		if(sum(grep("\\.y",names(res.tab)))) stop(paste("some column names are identical in the instruments file and the outcome database; ",names(res.tab)[grep(".y",names(res.tab))]," present in both files",sep=""))
	
		res.tab$effect.exposure<-as.numeric(res.tab$effect.exposure)
		res.tab$effect.outcome<-as.numeric(res.tab$effect.outcome)
		res.tab$eaf.trait<-as.numeric(res.tab$eaf.trait)
		res.tab$eaf.outcome[res.tab$eaf.outcome=="NR"]<-NA
		res.tab$eaf.outcome[res.tab$eaf.outcome=="NR "]<-NA
		res.tab$eaf.outcome<-as.numeric(res.tab$eaf.outcome)
		res.tab$eaf.outcome[res.tab$eaf.outcome==1]<-NA
		    
		#res.tab$p.value.outcome<-as.numeric(res.tab$p.value.outcome)
		#code mSNP effect so that effect allele is the allele that increases the trait
		#res.tab[,c("snp.exposure","effect.exposure","effect.outcome","allele1","allele2","effect_allele","other_allele","eaf.trait","eaf.outcome","info_s1","RSQ_s2","RSQ_s3","info_s4","q_p.value")]
		
		pos.change<-which(res.tab$effect.exposure<0)
		res.tab$eaf.trait[pos.change]<-1-res.tab$eaf.trait[pos.change]
		res.tab$effect.exposure[pos.change]<-res.tab$effect.exposure[pos.change]*-1
		eff.allele.change<-res.tab$allele1[pos.change]
		oth.allele.change<-res.tab$allele2[pos.change]
		res.tab$allele1[pos.change]<-oth.allele.change
		res.tab$allele2[pos.change]<-eff.allele.change
	
		res.tab$allele1<-toupper(res.tab$allele1) #convert alleles to upper case
		res.tab$allele2<-toupper(res.tab$allele2) #convert alleles to upper case
		res.tab$effect_allele<-toupper(res.tab$effect_allele) #convert alleles to upper case
		res.tab$other_allele<-toupper(res.tab$other_allele) #convert alleles to upper case
		
		strand1<-c("G","C","T","A")
		strand2<-c("C","G","A","T")
		
		#When only the disease effect allele is known, infer other allele from allele1 and allele2 for the exposure 
		for(outcome in res.tab$outcome){
			pos.keep<-which(res.tab$outcome==outcome & is.na(res.tab$other_allele))
			res.tab.test<-res.tab[pos.keep,]
			pos.all<-1:nrow(res.tab)
			res.tab.excl<-res.tab[pos.all[!pos.all %in% pos.keep],]
	
			if(!all(is.na(res.tab.test$effect_allele))){ #this line is necessary because for some studies both effect_allele and other_allele are missing; the script witihn this F statement does not work if all effect_alleles are missing
				res.tab.test$other_allele[which(res.tab.test$effect_allele==res.tab.test$allele1)]<-res.tab.test$allele2[which(res.tab.test$effect_allele==res.tab.test$allele1)]
				res.tab.test$other_allele[which(res.tab.test$effect_allele==res.tab.test$allele2)]<-res.tab.test$allele1[which(res.tab.test$effect_allele==res.tab.test$allele2)]
				#if other allele still missing then the strands are different; recode effect_allele to other strand
				effect_allele_other_strand<-strand2[unlist(lapply(res.tab.test$effect_allele[is.na(res.tab.test$other_allele)],FUN=function(x) which(strand1==x)))]
				res.tab.test$effect_allele[is.na(res.tab.test$other_allele)]<-effect_allele_other_strand
				res.tab.test$other_allele[which(res.tab.test$effect_allele==res.tab.test$allele1)]<-res.tab.test$allele2[which(res.tab.test$effect_allele==res.tab.test$allele1)]
				res.tab.test$other_allele[which(res.tab.test$effect_allele==res.tab.test$allele2)]<-res.tab.test$allele1[which(res.tab.test$effect_allele==res.tab.test$allele2)]
			}
			res.tab<-rbind(res.tab.excl,res.tab.test) #put datasets back together
		}
	#		res.tab.test[,c("outcome","allele1","allele2","effect_allele","other_allele","eaf.outcome","eaf.trait")]
	
		##################
		#Fix disease SNPs#
		##################
		#code effect.outcome so that is per copy of the allele that increases the exposure
		#But be careful that the alleles aren't different because of different strands or because palindromic 
		
	
		######################
		#Fix palindromic SNPs#
		######################
	
		pos.amb<-unlist(lapply(1:4,FUN=function(i) which(with(res.tab, effect_allele==strand1[i] & other_allele==strand2[i]))))
		pos.unamb<-unlist(lapply(1:4,FUN=function(i) which(with(res.tab, effect_allele==strand1[i] & other_allele!=strand2[i]))))
		
		if(eaf.outcome.notmiss.study){ #don't have eaf for some studies
			#For ambiguous/palindromic SNPs, correct effect allele in disease GWAS to be same as effect allele in trait GWAS (allele that increases the trait, using EAF columns to infer the effect allele
			#if the effect allele frequencies are different (eaf.outcome versus eaf.trait) then effect alleles are different
			amb.tab<-res.tab[pos.amb,]
			length(unique(res.tab$snp)) #number of SNPs
			length(unique(amb.tab$snp)) #number of ambiguous/palindromic SNPs
			#exclude ambiguous SNPs where eaf is 0.42-0.58 because the effect allele cannot be reliably inferred
			res.tab[,c("eaf.outcome","eaf.trait")]
			amb.tab.keep<-amb.tab[which(amb.tab$eaf.outcome<0.42 | amb.tab$eaf.outcome>0.58),] #keep ambiguous snps if eaf.outcome NOT 0.42-0.58
			amb.tab.excl<-amb.tab[which(amb.tab$eaf.outcome>0.42 & amb.tab$eaf.outcome<0.58 | is.na(amb.tab$eaf.outcome)),] #exclude ambiguous SNPs if eaf.outcome is 0.42-0.58
	
			length(unique(amb.tab.keep$snp)) #number of ambiguous SNPs kept because eaf.dos NOT 0.42-0.58 and eaf.outcome not missing
			
			length(unique(amb.tab.excl$snp)) #number of ambiguous SNPs excluded because eaf.outcome 0.42-0.58
			amb.tab.keep$eaf.trait[is.na(amb.tab.keep$eaf.trait)]<-amb.tab.keep$eaf.outcome[is.na(amb.tab.keep$eaf.trait)] #if eaf.trait is missing replace with eaf.outcome
			amb.tab.keep<-amb.tab.keep[amb.tab.keep$eaf.trait<0.42 | amb.tab.keep$eaf.trait>0.58,] #keep ambiguous SNPs if eaf.trait NOT 0.42-0.58
			#amb.tab.excl$eaf.trait[is.na(amb.tab.excl$eaf.trait)]<-amb.tab.excl$eaf.outcome[is.na(amb.tab.excl$eaf.trait)] #if eaf.trait is missing replace with eaf.outcome. This doesn't make sense to do 
			amb.tab.excl2<-amb.tab.excl[amb.tab.excl$eaf.trait>0.42 & amb.tab.excl$eaf.trait<0.58,] #exclude ambiguous SNPs if eaf.trait is 0.42-0.58
			amb.tab.excl<-rbind(amb.tab.excl,amb.tab.excl2) #ambigous SNPs excluded because eaf.outcome or eaf.trait is 0.42-0.58
			length(unique(amb.tab.keep$snp)) #number of ambiguous SNPs kept because eaf.outcome or eaf.trait NOT 0.42-0.58
			length(unique(amb.tab.excl$snp)) #number of ambiguous SNPs excluded because eaf.outcome or eaf.trait is 0.42-0.58
			length(unique(amb.tab.excl2$snp)) #number of SNPs excluded because eaf.trait is 0.42-0.58
			
			amb.tab1<-amb.tab.keep[(amb.tab.keep$eaf.outcome>0.5 & amb.tab.keep$eaf.trait>0.5) | (amb.tab.keep$eaf.outcome<0.5 & amb.tab.keep$eaf.trait<0.5),] #amb SNPs where effect alleles same
			amb.tab2<-amb.tab.keep[(amb.tab.keep$eaf.outcome>0.5 & amb.tab.keep$eaf.trait<0.5) | (amb.tab.keep$eaf.outcome<0.5 & amb.tab.keep$eaf.trait>0.5),] #amb SNPs where effect alleles different
			effect.allele<-amb.tab2$allele1
			other.allele<-amb.tab2$allele2
			amb.tab2$effect_allele<-effect.allele
			amb.tab2$other_allele<-other.allele
			amb.tab2$eaf.outcome <-1-amb.tab2$eaf.outcome
			amb.tab2$effect.outcome<-amb.tab2$effect.outcome*-1
		#	amb.tab2[,c("eaf.outcome","eaf.trait")]
		#	amb.tab2[,c("effect_allele","other_allele","allele1","allele2","eaf.outcome","eaf.trait","snp.exposure","effect")]
		}
		
		res.tab.unamb<-res.tab[pos.unamb,] #exclude palindromic SNPs from res.tab
	
		######################
		#non palindromic SNPs#
		######################
		#for non-palindromic SNPs, select SNPs coded using different strands between disease GWAS and trait GWAS
			
		strand.diff<-which(res.tab.unamb$effect_allele != res.tab.unamb$allele1 & res.tab.unamb$effect_allele!=res.tab.unamb$allele2) #get position of non-palindromic SNPs where alleles are coded using different strands
		strand.same<-which(res.tab.unamb$effect_allele == res.tab.unamb$allele1 & res.tab.unamb$other_allele==res.tab.unamb$allele2 
			| res.tab.unamb$effect_allele == res.tab.unamb$allele2 & res.tab.unamb$other_allele ==res.tab.unamb$allele1)   #get position of non-palindromic SNPs where alleles are coded using same strands
	
	#	res.tab.unamb[strand.diff,c("snp.exposure","effect_allele","other_allele","allele1","allele2","eaf.outcome","eaf.trait","effect.outcome","chr")]
	#	res.tab.unamb[strand.same,c("snp.exposure","effect_allele","other_allele","allele1","allele2","eaf.outcome","eaf.trait","effect.outcome","chr")]
	
		#get position of different effect alleles when effects also coded using different strands
	
		#if coded using different strands can use this code to correct
		strdiff.tab<-res.tab.unamb[strand.diff,]
		
		ref.strand<-strdiff.tab$effect_allele
		oth.strand<-strdiff.tab$other_allele
		ref.pos<-unlist(lapply(ref.strand,FUN=function(i) which(strand1==i)))
		oth.pos<-unlist(lapply(oth.strand,FUN=function(i) which(strand1==i)))
	
		#create mirror copies of reference and other allele, ie what they look like on other strand
		reference.allele.other.strand<-unlist(lapply(1:length(ref.pos),FUN=function(i) strand2[ref.pos[i]]))
		other.allele.other.strand<-unlist(lapply(1:length(oth.pos),FUN=function(i) strand2[oth.pos[i]]))
		
		effect.allele<-reference.allele.other.strand
		other.allele<-other.allele.other.strand
		strdiff.tab$effect_allele<-effect.allele
		strdiff.tab$other_allele<-other.allele
		#find SNPs where effect allele in disease GWAS different from effect allele in trait GWAS
		pos.eadiff<-strdiff.tab$effect_allele!=strdiff.tab$allele1 & strdiff.tab$effect_allele==strdiff.tab$allele2 #effect alleles different
		pos.easame<-strdiff.tab$effect_allele==strdiff.tab$allele1 & strdiff.tab$effect_allele!=strdiff.tab$allele2 #effect alleles same 
		ea.diff.tab1<-strdiff.tab[pos.eadiff ,]
		effect.allele<-ea.diff.tab1$allele1
		other.allele<-ea.diff.tab1$allele2
		ea.diff.tab1$effect_allele<-effect.allele #recode effect allele
		ea.diff.tab1$other_allele<-other.allele #recode other allele
		if(eaf.outcome.notmiss.study){ #don't have eaf for studies listed in database 
			ea.diff.tab1$eaf.outcome <-1-ea.diff.tab1$eaf.outcome #recode eaf
		}
		ea.diff.tab1$effect.outcome<-ea.diff.tab1$effect.outcome*-1 #recode odds ratio
		#ea.diff.tab1[,c("eaf.outcome","eaf.trait")]
		
		ea.same.tab1<-strdiff.tab[pos.easame ,]
				
		#exclude non palindromic SNPs where alleles coded using different strands
		res.tab.unamb.same.strand<-res.tab.unamb[strand.same,]
		
		#for non-palindromic SNPs, where strand is same between two databases, recode disease GWAS effect allele to the trait effect allele
		pos.eadiff<-which(with(res.tab.unamb.same.strand, effect_allele != allele1 & effect_allele==allele2))
		pos.easame<-which(with(res.tab.unamb.same.strand, effect_allele == allele1 & effect_allele!=allele2))
		
	#		res.tab.unamb.same.strand[,c("allele1","allele2","effect_allele","other_allele","eaf.trait")]
	#		res.tab[,c("allele1","allele2","effect_allele","other_allele","eaf.trait")]
	
		ea.diff.tab2<-res.tab.unamb.same.strand[pos.eadiff,]	
		effect.allele<-ea.diff.tab2$allele1
		other.allele<-ea.diff.tab2$allele2
		ea.diff.tab2$effect_allele<-effect.allele #recode effect allele
		ea.diff.tab2$other_allele<-other.allele #recode other allele
		if(eaf.outcome.notmiss.study){ #don't have eaf for studies listed in database 
			ea.diff.tab2$eaf.outcome <-1-ea.diff.tab2$eaf.outcome #recode eaf
		}
		ea.diff.tab2$effect.outcome<-ea.diff.tab2$effect.outcome*-1 #recode odds ratio
	
		#exclude non palindromic SNPs where coded using same strand and where effect alleles are different	
		res.tab_unamb_same.strand_same.ea<-res.tab.unamb.same.strand[pos.easame,]
	
		#rbind fixed tables to res.tab
		#amb.tab1 #palindromic SNPs where effect alleles were same
		#amb.tab2 #fixed palindromic SNPs where effect alleles different
		#ea.same.tab1 #non-palindromic SNPs where strands different but effect alleles same
		#ea.diff.tab1 #fixed non-palindromic SNPs where strands different and effect alleles different
		#ea.diff.tab2 #fixed non-palindromic SNPs where strands same and effect alleles different
		
		if(eaf.outcome.notmiss.study){ 
			fix.tab<-rbind(res.tab_unamb_same.strand_same.ea,amb.tab1)
			fix.tab<-rbind(fix.tab,amb.tab2)
			fix.tab<-rbind(fix.tab,ea.same.tab1)
			fix.tab<-rbind(fix.tab,ea.diff.tab1)
			fix.tab<-rbind(fix.tab,ea.diff.tab2)
		}else{
			fix.tab<-rbind(res.tab_unamb_same.strand_same.ea,ea.same.tab1)
			fix.tab<-rbind(fix.tab,ea.diff.tab1)
			fix.tab<-rbind(fix.tab,ea.diff.tab2)
		}
			dim(fix.tab)
		
	
		
		if(build.plus.strand.same){ #if the build and coding strand are same, can use this script to ensure effect alleles are the same; handy if eaf is missing for dsease or trait GWAS, e.g. diagram has no eaf but I know it is coded using same ref strand as so youns metabolomic GWAS
			#for non-palindromic SNPs, where strand is same between two databases, recode disease GWAS effect allele to the trait effect allele
			pos.eadiff<-which(with(res.tab, effect_allele != allele1 & effect_allele==allele2))
			pos.easame<-which(with(res.tab, effect_allele == allele1 & effect_allele!=allele2))
			ea.diff.tab2<-res.tab[pos.eadiff,]	
			effect.allele<-ea.diff.tab2$allele1
			other.allele<-ea.diff.tab2$allele2
			ea.diff.tab2$effect_allele<-effect.allele #recode effect allele
			ea.diff.tab2$other_allele<-other.allele #recode other allele
			if(eaf.outcome.notmiss.study){ #don't have eaf for studies listed in database 
				ea.diff.tab2$eaf.outcome <-1-ea.diff.tab2$eaf.outcome #recode eaf
			}
			ea.diff.tab2$effect.outcome<-ea.diff.tab2$effect.outcome*-1 #recode odds ratio
			res.tab_same.ea<-res.tab[pos.easame,]
			fix.tab<-rbind(res.tab_same.ea,ea.diff.tab2)
		}
		dim(fix.tab)	
	
	#		coordinates in cardiogram are b36 (hg18), + strand 
	#		coordinates in diagram are b36 (hg18), forward strand 
	#		coordinates in ilcco are b37 (hg19) forward strand 
	#		coordinates in 1000genomes are hg19/GRCh37
		
		write.table(fix.tab,paste("data/analysis1_",trait,"_",database,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
		fix.tab
	}
}


outcome.database<-c("cardiogram",
		"cardiogramplusc4d",
		"diagram",
		"diagramplusmetabochip",
		"ilcco",
		"ilcco_adeno",
		"ilcco_squam",
		"smoking_ever",
		"smoking_former",
		"smoking_cpd",
		"smoking_onset",
		"giant_bmi",
		"giant_height",
		"giant_whr",
		"giant_wc",
		"giant_hip",
		"AMD_neovascular",
		"AMD_geographic",
		"AMD_Advanced",
		"ALS_disease",
		"magic_homaIR",
		"magic_homaB",
		"magic_fastProIns",
		"magic_fastIns",
		"magic_hba1c","magic_fastGluc","magic_2hrglucose",
		"icbp003591",
		"icbp003590",
		"icbp003588",
		"icbp003589",
		"glgc_hdl_metab_gwas","glgc_ldl_metab_gwas","glgc_tchol_metab_gwas","glgc_trig_metab_gwas","glgc_hdl_metab","glgc_ldl_metab","glgc_tchol_metab","glgc_trig_metab",
		"pgc_adhd",
		"pgc_scz",
		"pgc_mdd","pgc_bip","pgc_aut",
		"igap_stage1",
		 "gcan_anorexia",
		 "ibd_crohns",
		 "ibd_ucolitis","ibd_ichip",
		 "rheumatoid_arthritis",
		 "gefos_LSbmdpooled",
		 "gefos_FNbmdpooled",
		 "ssgac_educ_attain",
		 "ssgac_eduyrs",
		 "ssgac_college",
		 "egg_bw2",
		 "egg_bl",
		 "egg_hc",
		 "egg_obesity",
		 "gwasCatplusMetab_disease",
		 "gwasCatplusMetab_notdisease",
		"alcohol_dependence_dbGAP",
		"bone_ulstrasound_attenuation_dbGAP",
		"diabetic_nephropathy_dbGAP",
		"IBS_dbGAP",
		"ischemic_stroke_dbGAP",
		"multiple_sclerosis_dbGAP",
		"panscan_dbGAP",
		"neuroblastoma_dbGAP",
		"parkinsons_dbGAP",
		"prostate_cancer_dbGAP",
		"sle_dbGAP",
		"ug_cancers_dbGAP",
#		haematological traits
		"HaemGenRBC_MCHC",
		"HaemGenRBC_MCV",
		"HaemGenRBC_RBC", 
		"HaemGenRBC_Hb", 
		"HaemGenRBC_MCH",  
		"HaemGenRBC_PCV",
#		immunobase
		"atd_cooper_immunobase",
		"nar_faraco_immunobase",
		"cel_trynka_immunobase",
		"pbc_liu_immunobase",
		"jia_hinks_immunobase",	
		"pso_tsoi_immunobase",
		"ms_imsgc_immunobase",    
		"ra_eyre_immunobase",
#		t1dbase
		"t1d_bradfield_t1dbase",
#		melanoma
		"melanoma",
#		metabolome 
		"metabolome",
#		GTEx
		"GTEx",
#		ImmuneCellScience
		"ImmuneCellScience",
#		Uric acid
		"Uric_acid",
#		telomere outcomes, non GWAS SNP lookups for telomeres
		"telomere_outcomes")


names.database<-c("cardiogram","cardiogramplusc4d",
		"diagram","diagramplusmetabochip",
		"ilcco","ilcco_adeno","ilcco_squam",
		"smoking_ever","smoking_former","smoking_cpd","smoking_onset",
		"giant_bmi","giant_height","giant_whr","giant_wc","giant_hip",
		"AMD_neovascular","AMD_geographic","AMD_Advanced","ALS_disease",
		"magic_homaIR","magic_homaB","magic_fastProIns","magic_fastIns","magic_hba1c","magic_fastGluc","magic_2hrglucose",
		"icbp003591","icbp003590","icbp003588","icbp003589",
		"glgc_hdl_metab_gwas","glgc_ldl_metab_gwas","glgc_tchol_metab_gwas","glgc_trig_metab_gwas","glgc_hdl_metab","glgc_ldl_metab","glgc_tchol_metab","glgc_trig_metab",
		"pgc_adhd","pgc_scz","pgc_mdd","pgc_bip","pgc_aut",
		"igap_stage1",
		 "gcan_anorexia",
		 "ibd_crohns","ibd_ucolitis","ibd_ichip",
		 "rheumatoid_arthritis",
		 "gefos_LSbmdpooled","gefos_FNbmdpooled",
		 "ssgac_educ_attain","ssgac_eduyrs","ssgac_college",
		 "egg_bw2","egg_bl","egg_hc","egg_obesity",
		 "gwasCatplusMetab_disease","gwasCatplusMetab_notdisease",
		 #dbGAP
		"alcohol_dependence_dbGAP",
		"bone_ulstrasound_attenuation_dbGAP",
		"diabetic_nephropathy_dbGAP",
		"IBS_dbGAP",
		"ischemic_stroke_dbGAP",
		"multiple_sclerosis_dbGAP",
		"panscan_dbGAP",
		"neuroblastoma_dbGAP",
		"parkinsons_dbGAP",
		"prostate_cancer_dbGAP",
		"sle_dbGAP",
		"ug_cancers_dbGAP",
		#haematological traits
		"HaemGenRBC_MCHC",
		"HaemGenRBC_MCV",
		"HaemGenRBC_RBC", 
		"HaemGenRBC_Hb", 
		"HaemGenRBC_MCH",  
		"HaemGenRBC_PCV",
		#immunobase
		"atd_cooper_immunobase",
		"nar_faraco_immunobase",
		"cel_trynka_immunobase",
		"pbc_liu_immunobase",
		"jia_hinks_immunobase",	
		"pso_tsoi_immunobase",
		"ms_imsgc_immunobase",    
		"ra_eyre_immunobase",
		#t1dbase
		"t1d_bradfield_t1dbase",
		#melanoma
		"melanoma",
		#metabolome 
		"metabolome",
		#GTEx
		"GTEx",
		#ImmuneCellScience
		"ImmuneCellScience",
		#Uric acid
		"Uric_acid",
		#telomere outcomes, non GWAS SNP lookups for telomeres
		"telomere_outcomes")
				
	
	gwascat.files<-dir("/projects/MRC-IEU/users/ph14916/GWAS_summary_data/GWAScat")
	Dates<-unlist(strsplit(gwascat.files,spli="_"))[seq(2,length(gwascat.files)*2,2)]
	Dates<-Dates[!is.na(Dates)]
	Date<-sort(Dates)[length(Dates)]
	gwascat.file<-gwascat.files[grep(Date,gwascat.files)]

gwas.files<-c(
	#Cardiogram
	"cardiogram/CARDIoGRAM_GWAS_RESULTS.txt", 
	"cardiogram/cardiogramplusc4d_180814_update_data.txt",
	#Diagram
	"diagram/DIAGRAMv3.2012DEC17.txt",
	"diagram/DIAGRAM.metabochip.txt",
	 #ILCCO
	"ILCCO/ALL/*",
	"ILCCO/ADENO/*",
	"ILCCO/SQUAM/*", 	
	#TAG
	"TAG/tag.evrsmk.tbl",
	"TAG/tag.former.tbl",
	"TAG/tag.cpd.tbl",
	"TAG/tag.logonset.tbl",
	#GIANT
	"GIANT_2014_2015/BMI2015/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq",
	"GIANT_2014_2015/Height2014/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt", 
	"GIANT_2014_2015/Waist2015/GIANT_2015_WHR_COMBINED_AllAncestries.txt",
	"GIANT_2014_2015/Waist2015/GIANT_2015_WC_COMBINED_AllAncestries.txt",
	"GIANT_2014_2015/Waist2015/GIANT_2015_HIP_COMBINED_AllAncestries.txt",
	#AMD
	"AMDGC/AMDGene2013_Neovascular_v_Controls.txt",
	"AMDGC/AMDGene2013_GeographicAtropy_v_Controls.txt",
	"AMDGC/AMDGene2013_Advanced_v_Controls.txt",
	#ALS
	"ALS/EURO_ALS_SUMMARY_STAT_Al_B",
	#MAGIC
		"magic/MAGIC_ln_HOMA-IR.txt",	
		"magic/MAGIC_ln_HOMA-B.txt",
		"magic/MAGIC_ln_fastingProinsulin.txt",
		"magic/MAGIC_ln_FastingInsulin.txt",
		"magic/MAGIC_HbA1C.txt",
		"magic/MAGIC_FastingGlucose.txt",
		"magic/MAGIC_2hrGlucose_AdjustedForBMI.txt",
	#ICBP
	"icbp/phs000585.pha003591.txt",
	"icbp/phs000585.pha003590.txt",
	"icbp/phs000585.pha003588.txt",
	"icbp/phs000585.pha003589.txt",
	#GLGC
	"glgc/jointGwasMc_HDL.txt",
	"glgc/jointGwasMc_LDL.txt",
	"glgc/jointGwasMc_TC.txt",
	"glgc/jointGwasMc_TG.txt",
	"glgc/Mc_HDL.txt",
	"glgc/Mc_LDL.txt",
	"glgc/Mc_TC.txt",
	"glgc/Mc_TG.txt",
	#PGC
	"PGC/pgc.adhd.full.2012-10.txt",
	"PGC/ckqny.scz2snpres",
	"PGC/pgc.mdd.full.2012-04.txt",
	"PGC/pgc.bip.full.2012-04.txt",
	"PGC/pgc.cross.AUT8.2013-05.txt",
	#IGAP
	"IGAP/IGAP_stage_1.txt",
	#GCAN
	"gcan/gcan_meta.out",
	#IBD
	"IBD/cd-meta.txt",
	"IBD/ucmeta-sumstats.txt",
	"IBD/gwas_ichip_meta_release.txt",
	#RA
	"RA/RA_GWASmeta2_20090505-results.txt",
	#GEFOS
	"GEFOS/GEFOS2_LSBMD_POOLED_GC.txt",
	"GEFOS/GEFOS2_FNBMD_POOLED_GC.txt",
	#SSGAC
	"SSGAC/MA_EA_1st_stage.txt",
	"SSGAC/SSGAC_EduYears_Rietveld2013_publicrelease.txt",
	"SSGAC/SSGAC_College_Rietveld2013_publicrelease.txt",
	#EGG
	"EGG/EGG_BW2_DISCOVERY.txt",
	"EGG/EGG-GWAS-BL.txt",
	"EGG/EGG_HC_DISCOVERY.txt",
	"EGG/EGG_Obesity_Meta_Analysis_1.txt",
	#GWAScatalog
	paste("GWAScat/",gwascat.file,sep=""), #disease
	paste("GWAScat/",gwascat.file,sep=""), #not disease
	#dbGAP collections
	"dbGAP_collections/alcohol_dependence_pha002907/phs000092.pha002907.txt", #alcohol dependence
	"dbGAP_collections/boneUltraSoundAttenuation/phs000007.pha001753.txt", #bone ulstrasound attenuation
	"dbGAP_collections/diabetic_nephropathy/phs000018.pha002866.txt", #diabetic nephropathy
	"dbGAP_collections/IBS/phs000130.pha002847.txt", #IBS
	"dbGAP_collections/ischaemic_stroke/phs000102.pha002844.txt", #ischemic stroke
	"dbGAP_collections/Multiple_sclerosis/phs000171.pha002861.txt", #multiple sclerosis
	"dbGAP_collections/PANSCAN/phs000206.pha002874.txt", #panscan
	"dbGAP_collections/Neuroblastoma/phs000124.pha002895.txt", #neuroblastoma
	"dbGAP_collections/Parkinsons/phs000089.pha002868.txt", #parkinsons
	"dbGAP_collections/Prostate_cancer/phs000207.pha002877.txt", #prostate cancer
	"dbGAP_collections/Systemic_lupus_erythematosus/phs000122.pha002848.txt", #Systemic_lupus_erythematosus
	"dbGAP_collections/Upper_gastrointestinal_cancers/phs000361.pha003106.txt", #Upper_gastrointestinal_cancers
	#Haematological traits
	"haematological_traits/HaemGenRBC_MCHC.txt",
	"haematological_traits/HaemGenRBC_MCV.txt",
	"haematological_traits/HaemGenRBC_RBC.txt",
	"haematological_traits/HaemGenRBC_Hb.txt" ,
	"haematological_traits/HaemGenRBC_MCH.txt",
	"haematological_traits/HaemGenRBC_PCV.txt",
	#Immunobase
	"Immunobase/hg18_gwas_ic_atd_cooper_4_18_0", 
	"Immunobase/hg18_gwas_ic_nar_faraco_4_18_0",
	"Immunobase/hg18_gwas_ic_cel_trynka_4_18_0",  
	"Immunobase/hg18_gwas_ic_pbc_liu_4_18_0",
	"Immunobase/hg18_gwas_ic_jia_hinks_UK_4_18_0",
	"Immunobase/hg18_gwas_ic_pso_tsoi_4_18_0",
	"Immunobase/hg18_gwas_ic_ms_imsgc_4_18_0",     
	"Immunobase/hg18_gwas_ic_ra_eyre_4_18_0",
	#t1dbase 
	"T1DBase/hg18_gwas_t1d_bradfield_4_18_0",
	#melanoma
	"melanoma/melanoma_gwas",
	#metabolome
	"metabolome/metabolites_meta/*",
	#GTEx
#	"GTEx/sigonly/*",
	"GTEx/alldata/*",
	#ImmuneCellScience
	"ImmuneCellScience/2-GWASResults/*",
	#Uric acid
	"Uric_acid/tempPONE/METAserumurate_BMIbySNPinter_all.txt",
	#telomere outcomes, non GWAS SNP lookups in multiple outcomes
	"telomere_outcomes/harmonized/*")

#ensembl function
ensembl_get_position<-function(snp){
	library(biomaRt)
	Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
	Attr<-listAttributes(Mart)
	ensembl<-getBM(attributes=c("refsnp_id","chr_name","chrom_start"),filters="snp_filter",values=snp,mart=Mart)
	ensembl<-ensembl[nchar(ensembl$chr_name)<=2,]
	ensembl 
}
