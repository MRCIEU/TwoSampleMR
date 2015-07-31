find_r2_1k_function<-function(trait,instruments.file,chr_id,chr_pos,snp.exposure,genome.build.hg,liftover.dir,tabix.dir,vcftools.dir,sample.1k,super.pop){
	if(!"data" %in% dir()) dir.create("data")
	in.file<-instruments.file
	fix.tab<-read.table(in.file,sep="\t",head=T,colClasses="character",quote="",fill=T)
	head(fix.tab)
	fix.tab<-fix.tab[fix.tab[,chr_id] != "" ,]
	fix.tab<-fix.tab[!is.na(fix.tab[,chr_id]),]
	chr<-fix.tab[,chr_id][!duplicated(fix.tab[,snp.exposure])]
	chr[which(chr==23)]<-"X"
	pos<-as.numeric(fix.tab[,chr_pos][!duplicated(fix.tab[,snp.exposure])])
	snps<-fix.tab[,snp.exposure][!duplicated(fix.tab[,snp.exposure])]
	
	bed.matrix<-unique(data.frame(matrix(c(paste("chr",chr,sep=""),pos,pos+1,snps),nrow=length(snps),ncol=4),stringsAsFactors=F))
	bed.matrix$X2<-as.numeric(bed.matrix$X2)
	bed.matrix$X3<-as.numeric(bed.matrix$X3)
	write.table(bed.matrix,"data/bed.txt",sep=" ",col.names=F,row.names=F,quote=F)
	if(genome.build.hg=="hg38"){
		lift.cmd<-paste(liftover.dir,"liftOver data/bed.txt ",liftover.dir, "hg38ToHg19.over.chain.gz data/output.bed data/unlifted.bed",sep="")
	}
	if(genome.build.hg=="hg18"){
		lift.cmd<-paste(liftover.dir,"liftOver data/bed.txt ",liftover.dir, "hg18ToHg19.over.chain.gz data/output.bed data/unlifted.bed",sep="")
	}
	system(lift.cmd)
	bed<-read.table("data/output.bed",sep="\t",head=F,stringsAsFactors=F)
	dim(bed) #can happen that get error because numbers sometimes appear as scientific notation e.g. 120000000 -> 1.2e8. liftover doesnt like. may have to use unix commands to fix 
	dim(bed.matrix)
	head(bed)
	fix.tab<-merge(fix.tab,bed,by.x=snp.exposure,by.y="V4")
	names(fix.tab)[names(fix.tab)=="V2"]<-"Chr_pos_hg19"
	fix.tab<-fix.tab[,names(fix.tab)!="V3"]
	fix.tab<-fix.tab[,names(fix.tab)!="V1"]
	
	Region<-unique(paste(fix.tab[,chr_id],":",fix.tab$Chr_pos_hg19,"-",fix.tab$Chr_pos_hg19,sep=""))
	write.table(Region,"data/tabix.regions.txt",col.names=F,row.names=F,quote=F) #I think this line is possibly redundant and undeeded
	chr.index<-unique(fix.tab[,chr_id][1])
	tabix.cmd<-paste(tabix.dir,"tabix -fH ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",chr.index,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > data/",trait,"_genotypes.vcf",sep="")
	system(tabix.cmd)
	pos<-regexpr(":",Region)-1
	chr.test<-substring(Region,1,pos)
	
	for(i in unique(fix.tab[,chr_id])){
		print(i)
		region.test<-Region[chr.test==i]
		for(j in region.test){
			print(j)
			tab.cmd<-paste(tabix.dir,"tabix -f  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",
				i,
				".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ",
				j,
				" >> data/",trait,"_genotypes.vcf ",
				sep="")
			system(tab.cmd)
			system("sleep 1s")
		}
	}
		
	samples.1k<-read.table(sample.1k,sep="\t",head=T,colClasses="character",fill=T) #include only EUROPEAN samples
	samples.1k<-samples.1k[samples.1k$super_pop==super.pop,"sample"] #include only EUROPEAN samples
	#samples.1k<-data.frame(matrix(c(samples.1k,samples.1k),nrow=length(samples.1k)),stringsAsFactors=F)
	write.table(samples.1k,"data/eur.1k.txt",sep="\t",col.names=F,row.names=F,quote=F)
	
	#for specific trait SNPs
	system(paste("bgzip data/",trait,"_genotypes.vcf",sep=""))
	system(paste(tabix.dir,"tabix -p vcf data/",trait,"_genotypes.vcf.gz",sep=""))
	system(paste(vcftools.dir,"bin/vcf-sort data/",trait,"_genotypes.vcf.gz > data/",trait,"_genotypes_sorted.vcf",sep=""))
	system(paste(vcftools.dir,"bin/vcftools --vcf data/",trait,"_genotypes_sorted.vcf  --keep data/eur.1k.txt --max-alleles 2 --min-alleles 2 --remove-indels --geno-r2 --out data/",trait,sep="")) #calculate R2 on same chromosome only
	
	#create geno.ld file with rsids
	geno.ld<-read.table(paste("data/",trait,".geno.ld",sep=""),sep="\t",head=T,colClasses="character")
	out.bed<-read.table("data/output.bed",sep="\t",head=F,colClasses="character")
	chr<-unlist(strsplit(out.bed$V1,"chr"))
	out.bed$chr<-chr[seq(2,length(chr),2)]
	merge.tab<-merge(geno.ld,out.bed[,c("chr","V2","V4")],by.x=c("CHR","POS1"),by.y=c("chr","V2"))
	names(merge.tab)[names(merge.tab)=="V4"]<-"SNP_POS1"
	merge.tab<-merge(merge.tab,out.bed[,c("chr","V2","V4")],by.x=c("CHR","POS2"),by.y=c("chr","V2"))
	names(merge.tab)[names(merge.tab)=="V4"]<-"SNP_POS2"
	write.table(merge.tab,"data/geno.ld_plusRSid.txt",sep="\t",col.names=T,row.names=F,quote=F)
}

#ensembl function
ensembl_get_position<-function(snp){
	library(biomaRt)
	Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
	Attr<-listAttributes(Mart)
	ensembl<-getBM(attributes=c("refsnp_id","chr_name","chrom_start"),filters="snp_filter",values=snp,mart=Mart)
	ensembl<-ensembl[nchar(ensembl$chr_name)<=2,]
	ensembl 
}


#SNPlookups function
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
		pos<-which(names.database==database)
		gwas.file<-paste(gwas.dir,gwas.files[pos],sep="") #file path for GWAS database
		
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
			
		if(database=="metabolome"){
			head.cmd<-paste("head -1 ",gwas.dir,"Shin/M00053.metal.pos.txt.gz.uniform.af.txt > data/filehead.txt",sep="") 
		}		
		
		system(head.cmd)
		
		app.cmd<-"cat data/filehead.txt  data/gwasplusSNPs.txt > data/temp ; mv data/temp  data/gwasplusSNPs.txt"
		system(app.cmd)
		gwasplusSNPs<-read.table("data/gwasplusSNPs.txt",head=T,colClasses="character",fill=T)
		print(nrow(gwasplusSNPs))
	
		#nrows in gwasplusSNPs
		if(nrow(gwasplusSNPs)==0){		
			write.table(paste("the selected trait SNPs are not present in ",database," GWAS",sep=""),
				paste("log/",trait,"_",database,"_SNPlookups_log.txt",sep=""),col.names=F,row.names=F,quote=F)
				cat("the selected exposure SNPs are not present in the",database," database \n")
			next
		}
	
		if(database=="metabolome"){ 			
			tempvar<-unlist(strsplit(gwasplusSNPs$snp,split=":"))
			gwasplusSNPs$metabolonID<-tempvar[seq(1,length(tempvar),2)]
			gwasplusSNPs$snp.outcome<-tempvar[seq(2,length(tempvar),2)]
			pos.trait.start<-regexpr("Shin/",gwasplusSNPs$metabolonID)+5
			pos.trait.end<-regexpr(".metal",gwasplusSNPs$metabolonID)-1
			gwasplusSNPs$outcome<-substr(gwasplusSNPs$metabolonID,pos.trait.start,pos.trait.end)
			gwasplusSNPs<-gwasplusSNPs[,names(gwasplusSNPs)!="snp"]
			names(gwasplusSNPs)[names(gwasplusSNPs)=="snp.outcome"]<-"snp"

		}
			
		if(!database %in% c("ibd_ichip","telomere_outcomes","metabolome","gwasCatplusMetab_notdisease","gwasCatplusMetab_disease","GTEx","ImmuneCellScience")){ # these are databases that have multiple outcome traits
			gwasplusSNPs$outcome<-database
		}	   
		
		gwasplusSNPs$database.outcome<-database
		
		names(gwasplusSNPs)[names(gwasplusSNPs)=="all_freq"]<-"eaf.outcome"
		names(gwasplusSNPs)[names(gwasplusSNPs)=="beta"]<-"effect.outcome"
		names(gwasplusSNPs)[names(gwasplusSNPs)=="se"]<-"se.outcome"


		if(all(!names(gwasplusSNPs) %in% "database.outcome")) stop("database.outcome column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "outcome")) stop("outcome column missing from gwasplusSNPs table")
#		if(all(!names(gwasplusSNPs) %in% "snp.outcome")) stop("snp.outcome column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "eaf.outcome")) stop("eaf.outcome column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "effect_allele")) stop("effect_allele column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "other_allele")) stop("other_allele column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "effect.outcome")) stop("effect.outcome column missing from gwasplusSNPs table")
		if(all(!names(gwasplusSNPs) %in% "se.outcome")) stop("se.outcome column missing from gwasplusSNPs table")
		
		gwasplusSNPs<-gwasplusSNPs[,c("database.outcome","outcome","snp","effect_allele","other_allele","eaf.outcome","effect.outcome","se.outcome")]
		
		eaf.outcome.notmiss.study<-T
		if(all(is.na(gwasplusSNPs$eaf.outcome))){
			eaf.outcome.notmiss.study<-F #don't have eaf for some studies, e.g these studies c("diagram","icbp","pgc_scz","gcan_anorexia","ibd_ucolitis","rheumatoid_arthritis") 
		}
		
		res.tab<-merge(trait.assoc,gwasplusSNPs,by.x="snp.exposure",by.y="snp")
		if(sum(grep("\\.y",names(res.tab)))) stop(paste("some column names are identical in the instruments file and the outcome database; ",names(res.tab)[grep(".y",names(res.tab))]," present in both files",sep=""))
	
		res.tab$effect.exposure<-as.numeric(res.tab$effect.exposure)
		res.tab$effect.outcome<-as.numeric(res.tab$effect.outcome)
		res.tab$eaf.trait<-as.numeric(res.tab$eaf.trait)
		res.tab$eaf.outcome[res.tab$eaf.outcome=="NR"]<-NA
		res.tab$eaf.outcome[res.tab$eaf.outcome=="NR "]<-NA
		res.tab$eaf.outcome<-as.numeric(res.tab$eaf.outcome)
		res.tab<-res.tab[res.tab$eaf.outcome!=1,] # exclude SNPS where eaf is 1 in outcome database
		    
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
		
		#When only the effect allele in outcome database is known, infer other allele from allele1 and allele2 for the exposure 
		if(all(is.na(gwasplusSNPs$other_allele))){
			for(outcome in res.tab$outcome){
				pos.keep<-which(res.tab$outcome==outcome & is.na(res.tab$other_allele))
				res.tab.test<-res.tab[pos.keep,]
				pos.all<-1:nrow(res.tab)
				res.tab.excl<-res.tab[pos.all[!pos.all %in% pos.keep],]
		
				if(!all(is.na(res.tab.test$effect_allele))){ #this line is necessary because for some studies both effect_allele and other_allele are missing; the script witihn this if statement does not get executed if all effect_alleles are missing
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
		ea.diff.tab2$effect.outcome<-ea.diff.tab2$effect.outcome*-1 #recode effect.outcome
	
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

outcome.database<-c("metabolome")
names.database<-c("metabolome")
gwas.files<-c("Shin/*")

	

