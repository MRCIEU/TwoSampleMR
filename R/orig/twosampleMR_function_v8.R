#source("/home/ph14916/MR/R/twosampleMR_function_v7.R")

twosampleMR<-function(trait,res.tab,databases,genome.build,liftover.dir,out.file) { 
	for(database in databases){
		res.tab.database<-res.tab[res.tab$database==database,]
		for(outcome in unique(res.tab.database$outcome)){
		res.tab.outcome<-res.tab.database[res.tab.database$outcome==outcome,]
		for(exposure in unique(res.tab.outcome$exposure)){
			res.tab.test<-res.tab.outcome[res.tab.outcome$exposure==exposure,]
				#exclude rows with missing standard errors and missing effects
				res.tab.test<-res.tab.test[!is.na(as.numeric(res.tab.test$se.outcome)),]
				res.tab.test<-res.tab.test[!is.na(as.numeric(res.tab.test$effect.outcome)),]

				if(!"results" %in% dir()) dir.create("results")
				if(!"MRresults" %in% dir("results")) dir.create("results/MRresults")
				if(!"Correlation_matrices" %in% dir("results")) dir.create("results/Correlation_matrices")
				if(!"Sensitivity_analyses" %in% dir("results")) dir.create("results/Sensitivity_analyses")
				if(!"log" %in% dir()) dir.create("log")
				
				if(all(is.na(as.numeric(res.tab.test$effect.outcome)))) {
					write.table(paste("MR not possible with the ",database," database & ",outcome," outcome & ",exposure," exposure because effect.outcome is missing for all SNPs",sep=""),
						paste("log/",exposure,"_",database,"_",outcome,"_",out.file,"_log.txt",sep=""),col.names=F,row.names=F,quote=F)
					cat("MR not possible with database(",database,") & outcome(",outcome,"), & exposure(",exposure,") because effect.outcome is missing for all SNPs \n")
					next
				}
				cat("database=",database," & outcome= ",outcome," & exposure=",exposure,"\n")
				cat(dim(res.tab.test),"\n")
				outcome<-gsub(" ","_",outcome)
				if(outcome == trait | outcome == "telomere_length" & trait=="telomeres") next
				if(any(duplicated(res.tab.test$snp))){ #when GWAScatalog is outcome dataset, get duplicate SNP-outcome assocations
					if(sum(grep("gwasCat",database,invert=T))==1) stop("duplicate SNPs present in non-gwasCat database") 
					snp.dup<-unique(res.tab.test$snp[duplicated(res.tab.test$snp)])
					res.tab.test.notdup<-res.tab.test[!res.tab.test$snp %in% snp.dup,] #create a table of results with non duplicated SNPs
					for(i in 1:length(snp.dup)){
						res.tab.test.dup<-res.tab.test[res.tab.test$snp %in% snp.dup[i],] #table of duplicated results
						res.tab.test.dup.keep<-res.tab.test.dup[as.numeric(res.tab.test.dup$N_total.outcome)==max(as.numeric(res.tab.test.dup$N_total.outcome)),] #choose the result with the largests combined discovery and replication sample size
						res.tab.notdup<-rbind(res.tab.test.notdup,res.tab.test.dup.keep) #join the cleaned results with the other unduplicated results
					}
					res.tab.test<-res.tab.notdup
				}
				
				#individual SNP estimates
				effect.allele<-res.tab.test$effect_allele
				other.allele<-res.tab.test$other_allele
				eaf.exposure<-res.tab.test$eaf.trait
				eaf.outcome<-res.tab.test$eaf.outcome
#				Genes<-res.tab.test$genes
				chr.names<-res.tab.test$chr_name
				chr.pos<-res.tab.test$chr_pos
				snps<-res.tab.test$snp
				trait.test<-rep(exposure,length(snps))
				outcome.var<-rep(unique(outcome),length(snps))
				ratio <- as.numeric(res.tab.test$effect.outcome)/as.numeric(res.tab.test$effect.exposure)
				segd<-as.numeric(res.tab.test$se.outcome)
				gp<-as.numeric(res.tab.test$effect.exposure)
				gd<-as.numeric(res.tab.test$effect.outcome)
				segp<-as.numeric(res.tab.test$se.exposure)
				Cov<-0 #required when estimated in same participants. Because the sign of the term is negative, it implies that variance will be smaller for estimates that come from same study participants
				ratio.se<-sqrt((segd^2/gp^2) + (gd^2/gp^4)*segp^2 - 2*(gd/gp^3)*Cov)
				z.snp<-ratio/ratio.se
				p.snp<-2*pnorm(abs(z.snp) ,lower.tail=F)
				names.snp.matrix<-c("snps","trait.test","effect.allele","other.allele","eaf.exposure","eaf.outcome","chr.names","chr.pos","outcome.var","ratio","ratio.se","z.snp","p.snp","gd","segd","gp","segp")
				vars.snp.matrix<-unlist(lapply(1:length(names.snp.matrix),FUN=function(x) eval(parse(text=names.snp.matrix[x]))))
				results.snps<-data.frame(matrix(vars.snp.matrix,nrow=length(snps),ncol=length(names.snp.matrix)),stringsAsFactors=F)
				names(results.snps)<-c("effect.estimate","trait.test","effect_allele","other_allele","eaf_exposure","eaf_outcome","chr.names","chr.pos","outcome","effect.outcome","se.outcome","z","p","gd","segd","gp","segp")
				results.snps$n.snps<-1
				results.snps$snps.test<-results.snps$effect.estimate
				results.snps$OR<-exp(as.numeric(results.snps$effect.outcome))
				results.snps$LCI<-exp(as.numeric(results.snps$effect.outcome)-as.numeric(results.snps$se.outcome)*1.96)
				results.snps$UCI<-exp(as.numeric(results.snps$effect.outcome)+as.numeric(results.snps$se.outcome)*1.96)
				results.snps$lnlci<-as.numeric(results.snps$effect.outcome)-as.numeric(results.snps$se.outcome)*1.96
				results.snps$lnuci<-as.numeric(results.snps$effect.outcome)+as.numeric(results.snps$se.outcome)*1.96
				results.snps$database<-database
				results.snps$outcome<-outcome
								
				if(nrow(res.tab.test)==1){
					results.snps$effect.estimate<-"single_SNP_available"
					write.table(results.snps,paste("results/MRresults/MR_results_",exposure,"_",database,"_",outcome,"_singleSNP_",out.file,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
					next
				}			
				
				###############################
				#create SNP correlation matrix#
				###############################
				
				if(nrow(res.tab.test)>1){
					
					#load the ld.gen file and then obtain SNP rs ids for the snps in ld.gen file 
					gen.table<-read.table(paste("data/",trait,".geno.ld",sep=""),sep="\t",head=T,colClasses="character",row.names=NULL) #note that gene.ld file only contains one half of the correlation matrix diagonal and does not contain r2 for snp pairs comprised of identical SNPs
					gen.table<-gen.table[gen.table$R.2!="nan",] 
					Bed<-liftOver_function(data.table=res.tab.test,chr_name="chr_name",chr_pos="chr_pos",snp="snp.exposure",genome.build,liftover.dir) #positional information downloaded from ensembl which is hg38
					res.tab.test<-merge(res.tab.test,Bed,by.x="snp.exposure",by.y="snp")
					gen.table<-merge(gen.table,res.tab.test[,c("snp.exposure","Chr_pos_hg19")],by.x="POS1",by.y="Chr_pos_hg19")
					names(gen.table)[names(gen.table)=="snp.exposure"]<-"snp1"
					gen.table<-merge(gen.table,res.tab.test[,c("snp.exposure","Chr_pos_hg19")],by.x="POS2",by.y="Chr_pos_hg19")
					names(gen.table)[names(gen.table)=="snp.exposure"]<-"snp2"
					snps<-unique(res.tab.test$snp)
					snps<-snps[snps %in% gen.table$snp1]
					if(all(!snps %in% gen.table$snp2)) next
					snps<-snps[snps %in% gen.table$snp2]
					#create correlation matrix
					rho<-create_correlation_matrix_function(res.tab.test=res.tab.test,trait=trait,snps=snps,gen.table=gen.table)
					snps<-rownames(rho)
					res.tab.test<-res.tab.test[res.tab.test$snp.exposure%in% snps,]
					write.table(rho,paste("results/Correlation_matrices/Correlation_matrix_",exposure,"_",database,"_",outcome,"_allSNPs_",out.file,".txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
				
					#####################################################
					#SNP pruning, using threshold of r2 >0.90 / VIF >10 #
					#####################################################
	#				rho2<-rho
	#				rho<-rho2
					#create ped and map files for plink
					write.table(res.tab.test$snp.exposure,"data/snplist",sep="\t",col.names=F,row.names=F,quote=F)
					system(paste("/projects/MRC-IEU/programs/twosampleMR/vcftools_v0.1.13/bin/vcftools --vcf data/",trait,"_genotypes_sorted.vcf --snps data/snplist --keep data/eur.1k.txt --max-alleles 2 --min-alleles 2 --remove-indels --plink  --out data/",trait,sep="")) 
					system(paste("plink --file data/",trait," --indep 50 5 10 --noweb --out data/",trait,sep="")) #The parameters for --indep are: window size in SNPs (e.g. 50), the number of SNPs to shift the window at each step (e.g. 5), the VIF threshold. The VIF is 1/(1-R^2)
					snp.incl<-readLines(paste("data/",trait,".prune.in",sep=""))
					rho<-rho[rownames(rho) %in% snp.incl,rownames(rho) %in% snp.incl]
					res.tab.test<-res.tab.test[res.tab.test$snp.exposure %in% snp.incl,]

#					rho<-prune_snps_function(rho=rho,r=0.9486833,database,outcome,exposure,out.file) #rho = correlation matrix r=0.9486833
					
					if(nrow(res.tab.test)==1){ #all SNP pairs are >0.90, so effectively only one SNP needed
						results.snps2<-results.snps[results.snps$effect.estimate %in% snp.incl,]
						results.snps2$effect.estimate<-"single_independent_SNP_randomly_selected"
						results.snps2<-rbind(results.snps,results.snps2)
						write.table(results.snps2,paste("results/MRresults/MR_results_",exposure,"_",database,"_",outcome,"_singleSNP_",out.file,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
						next
					}
					
					snps<-paste(unique(res.tab.test$snp),collapse=",")
					if(any(!is.na(as.numeric(res.tab.test$effect.outcome)))){
						if(all(is.na(as.numeric(res.tab.test$se.outcome)))) stop(paste("se.outcome, BUT NOT effect.outcome, missing for all SNPs for ",database," database and ",outcome," outcome",sep="")) 
						#loglikelihood function for correlated SNPs
            
						opt<-loglikelihood_corr_function(x=as.numeric(res.tab.test$effect.exposure),sigmax=as.numeric(res.tab.test$se.exposure),
							y=as.numeric(res.tab.test$effect.outcome),sigmay=as.numeric(res.tab.test$se.outcome),rho=rho)
					}

					#likelihood results			
					trait.test<-exposure
					effect.outcome.lh<-opt$par[length(res.tab.test$effect.exposure)+1]
					if(all(opt$hessian==0)) stop(print("Lapack routine dgesv: system is exactly singular: U[1,1] = 0"))
					se.outcome.lh<-sqrt(solve(opt$hessian)[length(res.tab.test$effect.exposure)+1,length(res.tab.test$effect.exposure)+1])
					if(paste(sqrt(solve(opt$hessian)[length(res.tab.test$effect.exposure)+1,length(res.tab.test$effect.exposure)+1]))=="NaN") stop(paste("corlikelihood variance is negative for ",database," database & ",outcome," outcome",sep=""))
					z.lh<-effect.outcome.lh/se.outcome.lh
					p.lh<-2*pnorm(abs(z.lh) ,lower.tail=F)
					Q.lh<-2*opt$value
					Q.df.lh<-length(res.tab.test$effect.exposure)-1
					phet.lh<-round(pchisq(2*opt$value, df=length(res.tab.test$effect.exposure)-1, lower.tail=FALSE),3)
					snps.test.lh<-paste(snps,collapse=",")
					n.snps.lh<-length(unlist(strsplit(snps,split=",")))
					effect.estimate.lh<-"likelihood_corrSNPs" 
				
										
					######################################################################################
					#Fixed effects meta analysis of independent SNPs, include only SNPs in LD equilibrium#
					########################################################################################
					#SNP pruning, when r2 >0.1 / VIF>1.11#
					######################################
				
					write.table(res.tab.test$snp.exposure,"data/snplist",sep="\t",col.names=F,row.names=F,quote=F)
					#create ped and map files for plink
					system(paste("/projects/MRC-IEU/programs/twosampleMR/vcftools_v0.1.13/bin/vcftools --vcf data/",trait,"_genotypes_sorted.vcf --snps data/snplist --keep data/eur.1k.txt --max-alleles 2 --min-alleles 2 --remove-indels --plink  --out data/",trait,sep="")) 
					#plink code for pruning
					system(paste("plink --file data/",trait," --indep 50 5 1.11 --noweb --out data/",trait,sep="")) #The parameters for --indep are: window size in SNPs (e.g. 50), the number of SNPs to shift the window at each step (e.g. 5), the VIF threshold. The VIF is 1/(1-R^2)
					snp.incl<-readLines(paste("data/",trait,".prune.in",sep=""))
					rho<-rho[rownames(rho) %in% snp.incl,rownames(rho) %in% snp.incl]
					#my function doesn't always work so don't use
#					rho<-prune_snps_function(rho=rho,r=0.1,database,outcome,exposure,out.file) # r2 >0.01
	
					if(length(snp.incl)>1){ # if length = 1 then only 1 SNP left after further LD pruning
						res.tab.test<-res.tab.test[res.tab.test$snp %in% snp.incl,]
						snps<-unique(res.tab.test$snp)
						effect_allele<-res.tab.test$effect_allele
						other_allele<-res.tab.test$other_allele
						eaf.exposure<-res.tab.test$eaf.trait
						eaf.outcome<-res.tab.test$eaf.outcome
						
						#fixed and random effect models
						ratio <- as.numeric(res.tab.test$effect.outcome)/as.numeric(res.tab.test$effect.exposure)
						segd<-as.numeric(res.tab.test$se.outcome)
						gp<-as.numeric(res.tab.test$effect.exposure)
						gd<-as.numeric(res.tab.test$effect.outcome)
						segp<-as.numeric(res.tab.test$se.exposure)
						Cov<-0 #required when estimated in same participants. Because the sign of the term is negative, it implies that variance will be smaller for estimates that come from same study participants
						ratio.se<-sqrt((segd^2/gp^2) + (gd^2/gp^4)*segp^2 - 2*(gd/gp^3)*Cov)
						simple.se<-as.numeric(res.tab.test$se.outcome)
						if(any(as.numeric(res.tab.test$se.outcome)<0)) stop(paste("se.outcome is negative for",outcome))
						z.snp<-ratio/ratio.se
						p.snp<-2*pnorm(abs(z.snp) ,lower.tail=F)
						z.snp.simple<-ratio/simple.se
						p.snp.simple<-2*pnorm(abs(z.snp.simple) ,lower.tail=F)
						meta.results <-metagen(ratio,ratio.se,comb.fixed=T,sm="MD")
						meta.results.simple<-metagen(ratio,simple.se,comb.fixed=T,sm="MD")

						#fixed effects resuts delta method for standard error
						effect.outcome.fe<-meta.results$TE.fixed
						se.outcome.fe<-meta.results$seTE.fixed
						z.fe<-meta.results$TE.fixed/meta.results$seTE.fixed
						p.fe<-meta.results$pval.fixed
						Q.fe<-meta.results$Q 
						Q.df.fe<-meta.results$df.Q
						phet.fe<-round(pchisq(meta.results$Q, df=meta.results$df.Q, lower.tail=FALSE),3)
						snps.test.fe<-paste(snps,collapse=",")
						n.snps.fe<-length(snps)
						n.snps<-length(snps)
						effect.estimate.fe<-"fixed_effect_IND_PRUN" #fixed effect with independent SNPs
						
						#random effect results delta method for standard error
						effect.outcome.re<-meta.results$TE.random
						se.outcome.re<-meta.results$seTE.random
						z.re<-meta.results$TE.random/meta.results$seTE.random
						p.re<-meta.results$pval.random
						Q.re<-meta.results$Q 
						Q.df.re<-meta.results$df.Q
						phet.re<-round(pchisq(meta.results$Q, df=meta.results$df.Q, lower.tail=FALSE),3)
						snps.test.re<-paste(snps,collapse=",")
						n.snps.re<-length(snps)
						effect.estimate.re<-"random_effect_IND_PRUN" #random effect with independent SNPs
						
						#fixed effect results simple standard error
						effect.outcome.fe.simpleSE<-meta.results.simple$TE.fixed
						se.outcome.fe.simpleSE<-meta.results.simple$seTE.fixed
						z.fe.simpleSE<-meta.results.simple$TE.fixed/meta.results.simple$seTE.fixed
						p.fe.simpleSE<-meta.results.simple$pval.fixed
						Q.fe.simpleSE<-meta.results.simple$Q 
						Q.df.fe.simpleSE<-meta.results.simple$df.Q
						phet.fe.simpleSE<-round(pchisq(meta.results.simple$Q, df=meta.results.simple$df.Q, lower.tail=FALSE),3)
						snps.test.fe.simpleSE<-paste(snps,collapse=",")
						n.snps.fe.simpleSE<-length(snps)
						effect.estimate.fe.simpleSE<-"fixed_effect_IND_PRUN_simpleSE" #fixed effect with independent SNPs using simple standard errors
						
						#random effect results simple standard error
						effect.outcome.re.simpleSE<-meta.results.simple$TE.random
						se.outcome.re.simpleSE<-meta.results.simple$seTE.random
						z.re.simpleSE<-meta.results.simple$TE.random/meta.results.simple$seTE.random
						p.re.simpleSE<-meta.results.simple$pval.random
						Q.re.simpleSE<-meta.results.simple$Q 
						Q.df.re.simpleSE<-meta.results.simple$df.Q
						phet.re.simpleSE<-round(pchisq(meta.results.simple$Q, df=meta.results.simple$df.Q, lower.tail=FALSE),3)
						snps.test.re.simpleSE<-paste(snps,collapse=",")
						n.snps.re.simpleSE<-length(snps)
						effect.estimate.re.simpleSE<-"random_effect_IND_PRUN_simpleSE" #random effect with independent SNPs using simple standard errors
						
						
						
						#likelihood results using uncorrelated genetic variants		
						x<-as.numeric(res.tab.test$effect.exposure) # genetic associations with risk factor 
						sigmax<-as.numeric(res.tab.test$se.exposure) # standard errors 
						y<-as.numeric(res.tab.test$effect.outcome)  # genetic associations with outcome (log odds ratios)
						sigmay<-as.numeric(res.tab.test$se.outcome) # standard errors
										
						# log-likelihood function
						loglikelihood <- function(param) { 
							return(1/2*sum((x-param[1:length(x)])^2/sigmax^2)+1/2* sum((y-param[length(x)+1]*param[1:length(x)])^2/sigmay^2)) 
						}
						opt = optim(c(x, sum(x*y/sigmay^2)/sum(x^2/sigmay^2)), loglikelihood, hessian=TRUE, control = list(maxit=25000)) # optimization command
									
						#likelihood results for uncorrelated genetic variants			
						effect.outcome.ulh<-opt$par[length(x)+1]
						se.outcome.ulh<-sqrt(solve(opt$hessian)[length(x)+1,length(x)+1])
						if(paste(sqrt(solve(opt$hessian)[length(x)+1,length(x)+1]))=="NaN") stop(paste("uncorlikelihood variance is negative for ",database," database & ",outcome," outcome & ",exposure,"exposure",sep=""))
						z.ulh<-effect.outcome.ulh/se.outcome.ulh
						p.ulh<-2*pnorm(abs(z.ulh) ,lower.tail=F)
						Q.ulh<-2*opt$value
						Q.df.ulh<-length(x)-1
						phet.ulh<-round(pchisq(2*opt$value, df=length(x)-1, lower.tail=FALSE),3)
						snps.test.ulh<-paste(snps,collapse=",")
						n.snps.ulh<-length(unlist(strsplit(snps,split=",")))
						effect.estimate.ulh<-"likelihood_uncorrSNPs" 
						
						
						######################
						#sensitivity analyses#
						######################
	 
						#Egger regression 
						BYG<-as.numeric(res.tab.test$effect.outcome)
						BXG<-as.numeric(res.tab.test$effect.exposure)
						se.BetaYG<-as.numeric(res.tab.test$se.outcome)
						se.BetaXG<-as.numeric(res.tab.test$se.exposure)
						
						if(length(BYG)>2){ # egger regression only possible with more than 2 genetic variants
							snps.test<-paste(snps,collapse=",")				
							egger.res<-summary(lm(BYG~BXG,weights=1/se.BetaYG^2))
									   
							# Inference with correct standard errors
							#intercept egger
							egg.intercept.beta<-egger.res$coef["(Intercept)","Estimate"]
							egg.intercept.se<-egger.res$coef["(Intercept)","Std. Error"]/min(1,egger.res$sigma) #corrected 31 July 2015; sigma is the residual standard error
							DF<-length(BYG)-2
							egg.intercept.LCI<-egg.intercept.beta -1*qt(df=DF, 0.975)*egg.intercept.se
							egg.intercept.UCI<-egg.intercept.beta +1*qt(df=DF, 0.975)*egg.intercept.se
							egg.intercept.t<-abs(egg.intercept.beta/egg.intercept.se)
							egg.intercept.p<-2*(1-pt(abs(egg.intercept.beta/egg.intercept.se),DF))
							
							effect.estimate.egg.int<-"Intercept_egger"
							
							#Slope egger 
							egg.bxg.beta<-egger.res$coef["BXG","Estimate"]
							egg.bxg.se<-egger.res$coef["BXG","Std. Error"]/min(1,egger.res$sigma) #code updated on 31 July 2015; sigma is the residual standard error
							egg.bxg.LCI<-egg.bxg.beta-1*qt(df=DF, 0.975)*egg.bxg.se
							egg.bxg.UCI<-egg.bxg.beta+1*qt(df=DF, 0.975)*egg.bxg.se
							egg.bxg.t<-abs(egg.bxg.beta/egg.bxg.se)
							egg.bxg.p= 2*(1-pt(abs(egg.bxg.beta/egg.bxg.se),DF))
							effect.estimate.egg.bxg<-"Slope_egger" 
							
							#parametric bootstrap for confidence intervals and standard error for Egger slope
							egger.boot<-bootCI_egger_slope(BYG,BXG,se.BetaYG,se.BetaXG)    
							egg.bxg.UCI  = unlist(egger.boot[1])
							egg.bxg.LCI  = unlist(egger.boot[2])
							egg.bxg.se   = unlist(egger.boot[3])
	
							#Egger intercept
							names.results<-c("trait.test","n.snps","egg.intercept.beta","egg.intercept.LCI","egg.intercept.UCI","egg.intercept.se","egg.intercept.t","egg.intercept.p","snps.test","effect.estimate.egg.int") 
							vars<-unlist(lapply(1:length(names.results),FUN=function(x) eval(parse(text=names.results[x]))))
							results.table.eggint<-data.frame(matrix(vars,nrow=length(trait.test),ncol=length(names.results)),stringsAsFactors=F)
							names(results.table.eggint)<-c("trait.test","n.snps","effect.outcome","lnlci","lnuci","se.outcome","t","p","snps.test","effect.estimate")
							
							#Egger slope
							names.results<-c("trait.test","n.snps","egg.bxg.beta","egg.bxg.LCI","egg.bxg.UCI","egg.bxg.se","egg.bxg.t","egg.bxg.p","snps.test","effect.estimate.egg.bxg") 
							vars<-unlist(lapply(1:length(names.results),FUN=function(x) eval(parse(text=names.results[x]))))
							results.table.eggslope<-data.frame(matrix(vars,nrow=length(trait.test),ncol=length(names.results)),stringsAsFactors=F)
							names(results.table.eggslope)<-c("trait.test","n.snps","effect.outcome","lnlci","lnuci","se.outcome","t","p","snps.test","effect.estimate")
						}
																		
						#IVW regression approach from Jack Bowden paper, regressing gene-outcome over gene-exposure, intercept constrained to zero
						ivw.res<-summary(lm(BYG~ -1 + BXG,weights=1/se.BetaYG^2))
						ivw.reg.beta<-ivw.res$coef["BXG","Estimate"]
						DF<-length(BYG)-1
						ivw.reg.se<-ivw.res$coef["BXG","Std. Error"]/ivw.res$sigma #sigma is the residual standard error
						ivw.reg.p<-2*(1-pt(abs(ivw.reg.beta/ivw.reg.se),DF))
						ivw.reg.t<-pt(abs(ivw.reg.beta/ivw.reg.se),DF)
						ivw.reg.LCI<-ivw.reg.beta-1*qt(df=DF, 0.975)*ivw.reg.se
						ivw.reg.UCI<-ivw.reg.beta+1*qt(df=DF, 0.975)*ivw.reg.se
	#					IVW_CI  = ivw.reg.beta + c(-1,1)*qt(df=DF, 0.975)*SE
						effect.estimate.ivwreg<-"IVW_regression" 
						
						#jack bowden weighted median estimator
									
						BetaIV<-as.numeric(res.tab.test$effect.outcome)/as.numeric(res.tab.test$effect.exposure)
						BetaXG<-as.numeric(res.tab.test$effect.exposure)
						BetaYG<-as.numeric(res.tab.test$effect.outcome)
						se.BetaXG<-as.numeric(res.tab.test$se.exposure)
						se.BetaYG<-as.numeric(res.tab.test$se.outcome)
						VBj<-((se.BetaYG)^2)/(BetaXG)^2 + (BetaYG^2)*((se.BetaXG^2))/(BetaXG)^4 #delta approximated variance, gives same variance as my script above
						Med.wt.beta<-weighted.median(BetaIV,w=1/VBj)
						res.boot<-BootCI(BXG=BetaXG,BYG = BetaYG,seXG = se.BetaXG,seYG = se.BetaYG)                    
						Med.wt.CIs<-unlist(res.boot[1])
						Med.wt.LCI<-Med.wt.CIs[1]
						Med.wt.UCI<-Med.wt.CIs[2]
						Med.wt.se<-unlist(res.boot[2])
						Med.wt.t<-unlist(res.boot[3])
						Med.wt.p<-unlist(res.boot[4])
						effect.estimate.wtmed<-"Weighted_median_estimator" 
						
						n.snps<-length(BXG)
						snps.test<-paste(res.tab.test$snp.exposure,collapse=",")						
						
						#collate the various results for SNPs in complete equilibrium 
						#sensitivity analyses
						#weighted median estimator
						names.results<-c("trait.test","n.snps","Med.wt.beta","Med.wt.LCI","Med.wt.UCI","Med.wt.se","Med.wt.t","Med.wt.p","snps.test","effect.estimate.wtmed") 
						vars<-unlist(lapply(1:length(names.results),FUN=function(x) eval(parse(text=names.results[x]))))
						results.table.wtmed<-data.frame(matrix(vars,nrow=length(trait.test),ncol=length(names.results)),stringsAsFactors=F)
						names(results.table.wtmed)<-c("trait.test","n.snps","effect.outcome","lnlci","lnuci","se.outcome","t","p","snps.test","effect.estimate")
						
						#IVW regression approach		
						names.results<-c("trait.test","n.snps","ivw.reg.beta","ivw.reg.LCI","ivw.reg.UCI","ivw.reg.se","ivw.reg.t","ivw.reg.p","snps.test","effect.estimate.ivwreg") 
						vars<-unlist(lapply(1:length(names.results),FUN=function(x) eval(parse(text=names.results[x]))))
						results.table.ivwreg<-data.frame(matrix(vars,nrow=length(trait.test),ncol=length(names.results)),stringsAsFactors=F)
						names(results.table.ivwreg)<-c("trait.test","n.snps","effect.outcome","lnlci","lnuci","se.outcome","t","p","snps.test","effect.estimate")
								
						#collate results of sensitivity analyses
						results.all.sens<-list(results.table.ivwreg,results.table.wtmed)
						if("results.table.eggslope" %in% ls()){
							results.all.sens<-list(results.table.ivwreg,results.table.eggslope,results.table.eggint,results.table.wtmed)
						}
						results.all.sens<-do.call(rbind,results.all.sens)
						results.all.sens$outcome<-unique(res.tab.test$outcome)
						results.all.sens$database<-unique(res.tab.test$database)
						results.all.sens$OR<-exp(as.numeric(results.all.sens$effect.outcome))
						results.all.sens$LCI<-exp(as.numeric(results.all.sens$lnlci))
						results.all.sens$UCI<-exp(as.numeric(results.all.sens$lnuci))
						write.table(results.all.sens,paste("results/Sensitivity_analyses/Sensitivity_analyses_",exposure,"_",database,"_",outcome,"_",out.file,".txt",sep=""),sep="\t",col.name=T,row.name=F,quote=F)
					
						#likelihood effect uncorrelated SNPs
						names.results.ulh<-c("trait.test","n.snps.ulh","effect.outcome.ulh","se.outcome.ulh","z.ulh","p.ulh","Q.ulh","Q.df.ulh","phet.ulh","snps.test.ulh","effect.estimate.ulh") 
						vars.ulh<-unlist(lapply(1:length(names.results.ulh),FUN=function(x) eval(parse(text=names.results.ulh[x]))))
						results.table.ulh<-data.frame(matrix(vars.ulh,nrow=length(trait.test),ncol=length(names.results.ulh)),stringsAsFactors=F)
						names(results.table.ulh)<-c("trait.test","n.snps","effect.outcome","se.outcome","z","p","Q","Q.df","phet","snps.test","effect.estimate")
						
						#fixed effect delta method standard error
						names.results.fe<-c("trait.test","n.snps.fe","effect.outcome.fe","se.outcome.fe","z.fe","p.fe","Q.fe","Q.df.fe","phet.fe","snps.test.fe","effect.estimate.fe") 
						vars.fe<-unlist(lapply(1:length(names.results.fe),FUN=function(x) eval(parse(text=names.results.fe[x]))))
						results.table.fe<-data.frame(matrix(vars.fe,nrow=length(trait.test),ncol=length(names.results.fe)),stringsAsFactors=F)
						names(results.table.fe)<-c("trait.test","n.snps","effect.outcome","se.outcome","z","p","Q","Q.df","phet","snps.test","effect.estimate")
						
						#random effect delta method standard error
						names.results.re<-c("trait.test","n.snps.re","effect.outcome.re","se.outcome.re","z.re","p.re","Q.re","Q.df.re","phet.re","snps.test.re","effect.estimate.re") 
						vars.re<-unlist(lapply(1:length(names.results.re),FUN=function(x) eval(parse(text=names.results.re[x]))))
						results.table.re<-data.frame(matrix(vars.re,nrow=length(trait.test),ncol=length(names.results.re)),stringsAsFactors=F)
						names(results.table.re)<-c("trait.test","n.snps","effect.outcome","se.outcome","z","p","Q","Q.df","phet","snps.test","effect.estimate")
								
						#fixed effect simple standard errors
						names.results.fe.simpleSE<-c("trait.test","n.snps.fe.simpleSE","effect.outcome.fe.simpleSE","se.outcome.fe.simpleSE","z.fe.simpleSE","p.fe.simpleSE","Q.fe.simpleSE","Q.df.fe.simpleSE","phet.fe.simpleSE","snps.test.fe.simpleSE","effect.estimate.fe.simpleSE") 
						vars.fe.simpleSE<-unlist(lapply(1:length(names.results.fe.simpleSE),FUN=function(x) eval(parse(text=names.results.fe.simpleSE[x]))))
						results.table.fe.simpleSE<-data.frame(matrix(vars.fe.simpleSE,nrow=length(trait.test),ncol=length(names.results.fe.simpleSE)),stringsAsFactors=F)
						names(results.table.fe.simpleSE)<-c("trait.test","n.snps","effect.outcome","se.outcome","z","p","Q","Q.df","phet","snps.test","effect.estimate")
						
						#random effect simple standard errors
						names.results.re.simpleSE<-c("trait.test","n.snps.re.simpleSE","effect.outcome.re.simpleSE","se.outcome.re.simpleSE","z.re.simpleSE","p.re.simpleSE","Q.re.simpleSE","Q.df.re.simpleSE","phet.re.simpleSE","snps.test.re.simpleSE","effect.estimate.re.simpleSE") 
						vars.re.simpleSE<-unlist(lapply(1:length(names.results.re.simpleSE),FUN=function(x) eval(parse(text=names.results.re.simpleSE[x]))))
						results.table.re.simpleSE<-data.frame(matrix(vars.re.simpleSE,nrow=length(trait.test),ncol=length(names.results.re.simpleSE)),stringsAsFactors=F)
						names(results.table.re.simpleSE)<-c("trait.test","n.snps","effect.outcome","se.outcome","z","p","Q","Q.df","phet","snps.test","effect.estimate")
					}
					
					#results for correlated SNPs
					#likelihood effect correlated SNPs
					names.results.lh<-c("trait.test","n.snps.lh","effect.outcome.lh","se.outcome.lh","z.lh","p.lh","Q.lh","Q.df.lh","phet.lh","snps.test.lh","effect.estimate.lh") 
					vars.lh<-unlist(lapply(1:length(names.results.lh),FUN=function(x) eval(parse(text=names.results.lh[x]))))
					results.table.lh<-data.frame(matrix(vars.lh,nrow=length(trait.test),ncol=length(names.results.lh)),stringsAsFactors=F)
					names(results.table.lh)<-c("trait.test","n.snps","effect.outcome","se.outcome","z","p","Q","Q.df","phet","snps.test","effect.estimate")
										
					results.all.overall<-list(results.table.lh)
					if(any(ls()=="results.table.fe")){	
						results.all.overall<-list(results.table.fe,results.table.re,results.table.fe.simpleSE,results.table.re.simpleSE,results.table.lh,results.table.ulh)
					}
					results.all.overall<-do.call(rbind,results.all.overall)
					results.all.overall$OR<-exp(as.numeric(results.all.overall$effect.outcome))
					results.all.overall$LCI<-exp(as.numeric(results.all.overall$effect.outcome)-1.96*as.numeric(results.all.overall$se.outcome))
					results.all.overall$UCI<-exp(as.numeric(results.all.overall$effect.outcome)+1.96*as.numeric(results.all.overall$se.outcome))
					
					#rbind all results together
					results.all<-rbind.fill(results.all.overall,results.snps)
					if(any(ls()=="results.all.sens")){	
						results.all<-rbind.fill(results.all,results.all.sens)
						if(any(ls()=="egg.intercept.p")) results.all$p_pleiotropy<-egg.intercept.p
						results.all$lnlci[is.na(results.all$lnlci)]<-as.numeric(results.all$effect.outcome[is.na(results.all$lnlci)])-1.96*as.numeric(results.all$se.outcome[is.na(results.all$lnlci)])
						results.all$lnuci[is.na(results.all$lnuci)]<-as.numeric(results.all$effect.outcome[is.na(results.all$lnuci)])+1.96*as.numeric(results.all$se.outcome[is.na(results.all$lnuci)])
					}
					
					results.all$outcome <-outcome
					results.all$database <-database
					if(all(names(results.all) != "lnlci")){
						results.all$lnlci<-as.numeric(results.all$effect.outcome)-1.96*as.numeric(results.all$se.outcome)
						results.all$lnuci<-as.numeric(results.all$effect.outcome)+1.96*as.numeric(results.all$se.outcome)
					}
					
					ens.genes<-ensembl_get_genes(snp=unique(res.tab$snp.exposure))
					results.all<-merge(results.all,ens.genes,by.x="effect.estimate",by.y="snp",all.x=T)
					write.table(results.all,paste("results/MRresults/MR_results_",exposure,"_",database,"_",outcome,"_multiSNPs_",out.file,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
				}
			}
		
		}
	}
}

# Bootstrap method to obtain confidence interval 
# and standard error for slope parameter in 
# MR-Egger regression (with MAF corrected weights)

bootCI_egger_slope(BYG,BXG,se.BetaYG,se.BetaXG,boot = NULL, straps = 10000){
	for (i in 1:straps) {
		BYG_boot = rnorm(length(BYG), mean=BYG, sd=se.BetaYG)
		BXG_boot = rnorm(length(BXG), mean=BXG, sd=se.BetaXG)
		BYG_boot = BYG_boot*sign(BXG_boot)   
		BXG_boot = abs(BXG_boot)      
		boot[i] = summary(lm(BYG_boot~BXG_boot,weights=se.BetaYG^-2))$coef[2,1]
	}
	boot_upper= sort(boot)[9751]
	boot_lower  = sort(boot)[250]
	boot_se    = sd(boot)
	return(list(boot_upper,boot_lower ,boot_se ))
}

#new weighted median estimator function, update 31 July 2015
weighted.median <- function(x=BetaIV, w=1/VBj) {
	N    = length(x)
	ord  = order(x);
	x    = x[ord];
	w    = w[ord];
	Sn   = cumsum(w) 
	S_N  = Sn[N] 
	Pn   = (100/S_N)*(Sn-w/2)
	if(sort(abs(Pn-50))[1] == 0){M = which(Pn==50); return(x[M])} 
	Q  = length(Pn[sign(Pn-50)==-1])
	V1 = Q;  V2 = Q+1
	M  = x[V1] + (50 - Pn[V1])*(x[V2]-x[V1])/(Pn[V2]-Pn[V1])
	return(M) 
}

#new script for weighted median confidence interval, update 31 July 2015
BootCI = function(BXG,BYG,seXG,seYG,Nsim=1000,alpha=0.05,W=1/VBj){
	med = NULL
	for(i in 1:Nsim){
		BYG_boot = rnorm(length(BYG), mean=BYG, sd=seYG)
		BXG_boot = rnorm(length(BXG), mean=BXG, sd=seXG)
		BIVboot  = BYG_boot/BXG_boot
		med[i] = weighted.median(x=BIVboot,w = W)
	}
	lower = Nsim*alpha/2
	upper = Nsim*(1-alpha/2)
	Sort  = sort(med)
	CI    = c(Sort[lower],Sort[upper])
	se    = sd(med)
	t     = mean(med)/se
	p     = 2*(1-pt(abs(t),length(BYG)-1))
	return(list(CI=CI,se=se,t=t,p=p))
}



#ensembl function
ensembl_get_genes<-function(snp){
	library(biomaRt)
	Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
	Attr<-listAttributes(Mart)
	ensembl<-getBM(attributes=c("refsnp_id","associated_gene","ensembl_gene_stable_id"),filters="snp_filter",values=snp,mart=Mart)
	snp.genes<-NULL
	for(snp in ensembl$refsnp_id){
		genes<-paste(unique(unlist(strsplit(paste(ensembl[ensembl$refsnp_id==snp,"associated_gene"],collapse=","),split=","))),collapse=",")
		snp.genes[[snp]]<-matrix(c(snp,genes),nrow=1,ncol=2)
	}
	snp.genes<-data.frame(do.call(rbind,snp.genes))
	names(snp.genes)<-c("snp","genes_ensembl")
	snp.genes
}

#liftover function
liftOver_function<-function(data.table,chr_name,chr_pos,snp,genome.build,liftover.dir){
#	chr_name<-"chr_name"
#	chr_pos<-"chr_pos"
#	snp<-"snp.exposure"
#	genome.build<-"hg38"
	chr<-as.numeric(data.table[,chr_name])
	pos<-as.numeric(data.table[,chr_pos])
	snps<-data.table[,snp] 
	bed.matrix<-data.frame(matrix(c(paste("chr",chr,sep=""),pos,pos+1,snps),nrow=length(snps),ncol=4))
	write.table(bed.matrix,"data/bed.txt",sep=" ",col.names=F,row.names=F,quote=F)
	if(genome.build=="hg38"){
		lift.cmd<-paste(liftover.dir,"liftOver data/bed.txt ",liftover.dir, "hg38ToHg19.over.chain.gz data/output.bed data/unlifted.bed",sep="")
	}
	if(genome.build=="hg18"){
		lift.cmd<-paste(liftover.dir,"liftOver data/bed.txt ",liftover.dir, "hg18ToHg19.over.chain.gz data/output.bed data/unlifted.bed",sep="")
	}
	system(lift.cmd)
	Bed<-read.table("data/output.bed",sep="\t",head=F,stringsAsFactors=F)
	names(Bed)[names(Bed)=="V2"]<-"Chr_pos_hg19" 
	names(Bed)[names(Bed)=="V1"]<-"chr_hg19"
	names(Bed)[names(Bed)=="V4"]<-"snp"
	snps<-paste(Bed$chr_hg19,":",Bed$Chr_pos_hg19,sep="")
	snps<-gsub("chr","",snps)
	Bed$bed_hg19<-snps
	Bed<-Bed[,c("chr_hg19","Chr_pos_hg19" ,"snp","bed_hg19")]
	Bed
}


create_correlation_matrix_function<-function(res.tab.test,trait,snps,gen.table){
	nsnps<-length(snps)
	snps1<-unlist(lapply(1:length(snps),FUN=function(i) rep(snps[i],nsnps)))
	snps2<-rep(snps,nsnps)
	rhosq.tab<-data.frame(snps1,snps2,stringsAsFactors=F)
	rhosq.tab$r2<-NA
	rhosq.tab$r2[rhosq.tab$snps1==rhosq.tab$snps2]<-1
	#which SNP pairs are on different chromosomes? set r2 of such SNPs to 0
	rhosq.tab<-merge(rhosq.tab,res.tab.test[,c("snp.exposure","chr_name")],by.x="snps1",by.y="snp.exposure")
	names(rhosq.tab)[names(rhosq.tab)=="chr_name"]<-"chr1"
	rhosq.tab<-merge(rhosq.tab,res.tab.test[,c("snp.exposure","chr_name")],by.x="snps2",by.y="snp.exposure")
	names(rhosq.tab)[names(rhosq.tab)=="chr_name"]<-"chr2"
	rhosq.tab$r2[rhosq.tab$chr1!=rhosq.tab$chr2]<-0
	
	#obtain R2 stats for rhosq.tab 
	rhosq.tab<-merge(rhosq.tab,gen.table[,c("R.2","snp1","snp2")],by.x=c("snps1","snps2"),by.y=c("snp1","snp2"),all.x=T) #first diagonal of matrix
	rhosq.tab<-merge(rhosq.tab,gen.table[,c("R.2","snp1","snp2")],by.x=c("snps1","snps2"),by.y=c("snp2","snp1"),all.x=T) #other diagonal of matrix
	rhosq.tab$r2[is.na(rhosq.tab$r2)]<-rhosq.tab$R.2.x[is.na(rhosq.tab$r2)]
	rhosq.tab$r2[is.na(rhosq.tab$r2)]<-rhosq.tab$R.2.y[is.na(rhosq.tab$r2)]
#	rhosq.tab<-rhosq.tab[!is.na(rhosq.tab$r2),]
	
	#create correlation matrix
	rho<-sqrt(as.numeric(rhosq.tab$r2)) #Stephen says it should be the sqrt of the R2
	nsnps<-length(unique(rhosq.tab$snps1))
	rho<-matrix(rho,ncol=nsnps,nrow=nsnps)
	colNames<-unique(as.character(rhosq.tab$snps1))
	rowNames<-unique(as.character(rhosq.tab$snps1))
	rownames(rho)<-rowNames
	colnames(rho)<-colNames
	rho
}

loglikelihood_corr_function<-function(x,sigmax,y,sigmay,rho){
	Sigmax = array(NA, dim=c(length(x),length(x))); 
	Sigmay = array(NA, dim=c(length(x),length(x)));
				
	for (k1 in 1:length(x)) { 
		for (k2 in 1:length(x)) { 
			Sigmax[k1,k2] = sigmax[k1]*sigmax[k2]*rho[k1,k2] 
			Sigmay[k1,k2] = sigmay[k1]*sigmay[k2]*rho[k1,k2] 
		}
	}
		
	Taux = solve(Sigmax); Tauy = solve(Sigmay)
		
	loglikelihoodcorrel <- function(param) { # log-likelihood function
		return(1/2*t(x-param[1:length(x)])%*%Taux%*%(x-param[1:length(x)])+1/2* 
			t(y-param[length(x)+1]*param[1:length(x)])%*%Tauy%*% 
				(y-param[length(x)+1]*param[1:length(x)])) }
				
	opt = optim(c(x, sum(x*y/sigmay^2)/sum(x^2/sigmay^2)), 
		loglikelihoodcorrel, hessian=TRUE, control = list(maxit=25000)) 
		# optimization command
	opt 		
}
