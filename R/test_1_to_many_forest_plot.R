library(devtools)
install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)


ao<-available_outcomes()

# unique(ao[ao$subcategory == "Anthropometric",c("trait","id","year")])

exp_dat <- extract_instruments(outcomes=c(2,100,1032,104,1,72,999))
table(exp_dat$exposure)
chd_out_dat <- extract_outcome_data(
	snps = exp_dat$SNP,
	outcomes = 7
)

dat <- harmonise_data(
	exposure_dat = exp_dat, 
	outcome_dat = chd_out_dat
)

# sorted by effect size one method per exp-out
res<-mr(dat) 
res$lci<-exp(res$b-1.96*res$se )
res$uci<-exp(res$b+1.96*res$se )
mr_res<-res
mr_res<-split_exposure(mr_res)
min(res$lci)
max(res$uci)
mr_res$ncase<-10000
mr_res<-subset_on_method(mr_res)
mr_res<-Sort.1.to.many(mr_res,b="b",Sort.action=4)

# source("~/TwoSampleMR/R/forest_plot_1-to-many.R")
forest_plot_1_to_many(mr_res,b="b",se="se",
    exponentiate=T,ao_slc=F,Lo=0.3,Up=2.5,
    TraitM="exposure",Col1_title="Risk factor",Col1_width=1.2,by=NULL,
    trans="log2",xlab="OR for CHD per SD increase in risk factor (95% confidence interval)",Addcols=c("nsnp","ncase","pval","id.exposure","method"),Addcol_widths=c(0.5,1.0,1.0,1.0,1.0),Addcol_titles=c("nsnp","ncase","pval","id.exposure","method"))

# Group results by Mr method
res<-mr(dat) 
mr_res<-res
mr_res$ncase<-10000
mr_res<-Sort.1.to.many(mr_res,Group="method",Sort.action=3,Priority="Inverse variance weighted")
mr_res<-split_exposure(mr_res) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 

# source("~/TwoSampleMR/R/forest_plot_1-to-many_v2.R")
forest_plot_1_to_many(mr_res,b="b",se="se",
    exponentiate=T,trans="log2",ao_slc=F,Lo=0.03,
    Up=22,Col1_title="Risk factor",Col1_width=1.7,by="exposure",TraitM="method",
    xlab="OR for CHD per SD increase in risk factor (95% confidence interval)",Addcols=c("nsnp","ncase","pval","id.exposure","method"),Addcol_widths=c(0.5,1.0,1.0,1.0,1.0),Addcol_titles=c("nsnp","ncase","pval","id.exposure","method"))

# stratify on grouping variable

res<-mr(dat) 
mr_res<-res
# Het<-mr_heterogeneity(dat)
# Plt<-mr_pleiotropy_test(dat)
# Sin<-mr_singlesnp(dat)
# res1<-combine_all_mrresults(res,Het,Plt,Sin,ao_slc=T,Exp=F,split.exposure=F,split.outcome=T)

mr_res$lci<-exp(mr_res$b-1.96*mr_res$se )
mr_res$uci<-exp(mr_res$b+1.96*mr_res$se )
# mr_res<-res
mr_res<-split_exposure(mr_res)
mr_res<-subset_on_method(mr_res)
mr_res$subcategory[mr_res$exposure %in% c("Adiponectin","Hip circumference","Waist circumference")]<-"Group 1"
mr_res$subcategory[is.na(mr_res$subcategory)]<-"Group 2"
mr_res<-Sort.1.to.many(mr_res,Sort.action=1,Group="subcategory")
min(mr_res$lci)
max(mr_res$uci)
mr_res$ncase<-10000
# source("~/TwoSampleMR/R/forest_plot_1-to-many_v2.R")
forest_plot_1_to_many(mr_res,b="b",se="se",exponentiate=T,trans="log2",ao_slc=F,Lo=0.3,Up=2.5,TraitM="Trait",Col1_title="Risk factor",Col1_width=1.2,by="subcategory",Addcols=c("nsnp","ncase","pval","id.exposure","method"),Addcol_widths=c(0.5,1.0,1.0,1.0,1.0),Addcol_titles=c("nsnp","ncase","pval","id.exposure","method"))


# res[,c("exposure","method","nsnp","b","se")]

# mr_res<-mr_res[mr_res$method %in% c("Inverse variance weighted","MR Egger"),]
# group by mr method
res<-mr(dat) 
res$lci<-exp(res$b-1.96*res$se )
res$uci<-exp(res$b+1.96*res$se )
min(res$lci)
max(res$uci)
mr_res<-res

mr_res<-Sort.1.to.many(mr_res,Group="method",Sort.action=3,Priority="Inverse variance weighted")
mr_res<-split_exposure(mr_res)

trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

mr_res1<-mr_res[1:10,]
mr_res2<-mr_res[11:22,]
pdf("temp2.pdf")
forest_plot_1_to_many(mr_res2,b="b",se="se",exponentiate=T,trans="log2",ao_slc=F,Lo=0.03,Up=22,width1=1.7,by="exposure",TraitM="method")
dev.off()

p5<-forest_plot_1_to_many(res,b="b",se="se",exponentiate=T,ao_slc=F,Lo=0.3,Up=2.5,TraitM="exposure",width1=1.2,by=NULL,trans="log2")

# mr_res<-res
b="b"
se="se"
exponentiate=T
ao_slc=F
Lo=0.3
Up=2.5
TraitM="exposure"
width1=1.2
by="exposure"
trans="log2"
Alt_name="method"

p5<-forest_plot_1_to_many(mr_res2,b="b",se="se",exponentiate=T,trans="log2",ao_slc=F,Lo=0.5,Up=2.4,width1=1.5,by="Trait",TraitM="Trait",Alt_name="study",trans="log2")
# mr_res<-Sort.1.to.many(mr_res,Group="study",Sort.action=3,Priority="Genetic study")


library(rmarkdown); library(knitr); render("index.rmd")

