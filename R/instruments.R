#' Find instruments for use in MR from the MR Base database
#'
#' This function searches for GWAS significant SNPs (for a given p-value) for a specified set of outcomes.
#' It then performs LD based clumping to return only independent significant associations.
#'
#' @param outcomes Array of outcome IDs (see [`available_outcomes`]).
#' @param p1 Significance threshold. The default is `5e-8`.
#' @param clump Logical; whether to clump results. The default is `TRUE`.
#' @param p2 Secondary clumping threshold. The default is `5e-8`.
#' @param r2 Clumping r2 cut off. The default is `0.001`.
#' @param kb Clumping distance cutoff. The default is `10000`.
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data. The default is `ieugwasr::check_access_token()`.
#' @param force_server Force the analysis to extract results from the server rather than the MRInstruments package.
#'
#' @export
#' @return data frame
extract_instruments <- function(outcomes, p1 = 5e-8, clump = TRUE, p2 = 5e-8, r2 = 0.001, kb = 10000, access_token = ieugwasr::check_access_token(), force_server=FALSE)
{
	# .Deprecated("ieugwasr::tophits()")
	outcomes <- ieugwasr::legacy_ids(unique(outcomes))

	d <- ?ieugwasr::tophits(outcomes, pval=p1, clump=clump, r2=r2, kb=kb, force_server=FALSE, access_token=access_token)

	# d$phenotype.deprecated <- paste0(d$trait, " || ", d$consortium, " || ", d$year, " || ", d$unit)
	if(nrow(d) == 0) return(NULL)
	d$phenotype <- paste0(d$trait, " || id:", d$id)
	d <- format_data(
		d,
		type="exposure",
		snps=NULL,
		phenotype_col="phenotype",
		snp_col="rsid",
		chr_col="chr",
		pos_col="position",
		beta_col="beta",
		se_col="se",
		eaf_col="eaf",
		effect_allele_col="ea",
		other_allele_col="nea",
		pval_col="p",
		samplesize_col="n",
		min_pval=1e-200,
		id_col="id"
	)
	d$data_source.exposure <- "igd"
	return(d)
}

#' Create instruments using the NHGRI-EBI GWAS Catalog
#'
#' This function extracts associations from the GWAS catalog and converts them into the expected exposure format for the harmonise_data function. It does not, however, check that the associations are independent. After running this function, consider using clump_data  to ensure the associations are independent. This function may take a few minutes to run, especially for traits with hundreds of genetic associations identified in multiple studies
#'
#' @param trait the trait of interest as reported in the GWAS catalog
#' @param efo_id identifier for the trait of interest in the experimental factor ontology 
#' @param efo trait of interest in the experimental factor ontology
#' @param map_association_to_study map the genetic association results to study level information in the GWAS catalog, including GWAS catalog study identifiers and reported ancestry. Setting this argument to FALSE will speed up the function but will also result in more biased instruments and biased Mendelian randomization results. Hence it should only be set to FALSE by expert users. Default = TRUE  
#'
#'
#' @importFrom magrittr %>%
#' @export
#' @return data frame

extract_instruments_gwas_catalog<-function(trait=NULL,efo=NULL,efo_id=NULL,map_association_to_study=TRUE)
{

	gwas_associations<-get_gwas_associations(reported_trait=trait,efo_trait=efo,efo_id=efo_id)
# gwas_associations<-gwas_associations_by_reported_trait
	# gwas_associations<-gwasrapidd::union(gwas_associations_by_reported_trait, gwas_associations_by_efo_id)

	if(class(gwas_associations) =="associations") 
	{
		
		if(nrow(gwas_associations@associations)!=0)	
		{
			
			associations<-data.frame(gwas_associations@associations,stringsAsFactors=F)
			risk_alleles<-data.frame(gwas_associations@risk_alleles)
			gwas_results<-merge(associations,risk_alleles,by="association_id")
			col.keep<-c("variant_id","association_id","risk_allele","beta_gc","se_gc","risk_frequency","pvalue","z_scores","z.x","beta_unit","beta_description")
			names_col.keep<-c("rsid","association_id","effect_allele","beta_gc","se_gc","eaf","p","test_statistic","z.x","beta_unit","beta_description")
			
						
			gwas_results$z_scores<-stats::qnorm(gwas_results$pvalue/2,lower.tail=F)
			gwas_results$beta_gc<-log(gwas_results$or_per_copy_number)
			Test<-FALSE 
			if(all(is.na(gwas_results$beta_gc)))
			{
				if(!all(is.na(gwas_results$beta_number))){
					gwas_results$beta_gc<-gwas_results$beta_number		
					Test<-any(is.na(gwas_results$beta_gc))
					if(Test){
						gwas_results2<-gwas_results[is.na(gwas_results$beta_gc),]
						gwas_results2$se_gc<-NA
					}
					gwas_results<-gwas_results[!is.na(gwas_results$beta_gc),]
					Pos1<-which(gwas_results$beta_direction == "increase")
					Pos2<-which(gwas_results$beta_direction == "decrease")
					if(!all(gwas_results$beta_gc>=0)) stop("direction of beta_gc not always positive")
					gwas_results$beta_gc[Pos2]<-gwas_results$beta_gc[Pos2]*-1
					if(!all(unique(gwas_results$beta_direction) %in% c("increase","decrease"))) stop("beta_direction in gwas catalog not always increase or decrease")
				}
			}
			Pos<-which(!is.na(gwas_results$standard_error))		
			gwas_results$se_gc<-NA
			gwas_results$se_gc[Pos]<-gwas_results$standard_error[Pos]
			Pos<-which(is.na(gwas_results$se_gc))
			gwas_results$se_gc[Pos]<-abs(gwas_results$beta_gc[Pos]/gwas_results$z_scores[Pos])
			if(Test){
				gwas_results<-rbind(gwas_results,gwas_results2) 
			}

			# beta_gc still missing. This happens when there is at least one odds ratio amongst the rows. Sometimes the odds ratio colmn is missing beta_number corresponds to signed Z scores. 
			gwas_results$z.x<-gwas_results$beta_gc / gwas_results$se_gc 
			gwas_results1<-gwas_results[!is.na(gwas_results$beta_gc),]
			gwas_results2<-gwas_results[is.na(gwas_results$beta_gc),]
			# update to make allowance for presence of any z scores amongst the beta_units? 
			if(nrow(gwas_results2)>0){
				if(all(!is.na(gwas_results2$beta_unit) & gwas_results2$beta_unit == "z score"))
				{
					gwas_results2$z.x<-gwas_results2$beta_number
					Pos<-which(gwas_results2$beta_direction=="decrease")
					gwas_results2$z.x[Pos]<-gwas_results2$z.x[Pos]*-1
					gwas_results<-rbind(gwas_results1,gwas_results2)
				}
			}
			
			if(!is.null(efo_id)) trait_efo<-efo_id
			if(!is.null(efo)) trait_efo<-efo
			if(!is.null(trait)) trait_efo<-trait
			
			if(is.null(gwas_results))
			{
				warning(paste0("no results found in GWAS catalog for ",trait_efo))
				return("no results found")
			}
			
			if(!is.null(gwas_results))			
			{
				if(nrow(gwas_results)>100) warning("more than 100 genetic associations in the GWAS catalog, which may cause this function to run slowly. To speed things up, you could consider searching on only the reported trait or EFO (not both).")
				# we're only intereted in results where z.x or risk allele frequency are not missing. 
				Pos<-which(!is.na(gwas_results$z.x) |  !is.na(gwas_results$risk_frequency))
				gwas_results<-gwas_results[Pos,]
				
				if(map_association_to_study)
				{
										
					# if(nrow(gwas_associations)>100){
					# 	gwas_associations<-gwas_associations[1:100,]
					# }
					assoc2study<-map_association_to_study_id(associations=gwas_associations)		
					# assoc2study<-association2study
					gwas_results<-merge(gwas_results,assoc2study,by="association_id")	

					ancestry_tab<-make_ancestry_table(association_id=gwas_results$association_id)		
					# ancestry_tab<-anc2

					gwas_results<-merge(gwas_results,ancestry_tab,by="study_id")
					col.keep<-c(col.keep,"study_id","ancestral_group","type")
					names_col.keep<-c(names_col.keep,"study_id","ancestral_group","gc_gwas_type")				
					gwas_results<-gwas_results[,col.keep]
					names(gwas_results)<-c(names_col.keep)				
					gwas_results2<-gwas_results[,c("rsid","beta_description","gc_gwas_type","study_id")]
					gwas_results<-TwoSampleMR::format_data(gwas_results,snp_col="rsid",beta_col="beta_gc",se_col="se_gc",eaf_col="eaf",effect_allele_col="effect_allele",pval_col="p",phenotype_col="study_id",type="exposure",units_col="beta_unit")
					study_id<-unlist(strsplit(gwas_results$exposure,split=" \\("))
					gwas_results$gc_study_id<-study_id[seq(1,length(study_id),by=2)]
					# if(all(gwas_results$gc_study_id=="exposure")) gwas_results$gc_study_id<-NA
					gwas_results$exposure<-paste0(trait_efo," (study: ",gwas_results$gc_study_id,")")
					gwas_results<-merge(gwas_results,gwas_results2,by.x=c("SNP","gc_study_id"),by.y=c("rsid","study_id")) 
					return(gwas_results)
				}	

				if(!map_association_to_study) {
					warning("setting map_association_to_study to FALSE could result in biased instruments and biased MR results.")	
					gwas_results<-gwas_results[,col.keep]
					names(gwas_results)<-c(names_col.keep)
					gwas_results<-TwoSampleMR::format_data(gwas_results,snp_col="rsid",beta_col="beta_gc",se_col="se_gc",eaf_col="eaf",effect_allele_col="effect_allele",pval_col="p",phenotype_col="study_id",type="exposure",units_col="beta_unit")
					study_id<-unlist(strsplit(gwas_results$exposure,split=" \\("))
					gwas_results$gc_study_id<-study_id[seq(1,length(study_id),by=2)]
					if(all(gwas_results$gc_study_id=="exposure")) gwas_results$gc_study_id<-NA
					gwas_results$exposure<-paste0(trait_efo," (study: ",gwas_results$gc_study_id,")")	
					return(gwas_results)
				}
			}
			
		}		
	}

	return(gwas_associations)	
}


map_association_to_study_id<-function(associations=NULL){
	association_ids <- associations@associations$association_id
	    names(association_ids) <- association_ids	
	
	studies <-purrr::map(association_ids, ~ gwasrapidd::get_studies(association_id = .x))	
	
	association2study <-
	      purrr::imap_dfr(
	        studies,
	        ~ tibble::tibble(
	          association_id = .y,
	          study_id = .x@studies$study_id
	        )
	      )

	return(association2study)
}

make_ancestry_table<-function(association_id=NULL){
				
	studies<-gwasrapidd::get_studies(association_id = association_id)
	
	ancestries <-
      dplyr::left_join(studies@ancestries,
                       studies@ancestral_groups,
                       by = c('study_id', 'ancestry_id')) %>%
      dplyr::left_join(studies@countries_of_origin, by = c('study_id', 'ancestry_id')) %>%
      dplyr::rename(
        "co_country_name"="country_name",
        "co_major_area"="major_area",
        "co_region"="region"
      ) %>%
      dplyr::left_join(studies@countries_of_recruitment,
                       by = c('study_id', 'ancestry_id')) %>%      
      dplyr::rename(
        "cr_country_name"="country_name",
        "cr_major_area" ="major_area",
        "cr_region"="region"
      )
	ancestry_tab<-data.frame(ancestries,stringsAsFactors=F)
	ancestry_tab$ancestral_group[ancestry_tab$ancestral_group=="NA"]<-NA
	# if ancestral group is missing replace with major geographic area of recruitment
	ancestry_tab$ancestral_group[is.na(ancestry_tab$ancestral_group)]<-	ancestry_tab$cr_major_area[is.na(ancestry_tab$ancestral_group)]
	# ancestry_tab<-unique(ancestry_tab[,c("study_id","ancestral_group")])
	# ancestry_tab<-ancestry_tab[!is.na(ancestry_tab$ancestral_group),]
	study_ids<-unique(ancestry_tab$study_id)
	
	# for when ancestry varies within study, 
	anc2<-NULL
	for(i in 1:length(study_ids)){
		anc1<-ancestry_tab[ancestry_tab$study_id == study_ids[i],]
		Names<-names(anc1)
		for(j in 1:length(Names)){
			anc1[,Names[j]]<-paste(unique(anc1[,Names[j]]),collapse="; ")	
		}		
		if(nrow(unique(anc1))!=1) stop("unexpected number of rows in ancestry table")
		anc2[[i]]<-unique(anc1)
	}
	anc2<-do.call(rbind,anc2)
	return(anc2)
}


get_gwas_associations<-function(reported_trait=NULL,efo_trait=NULL,efo_id=NULL,verbose = FALSE,warnings = TRUE) 
{
  
	if(!is.null(reported_trait)) 
	{

		
		gwas_studies <- gwasrapidd::get_studies(reported_trait = reported_trait)
		reported_trait<-NULL 
		if(class(unlist(gwas_studies)) != "character")		
		{
			if(nrow(gwas_studies@studies)!=0)
			{	
				gwas_associations_by_reported_trait <- gwasrapidd::get_associations(study_id = gwas_studies@studies$study_id,verbose = verbose,warnings = warnings)
				reported_trait<-"results found"
				# sometimes a study is in the GWAS catalog but has no associations. 
				if(nrow(gwas_associations_by_reported_trait@associations)==0)
				{
					reported_trait<-NULL
				}
			}
		}
	}

	if(!is.null(efo_trait)) 
	{
		gwas_associations_by_efo_trait <- gwasrapidd::get_associations(efo_trait = efo_trait,verbose = verbose,warnings = warnings)
		# association objects are retrieved even if there are no associations in the GWAS catalog 
		if(nrow(gwas_associations_by_efo_trait@associations)==0)
		{
			efo_trait<-NULL
		}
	}

	if(!is.null(efo_id)) 
	{
		gwas_associations_by_efo_id <- gwasrapidd::get_associations(efo_id = efo_id,verbose = verbose,warnings = warnings)
		# sometimes association objects are retrieved even if there are no associations in the GWAS catalog 
		if(nrow(gwas_associations_by_efo_id@associations)==0)
		{
			efo_id<-NULL
		}
	}

	# it is redundant to specify both efo_trait and efo_id. 

	# if(!is.null(reported_trait) && !is.null(efo_trait))
	if(!is.null(reported_trait) && !is.null(efo_trait))
		return(gwasrapidd::union(gwas_associations_by_reported_trait, gwas_associations_by_efo_trait))

	# if(!is.null(reported_trait) && !is.null(efo_id) && is.null(efo_trait))
	if(!is.null(reported_trait) && !is.null(efo_id) && is.null(efo_trait))
		return(gwasrapidd::union(gwas_associations_by_reported_trait, gwas_associations_by_efo_id))

	# if(!is.null(reported_trait) && is.null(efo_trait) && is.null(efo_id))
	if(!is.null(reported_trait) && is.null(efo_trait) && is.null(efo_id))
		return(gwas_associations_by_reported_trait)

	# if(is.null(reported_trait) && !is.null(efo_trait))
	if(is.null(reported_trait) && !is.null(efo_trait))
		return(gwas_associations_by_efo_trait)

	# if(is.null(reported_trait) && is.null(efo_trait) && !is.null(efo_id))
	if(is.null(reported_trait) && is.null(efo_trait) && !is.null(efo_id))
		return(gwas_associations_by_efo_id)

	warning('no results found for the trait/efo in the GWAS catalog. Or perhaps you failed to specify `trait`, `efo` or `efo_id`? At least one of the three must be specified')
	return("no results found")
}
