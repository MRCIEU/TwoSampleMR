## To Do

## Plot
# categorise traits
# print SNP names or gene names
# option to have no names
# automate volcano plot

## Scan
# Implement cooks distance as option for finding outliers
# Implement MR PRESSO as option for finding outliers

## Scan output
# calculate significant associations (currently in plot function - should remove from here)

## Analysis
# Simulate improvement in multiple testing correction when filtering by outlier associations

testing_commands <- function()
{

	library(devtools)
	load_all()
	library(RadialMR)
	library(dplyr)

	a <- extract_instruments(300)
	b <- extract_outcome_data(a$SNP, 7)
	dat <- harmonise_data(a, b)
	outlierscan <- outlier_scan(dat, mr_method="mr_ivw")
	outlierplot <- outlier_graph(outlierscan)


	urate <- extract_instruments(1055)
	egfr <- extract_outcome_data(a$SNP, 1105)
	urate_egfr <- harmonise_data(urate, egfr)
	outlierscan <- outlier_scan(dat, mr_method="mr_strategy1")
	outlierplot <- outlier_graph(outlierscan)
 

}

#' Outlier scan
#' 
#' A simple wrapper function.
#' Using a summary set, find outliers in the MR analysis between the pair of trais.
#' Find other 'candidate traits' associated with those outliers.
#' Perform MR of each of those candidate traits with the original exposure and outcome
#' 
#' @param dat Output from harmonise_data. Note - only the first id.exposure - id.outcome pair will be used
#' @param outliers Default is to use the RadialMR package to identify IVW outliers. Alternatively can providen an array of SNP names that are present in dat$SNP to use as outliers
#' @param use_proxies Whether to use proxies when looking up associations. FALSE by default for speed
#' @param search_threshold The p-value threshold for detecting an association between an outlier and a candidate trait. Default is 5e-8
#' @param id_list The list of trait IDs to search through for candidate associations. The default is the high priority traits in available_outcomes()
#' @param include_outliers When performing MR of candidate traits against exposures or outcomes, whether to include the original outlier SNP. Default is FALSE.
#' @param mr_method Method to use for candidate trait - exposure/outcome analysis. Default is mr_strategy1. Can also provide basic MR methods e.g. mr_ivw, mr_weighted_mode etc
#'
#' @export
#' @return List
#' dat  Cleaned dat input
#' radialmr  Results from RadialMR analysis
#' outliers  List of outliers used
#' id_list  List of GWAS IDs used
#' search  Result from search of outliers against GWAS IDs
#' candidate_instruments  Instruments for candidate traits
#' candidate_outcome  Extracted instrument SNPs from outcome
#' candidate_outcome_dat  Harmonised candidate - outcome dataset
#' candidate_outcome_mr  MR analysis of candidates against outcome
#' candidate_exposure   Extracted instrument SNPs from exposure
#' candidate_exposure_dat  Harmonised candidate - exposure dataset
#' candidate_exposure_mr  MR analysis of candidates against exposure
outlier_scan <- function(dat, outliers="RadialMR", use_proxies=FALSE, search_threshold=5e-8, id_list="default", include_outliers=FALSE, mr_method="mr_strategy1")
{
	# Get outliers

	output <- list()
	if(length(unique(dat$id.exposure)) > 1 | length(unique(dat$id.outcome)) > 1)
	{
		message("Warning! Multiple exposure/outcome combinations found")
		message("Only using first exposure / outcome combination")
	}
	dat <- subset(dat, id.exposure == id.exposure[1] & id.outcome == id.outcome[1] & mr_keep)
	output$dat <- dat

	stopifnot(length(mr_method) == 1)
	stopifnot(mr_method %in% mr_method_list()$obj | mr_method == "mr_strategy1")

	if(outliers[1] == "RadialMR")
	{
		message("Using RadialMR package to detect outliers")
		cpg <- require(RadialMR)
		if(!cpg)
		{
			stop("Please install the RadialMR package\ndevtools::install_github('WSpiller/RadialMR')")
		}
		cpg <- require(dplyr)
		if(!cpg)
		{
			stop("Please install the RadialMR package\ndevtools::install_github('WSpiller/RadialMR')")
		}
		radial <- RadialMR::RadialMR(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP, "IVW", "YES", "NO", 0.05/nrow(dat), "NO")
		outliers <- as.character(radial$outliers$SNP)
		message("Identified ", length(outliers), " outliers")
		output$radialmr <- radial
		if(length(outliers) == 0)
		{
			message("No outliers found. Exiting")
			return(output)
		}
	} else {
		nout <- length(outliers)
		outliers <- subset(dat, SNP %in% outliers)$SNP
		message(length(outliers), " of ", nout, " of specified outliers present in dat")
		if(length(outliers) == 0)
		{
			message("No outliers found. Exiting")
			return(output)
		}
	}
	output$outliers <- outliers

	# Find associations with outliers

	if(id_list[1] == "default")
	{
		ao <- suppressMessages(available_outcomes())
		ids <- subset(ao, priority == 1 & nsnp > 500000 & sample_size > 5000) %>%
			arrange(desc(sample_size)) %>%
			filter(!duplicated(trait), mr == 1) %>%
			filter(! id %in% c(dat$id.exposure[1], dat$id.outcome[1]))
		id_list <- ids$id
		message("Using default list of ", nrow(ids), " traits")
	}
	output$id_list <- id_list

	output$search <- extract_outcome_data(radial$outliers$SNP, ids$id, proxies=use_proxies)
	output$search$sig <- output$search$pval.outcome < search_threshold
	out2 <- subset(output$search, sig)
	if(nrow(out2) == 0)
	{
		message("Outliers do not associate with any other traits. Try relaxing the search_threshold")
		return(output)
	}
	message("Found ", length(unique(out2$id.outcome)), " candidate traits associated with outliers at p < ", search_threshold)

	message("Finding instruments for candidate traits")
	output$candidate_instruments <- extract_instruments(unique(out2$id.outcome))

	if(nrow(output$candidate_instruments) == 0)
	{
		message("No instruments available for the candidate traits")
		return(output)
	}

	if(!include_outliers)
	{
		message("Removing outlier SNPs from candidate trait instrument lists")
		output$candidate_instruments <- group_by(output$candidate_instruments, id.exposure) %>%
			do({
				x <- .
				y <- subset(out2, id.outcome == x$id.exposure[1])
				x <- subset(x, !SNP %in% y$SNP)
				x
			})
	}

	if(nrow(output$candidate_instruments) == 0)
	{
		message("No instruments available for the candidate traits")
		return(output)
	}

	message(length(unique(output$candidate_instruments$id.exposure)), " traits with at least one instrument")

######

	message("Looking up candidate trait instruments for ", dat$outcome[1])
	output$candidate_outcome <- extract_outcome_data(unique(output$candidate_instruments$SNP), dat$id.outcome[1], proxies=use_proxies)
	if(is.null(output$candidate_outcome))
	{
		message("None of the candidate trait instruments available for ", dat$outcome[1])
		return(output)
	}
	message(nrow(output$candidate_outcome), " instruments extracted for ", dat$outcome[1])
	
	output$candidate_outcome_dat <- suppressMessages(harmonise_data(output$candidate_instruments, output$candidate_outcome))
	output$candidate_outcome_dat <- subset(output$candidate_outcome_dat, mr_keep)
	if(nrow(output$candidate_outcome_dat) == 0)
	{
		message("None of the candidate trait instruments available for ", dat$outcome[1], " after harmonising")
		return(output)
	}

	message("Performing MR of ", length(unique(output$candidate_outcome_dat$id.exposure)), " candidate traits against ", dat$outcome[1])
	if(mr_method == "mr_strategy1")
	{
		temp <- mr_strategy1(output$candidate_outcome_dat)
		output$candidate_outcome_mr <- temp$res
		output$candidate_outcome_mr_full <- temp
	} else {
		output$candidate_outcome_mr <- suppressMessages(mr(output$candidate_outcome_dat, method_list=c("mr_wald_ratio", mr_method)))
	}


######

	message("Looking up candidate trait instruments for ", dat$exposure[1])
	output$candidate_exposure <- extract_outcome_data(unique(output$candidate_instruments$SNP), dat$id.exposure[1], proxies=use_proxies)
	if(is.null(output$candidate_exposure))
	{
		message("None of the candidate trait instruments available for ", dat$exposure[1])
		return(output)
	}
	message(nrow(output$candidate_exposure), " instruments extracted for ", dat$exposure[1])
	
	output$candidate_exposure_dat <- suppressMessages(harmonise_data(output$candidate_instruments, output$candidate_exposure))
	output$candidate_exposure_dat <- subset(output$candidate_exposure_dat, mr_keep)
	if(nrow(output$candidate_exposure_dat) == 0)
	{
		message("None of the candidate trait instruments available for ", dat$exposure[1], " after harmonising")
		return(output)
	}

	message("Performing MR of ", length(unique(output$candidate_exposure_dat$id.exposure)), " candidate traits against ", dat$exposure[1])
	if(mr_method == "mr_strategy1")
	{
		temp <- mr_strategy1(output$candidate_exposure_dat)
		output$candidate_exposure_mr <- temp$res
		output$candidate_exposure_mr_full <- temp
	} else {
		output$candidate_exposure_mr <- suppressMessages(mr(output$candidate_exposure_dat, method_list=c("mr_wald_ratio", mr_method)))
	}

	return(output)
}


#' MR Strategy 1
#' 
#' How to choose the result for a set of different MR analysies?
#' Simple strategy:
#' Use Wald ratio if only one SNP
#' Use IVW if more than one SNP and heterogeneity is low
#' Use weighted mode if more than some minimum number of SNPs and heterogeneity is high
#' 
#' @param dat Output from harmonise_data function
#' @param het_threshold The p-value threshold for Cochran's Q - if lower than this threshold then run weighted mode. Default p = 0.05
#' @param ivw_max_snp Maximum SNPs to allow IVW result even if heterogeneity is high. Default = 1
mr_strategy1 <- function(dat, het_threshold=0.05, ivw_max_snp=1)
{
	message("First pass: running ", length(unique(paste(dat$id.exposure, dat$id.outcome))), " analyses with Wald ratio or IVW")
	a <- suppressMessages(mr(dat, method_list=c("mr_ivw", "mr_wald_ratio")))
	b <- subset(a, method == "Wald ratio")
	message(nrow(b), " analyses can only run Wald ratio")
	c <- subset(a, method == "Inverse variance weighted")
	message(nrow(c), " analyses can run IVW or mode")
	if(nrow(c) > 0)
	{
		d <- suppressMessages(mr_heterogeneity(dat, method_list="mr_ivw"))
		rerun <- subset(d, Q_pval < het_threshold & Q_df >= (ivw_max_snp-1))
		rerun <- paste(rerun$id.exposure, rerun$id.outcome)
		if(length(rerun) < 0)
		{
			message("All eligible IVW results have low heterogeneity")
			return(list(res=a, all=NULL, heterogeneity=d))
		}
		e <- subset(c, !paste(id.exposure, id.outcome) %in% rerun)
		message("Rerunning MR with weighted modal estimator for ", length(rerun), " analysis - this may be quite slow")
		f <- suppressMessages(mr(
			subset(dat, paste(id.exposure, id.outcome) %in% rerun),
			method_list="mr_weighted_mode"
		))
		res <- rbind(b, e, f)
		return(list(res=res, all=a, heterogeneity=d))
	} else {
		return(list(res=a, all=NULL, heterogeneity=NULL))
	}
}

#' Plot results from outlier)scan
#' 
#' Creates a simple graph depicting the connections between outlier instruments, the original exposure and outcome traits, and the detected candidate associations
#' 
#' @param outlierscan Output from outlier_scan function
#' @export
#' @return Prints plot, and returns dataframe of the connections
outlier_graph <- function(outlierscan, mr_threshold_method = "fdr", mr_threshold = 0.05)
{

	a <- require(igraph)
	if(!a)
	{
		stop("Please install the igraph R package")
	}
	a <- require(dplyr)
	if(!a)
	{
		stop("Please install the dplyr R package")
	}


	stopifnot("candidate_outcome_mr" %in% names(outlierscan))
	stopifnot("candidate_exposure_mr" %in% names(outlierscan))


	# Use threshold to retain causal relationships
	stopifnot(length(mr_threshold_method) == 1)
	if(mr_threshold_method != "unadjusted")
	{
		message("Adjusting p-value")
		outlierscan$candidate_outcome_mr$pval_adj <- p.adjust(outlierscan$candidate_outcome_mr$pval, mr_threshold_method)
		outlierscan$candidate_outcome_mr$sig <- outlierscan$candidate_outcome_mr$pval_adj < mr_threshold
		outlierscan$candidate_exposure_mr$pval_adj <- p.adjust(outlierscan$candidate_exposure_mr$pval, mr_threshold_method)
		outlierscan$candidate_exposure_mr$sig <- outlierscan$candidate_exposure_mr$pval_adj < mr_threshold
	} else {
		outlierscan$candidate_outcome_mr$sig <- outlierscan$candidate_outcome_mr$pval < mr_threshold
		outlierscan$candidate_exposure_mr$sig <- outlierscan$candidate_exposure_mr$pval < mr_threshold
	}

	message("Number of candidate - outcome associations: ", sum(outlierscan$candidate_outcome_mr$sig))
	message("Number of candidate - exposure associations: ", sum(outlierscan$candidate_exposure_mr$sig))


	ao <- available_outcomes()
	grp <- rbind(
		# Instruments for exposure
		data.frame(
			from=outlierscan$outliers,
			to=ao$trait[ao$id == outlierscan$dat$id.exposure[1]],
			what="Outlier instruments"
		),
		# Outlier instruments associated with candidate traits
		data.frame(
			from=subset(outlierscan$search, sig)$SNP,
			to=subset(outlierscan$search, sig)$originalname.outcome,
			what="Search associations"
		),
		# candidate traits associated with exposure
		data.frame(
			from=ao$trait[ao$id %in% subset(outlierscan$candidate_exposure_mr, sig)$id.exposure],
			to=ao$trait[ao$id == outlierscan$dat$id.exposure[1]],
			what="Candidate traits to exposure"
		),
		# candidate traits associated with outcome
		data.frame(
			from=ao$trait[ao$id %in% subset(outlierscan$candidate_outcome_mr, sig)$id.exposure],
			to=ao$trait[ao$id == outlierscan$dat$id.outcome[1]],
			what="Candidate traits to outcome"
		),
		data.frame(
			from=ao$trait[ao$id == outlierscan$dat$id.exposure[1]],
			to=ao$trait[ao$id == outlierscan$dat$id.outcome[1]],
			what="Main hypothesis"
		)
	)


	l <- list()
	i <- 1
	temp <- ao$id %in% subset(outlierscan$candidate_exposure_mr, sig)$id.exposure & ao$id %in% subset(outlierscan$candidate_outcome_mr, sig)$id.exposure
	if(sum(temp) > 0)
	{
		l[[i]] <- data.frame(
			name=ao$trait[temp],
			id=ao$id[temp],
			what="both"
		)
		i <- i + 1
	}
	temp <- ao$id %in% subset(outlierscan$candidate_exposure_mr, sig)$id.exposure
	if(sum(temp) > 0)
	{
		l[[i]] <- data.frame(
			name=ao$trait[temp],
			id=ao$id[temp],
			what="exposure"
		)
		i <- i + 1
	}
	temp <- ao$id %in% subset(outlierscan$candidate_outcome_mr, sig)$id.exposure
	if(sum(temp) > 0)
	{
		l[[i]] <- data.frame(
			name=ao$trait[temp],
			id=ao$id[temp],
			what="outcome"
		)
		i <- i + 1
	}
	temp <- ao$id %in% subset(outlierscan$search, sig)$id.outcome
	if(sum(temp) > 0)
	{
		l[[i]] <- data.frame(
			name=ao$trait[temp],
			id=ao$id[temp],
			what="search"
		)
		i <- i + 1
	}
	sig <- bind_rows(l) %>% filter(!duplicated(id))


	nodes <- rbind(
		data_frame(
			name=c(
				ao$trait[ao$id %in% outlierscan$dat$id.exposure[1]],
				ao$trait[ao$id %in% outlierscan$dat$id.outcome[1]]),
			id=c(outlierscan$dat$id.exposure[1], outlierscan$dat$id.outcome[1]),
			what=c("original")
		),
		data_frame(
			name=unique(outlierscan$outliers),
			id=NA,
			what="Outlier instruments"
		),
		sig
	)

	nodes$size <- 15
	nodes$size[nodes$what != "original"] <- 3
	nodes$label <- nodes$name
	nodes$label[! nodes$what %in% c("original")] <- NA
	nodes$colour <- "#386cb0"
	nodes$colour[nodes$what == "original"] <- "white"
	nodes$colour[nodes$what == "Outlier instruments"] <- "#f0027f"


	# Make layout
	layoutd <- group_by(nodes, what) %>%
		do({
			x <- .
			a <- as.data.frame(t(combn(as.character(x$name), 2)))
			# a <- data.frame(from=x$name[1], to=x$name[-1])
			names(a) <- c("from", "to")
			a$what <- x$what[1]
			a
		})

	layoutd2 <- rbind(as.data.frame(layoutd), 
		data.frame(
			from=subset(nodes, what != "original")$name,
			to=sample(subset(nodes, what == "original")$name, length(subset(nodes, what != "original")$name), replace=TRUE),
			what="temp"
		)
	)

	# Different colours for different edge types
	# Original should be large
	# original nodes should be large
	# no labels for confounders

	grp$colour <- "black"
	grp$colour[grp$what == "Outlier instruments"] <- "#7fc97f"
	grp$colour[grp$what == "Search associations"] <- "#beaed4"
	grp$colour[grp$what == "Candidate traits to outcome"] <- "#fdc086"
	grp$colour[grp$what == "Candidate traits to exposure"] <- "#bf5b17"
	grp$size <- 0.1
	grp$size[grp$what == "Main hypothesis"] <- 0.5

	layoutg <- graph_from_data_frame(layoutd2, vertices=nodes)
	l <- layout_with_fr(layoutg)
	grl <- graph_from_data_frame(grp, directed=TRUE, vertices=nodes)
	plot(grl, 
		layout=l, 
		vertex.size=V(grl)$size, 
		vertex.label=V(grl)$label,
		# edge.arrow.size=E(grl)$size, 
		edge.arrow.size=0.3,
		edge.color=E(grl)$colour,
		vertex.label.cex=0.5, 
		vertex.label.family="sans", 
		vertex.label.color="black", 
		vertex.color=V(grl)$colour,
		edge.color="red"
	)
	return(sig=dplyr::as_data_frame(sig))
}


