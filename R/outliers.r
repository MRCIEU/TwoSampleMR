## To Do

## Plot
# categorise traits
# print SNP names or gene names
# option to have no names

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

	temp <- subset(outlierscan$search, pval.outcome < 5e-8)
	# remove any cholesterol related traits
	temp <- subset(temp, 
		!grepl("simvastatin", outcome, ignore.case = TRUE) &
		!grepl("atorvastatin", outcome, ignore.case = TRUE) &
		!grepl("cholesterol", outcome, ignore.case = TRUE) &
		!grepl("heart disease", outcome, ignore.case = TRUE) &
		!grepl("HDL", outcome, ignore.case = TRUE) &
		!grepl("LDL", outcome, ignore.case = TRUE) &
		!grepl("myocardial", outcome, ignore.case = TRUE) &
		!grepl("ischaemic", outcome, ignore.case = TRUE) &
		!grepl("angina", outcome, ignore.case = TRUE)
	)

temp$outcome[grepl("mass", temp$outcome, ignore.case=TRUE)]
temp$outcome[grepl("ischaemic", temp$outcome, ignore.case=TRUE)]

	outlierscan <- outlier_sig(outlierscan)
	adj <- outlier_adjustment(outlierscan)
	plot_outliers(outlierscan, adj %>% arrange(d) %>% filter(!duplicated(SNP)))


	outlier_network(outlierscan)


	urate <- extract_instruments(1055)
	egfr <- extract_outcome_data(urate$SNP, 1105)
	chd <- extract_outcome_data(urate$SNP, 7)
	urate_egfr <- harmonise_data(urate, egfr)
	urate_chd <- harmonise_data(urate, chd)
	outlierscan <- outlier_scan(urate_egfr, mr_method="mr_ivw")
	outlierscan <- outlier_sig(outlierscan)
	mr_volcano_plot(rbind(outlierscan$candidate_exposure_mr, outlierscan$candidate_outcome_mr))
	outlier_network(outlierscan)
	temp <- outlier_adjustment(outlierscan)

	outlierscan2 <- outlier_scan(urate_chd, mr_method="mr_ivw")
	outlierscan2 <- outlier_sig(outlierscan2)
	tem2 <- outlier_adjustment(outlierscan)

	tem2$d <- (tem2$adj.qi / tem2$adj.Q) / (tem2$qi / tem2$Q)


	res <- mr(dat)
	ind <- dat$beta.exposure < 0
	dat$beta.exposure[ind] <- abs(dat$beta.exposure[ind])
	dat$beta.outcome[ind] <- dat$beta.outcome[ind] * -1
	mr_scatter_plot(res, dat)[[1]] +
	geom_point(data=subset(dat, SNP %in% l$SNP), colour="red") +
	geom_point()


	plot_outliers(outlierscan, adj %>% arrange(d) %>% filter(!duplicated(SNP)))

	# is the outlier adjustment doing the right thing?
	dev.new()
	RadialMR(outlierscan$dat$beta.exposure,outlierscan$dat$beta.outcome,outlierscan$dat$se.exposure,outlierscan$dat$se.outcome,outlierscan$dat$SNP,"BOTH","YES","NO",0.05,"NO")$plot


	a <- extract_instruments(2)
	b <- extract_outcome_data(a$SNP, 7)
	dat <- harmonise_data(a, b)
	outlierscan <- outlier_scan(dat, mr_method="mr_ivw")
	outlierscan <- outlier_sig(outlierscan)
	adj <- outlier_adjustment(outlierscan)
	plot_outliers(outlierscan, adj %>% arrange(d) %>% filter(!duplicated(SNP)))



}


#' Plot volcano plot of many MR analyses
#' 
#' @param res Dataframe similar to output from mr() function, requiring id.exposure, id.outcome. Ideally only provide one MR estimate for each exposure-outcome hypothesis. Also provide a sig column of TRUE/FALSE to determine if that association is to be labelled
#' @param what Whether to plot many exposures against few outcomes (default: exposure) or few exposures against many outcomes (outcome). If e.g. 'exposure' and there are multiple outcomes then will facet by outcome
#' @export
#' @return ggplot of volcano plots
volcano_plot <- function(res, what="exposure")
{
	cpg <- require(ggplot2)
	if(!cpg)
	{
		stop("Please install the ggplot2 package")
	}
	cpg <- require(ggrepel)
	if(!cpg)
	{
		stop("Please install the ggrepel package")
	}

	stopifnot(all(c("outcome", "exposure", "b", "se", "pval") %in% names(res)))
	if(!"sig" %in% names(res))
	{
		warning("Significant associations not defined. Defaulting to p < 0.05 cutoff")
		res$sig <- res$pval < 0.05
	}

	if(what == "exposure")
	{
		form <- as.formula(". ~ outcome")
		nout <- length(unique(res$outcome))
		if(nout > 5)
			warning(nout, " might be too many outcomes to plot clearly")

	} else if(what == "outcome"){
		form <- as.formula(". ~ exposure")
		nout <- length(unique(res$exposure))
		if(nout > 5)
			warning(nout, " might be too many outcomes to plot clearly")
	} else {
		stop("'what' argument should be 'exposure' or 'outcome'")
	}

	ggplot(res, aes(x=b, y=-log10(pval))) +
	geom_vline(xintercept=0, linetype="dotted") +
	geom_errorbarh(aes(xmin=b-1.96*se, xmax=b+1.96*se)) +
	geom_point(aes(colour=sig)) +
	facet_grid(form, scale="free") +
	geom_label_repel(data=subset(res, sig), aes(label=exposure), colour="black", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
	# geom_text_repel(data=subset(res, sig), aes(label=outcome, colour=category)) +
	geom_point(aes(colour=sig)) +
	theme_bw() + 
	theme(panel.grid.major = element_blank(), 
	      panel.grid.minor = element_blank(),
	      strip.text.x=element_text(size=16)
	) +
	# scale_colour_brewer(type="qual", palette="Dark2") +
	# scale_fill_brewer(type="qual", palette="Dark2") +
	labs(x="MR effect size")
}


radialmr <- function(dat, outlier=NULL)
{
	library(ggplot2)
	beta.exposure <- dat$beta.exposure
	beta.outcome <- dat$beta.outcome
	se.outcome <- dat$se.outcome
	w <- sqrt(beta.exposure^2 / se.outcome^2)
	ratio <- beta.outcome / beta.exposure
	ratiow <- ratio*w
	dat <- data.frame(w=w, ratio=ratio, ratiow=ratiow)
	if(is.null(outlier))
	{
		dat2 <- dat
	} else {
		dat2 <- dat[-c(outlier), ]
	}
	mod <- lm(ratiow ~ -1 + w, dat2)$coefficients[1]
	print(mod)
	ggplot(dat, aes(x=w, y=ratiow)) +
	geom_point() +
	geom_abline(slope=mod, intercept=0)
}





#' Outlier adjustment estimation
#' 
#' How much of the heterogeneity due to the outlier can be explained by alternative pathways?
#' 
#' @param outlierscan Output from outlier_scan
#' @export
#' @return data frame of adjusted effect estimates and heterogeneity stats
outlier_adjustment <- function(outlierscan)
{
	# for each outlier find the candidate MR analyses
	# if only exposure then ignore
	# if only outcome then re-estimate the snp-outcome association
	# if exposure and outcome then re-estimate the snp-exposure and snp-outcome association
	# outlier
	# candidate
	# what 
	# old.beta.exposure
	# adj.beta.exposure
	# old.beta.outcome
	# adj.beta.outcome
	# candidate.beta.outcome
	# candidate.se.outcome
	# candidate.beta.exposure
	# candidate.se.exposure
	# old.deviation
	# new.deviation

	l <- list()
	sig <- subset(outlierscan$search, sig)
	sige <- subset(outlierscan$candidate_exposure_mr, sig)
	sigo <- subset(outlierscan$candidate_outcome_mr, sig)
	if(!is.null(outlierscan$exposure_candidate_mr) & is.data.frame(outlierscan$exposure_candidate_mr))
	{
		sigc <- subset(outlierscan$exposure_candidate_mr, sig)		
	} else {
		sigc <- sige[1,]
	}


	dat <- outlierscan$dat
	dat$qi <- cochrans_q(dat$beta.outcome / dat$beta.exposure, dat$se.outcome / abs(dat$beta.exposure))
	print(dat$qi)
	dat$Q <- sum(dat$qi)

	for(i in 1:nrow(sig))
	{
		a <- subset(dat, SNP == sig$SNP[i], select=c(SNP, beta.exposure, beta.outcome, se.exposure, se.outcome, qi, Q))
		if(sig$id.outcome[i] %in% sigo$id.exposure)
		{
			a$candidate <- sig$outcome[i]
			a$i <- i
			if(sig$id.outcome[i] %in% sige$id.exposure & sig$id.outcome[i] %in% sigc$id.exposure)
			{
				message("p->x; x->p: ", a$SNP, " - ", sig$outcome[1])
				message("Bi-directional effects cannot be handled correctly at the moment. Skipping")
				a$what <- "p->x; x->p"
				a$candidate.beta.exposure <- NA
				a$candidate.se.exposure <- NA
				a$candidate.beta.outcome <- NA
				a$candidate.se.outcome <- NA
				a$adj.beta.exposure <- NA
				a$adj.se.exposure <- NA
				a$adj.beta.outcome <- NA
				a$adj.se.outcome <- NA

			} else if(sig$id.outcome[i] %in% sige$id.exposure) {
				message("p->x; p->y: ", a$SNP, " - ", sig$outcome[i])
				a$what <- "p->x; p->y"
				a$candidate.beta.exposure <- sige$b[sige$id.exposure == sig$id.outcome[i]]
				a$candidate.se.exposure <- sige$se[sige$id.exposure == sig$id.outcome[i]]
				a$candidate.beta.outcome <- sigo$b[sigo$id.exposure == sig$id.outcome[i]]
				a$candidate.se.outcome <- sigo$se[sigo$id.exposure == sig$id.outcome[i]]
				b <- bootstrap_path(a$beta.exposure, a$se.exposure, sig$beta.outcome[i], sig$se.outcome[i], sige$b[sige$id.exposure == sig$id.outcome[i]], sige$se[sige$id.exposure == sig$id.outcome[i]])
				a$adj.beta.exposure <- b[1]
				a$adj.se.exposure <- b[2]
				b <- bootstrap_path(a$beta.outcome, a$se.outcome, sig$beta.outcome[i], sig$se.outcome[i], sigo$b[sigo$id.exposure == sig$id.outcome[i]], sigo$se[sigo$id.exposure == sig$id.outcome[i]])
				a$adj.beta.outcome <- b[1]
				a$adj.se.outcome <- b[2]
			} else if(sig$id.outcome[i] %in% sigc$id.exposure) {
				message("x->p; p->y: ", a$SNP, " - ", sig$outcome[i])
				a$what <- "x->p; p->y"
				a$candidate.beta.exposure <- sigc$b[sigc$id.exposure == sig$id.outcome[i]]
				a$candidate.se.exposure <- sigc$se[sigc$id.exposure == sig$id.outcome[i]]
				a$candidate.beta.outcome <- sigo$b[sigo$id.exposure == sig$id.outcome[i]]
				a$candidate.se.outcome <- sigo$se[sigo$id.exposure == sig$id.outcome[i]]
				b <- bootstrap_path(a$beta.exposure, a$se.exposure, sig$beta.outcome[i], sig$se.outcome[i], sigc$b[sigc$id.exposure == sig$id.outcome[i]], sigc$se[sigc$id.exposure == sig$id.outcome[i]])
				a$adj.beta.exposure <- b[1]
				a$adj.se.exposure <- b[2]
				b <- bootstrap_path(a$beta.outcome, a$se.outcome, sig$beta.outcome[i], sig$se.outcome[i], sigo$b[sigo$id.exposure == sig$id.outcome[i]], sigo$se[sigo$id.exposure == sig$id.outcome[i]])
				a$adj.beta.outcome <- b[1]
				a$adj.se.outcome <- b[2]
			} else {
				message("outcome: ", a$SNP, " - ", sig$outcome[i])
				a$what <- "outcome"
				a$candidate.beta.exposure <- NA
				a$candidate.se.exposure <- NA
				a$candidate.beta.outcome <- sigo$b[sigo$id.exposure == sig$id.outcome[i]]
				a$candidate.se.outcome <- sigo$se[sigo$id.exposure == sig$id.outcome[i]]
				b <- bootstrap_path(a$beta.outcome, a$se.outcome, sig$beta.outcome[i], sig$se.outcome[i], sigo$b[sigo$id.exposure == sig$id.outcome[i]], sigo$se[sigo$id.exposure == sig$id.outcome[i]])
				a$adj.beta.exposure <- a$beta.exposure
				a$adj.se.exposure <- a$se.exposure
				a$adj.beta.outcome <- b[1]
				a$adj.se.outcome <- b[2]
			}
			temp <- dat
			temp$beta.exposure[temp$SNP == a$SNP] <- a$adj.beta.exposure
			temp$beta.exposure.se[temp$SNP == a$SNP] <- a$adj.beta.exposure.se
			temp$beta.outcome[temp$SNP == a$SNP] <- a$adj.beta.outcome
			temp$beta.outcome.se[temp$SNP == a$SNP] <- a$adj.beta.outcome.se
			temp$qi <- cochrans_q(temp$beta.outcome / temp$beta.exposure, temp$se.outcome / abs(temp$beta.exposure))
			a$adj.qi <- temp$qi[temp$SNP == a$SNP]
			a$adj.Q <- sum(temp$qi)

			l[[i]] <- a
		}

	}
	l <- bind_rows(l)
	l$d <- (l$adj.qi / l$adj.Q) / (l$qi / l$Q)
	return(l)
}

plot_outliers <- function(outlierscan, adj)
{
	dat <- subset(outlierscan$dat, mr_keep, select=c(SNP, beta.exposure, beta.outcome, se.exposure, se.outcome))
	dat$ratio <- dat$beta.outcome / dat$beta.exposure
	dat$weights <- sqrt(dat$beta.exposure^2 / dat$se.outcome^2)
	dat$what <- "Unadjusted"
	dat$candidate <- "NA"
	temp <- subset(adj, select=c(SNP, adj.beta.exposure, adj.beta.outcome, adj.se.exposure, adj.se.outcome, candidate))
	ind <- grepl("\\|", temp$candidate)
	if(any(ind))
	{
		temp$candidate[ind] <- sapply(strsplit(temp$candidate[ind], split=" \\|"), function(x) x[[1]])
	}
	names(temp) <- gsub("adj.", "", names(temp))
	temp$what <- "Adjusted"
	temp$ratio <- temp$beta.outcome / temp$beta.exposure
	temp$weights <- sqrt(temp$beta.exposure^2 / temp$se.outcome^2)
	# ind <- dat$beta.exposure < 0
	# dat$beta.exposure[ind] <- dat$beta.exposure[ind] * -1
	# dat$beta.outcome[ind] <- dat$beta.outcome[ind] * -1
	# ind <- temp$beta.exposure < 0
	# temp$beta.exposure[ind] <- temp$beta.exposure[ind] * -1
	# temp$beta.outcome[ind] <- temp$beta.outcome[ind] * -1

	# est <- with(subset(dat, !SNP %in% temp$SNP), mr_ivw(beta.exposure, beta.outcome, se.exposure, se.outcome)) %>% as.data.frame

	est <- lm(ratio ~ -1 + weights, data=subset(dat, !SNP %in% temp$SNP))$coefficients[1]

	temp2 <- merge(dat, temp, by="SNP")

	labs <- rbind(
		data.frame(label=temp2$SNP, x=temp2$weights.x, y=temp2$weights.x * temp2$ratio.x),
		data.frame(label=temp$candidate, x=temp$weights, y=temp$weights * temp$ratio)
	)

	p <- ggplot(rbind(dat, temp), aes(y=ratio * weights, x=weights)) +
	geom_label_repel(data=labs, aes(x=x, y=y, label=label), size=2, segment.colour = "grey") +
	geom_point(aes(colour=what)) +
	geom_segment(data=temp2, aes(x=weights.x, xend=weights.y, y=ratio.x * weights.x, yend=ratio.y * weights.y), arrow = arrow(length = unit(0.01, "npc"))) +
	geom_abline(intercept=0, slope=est, colour="red") +
	labs(colour="")
	p
	return(p)
}


cochrans_q <- function(b, se)
{
	xw <- sum(b / se^2) / sum(1/se^2)
	qi <- (1/se^2) * (b - xw)^2
	return(qi)
}

bootstrap_path <- function(gx, gx.se, gp, gp.se, px, px.se, nboot=1000)
{
	res <- rep(0, nboot)
	for(i in 1:nboot)
	{
		res[i] <- rnorm(1, gx, gx.se) - rnorm(1, gp, gp.se) * rnorm(1, px, px.se)
	}
	pe <- gx - gp * px
	return(c(pe, sd(res)))
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
#' @param mr_method Method to use for candidate trait - exposure/outcome analysis. Default is strategy1. Can also provide basic MR methods e.g. mr_ivw, mr_weighted_mode etc
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
outlier_scan <- function(dat, outliers="RadialMR", use_proxies=FALSE, search_threshold=5e-8, id_list="default", include_outliers=FALSE, mr_method="strategy1")
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
	stopifnot(mr_method %in% mr_method_list()$obj | mr_method == "strategy1")

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
	if(mr_method == "strategy1")
	{
		temp <- strategy1(output$candidate_outcome_dat)
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
	if(mr_method == "strategy1")
	{
		temp <- strategy1(output$candidate_exposure_dat)
		output$candidate_exposure_mr <- temp$res
		output$candidate_exposure_mr_full <- temp
	} else {
		output$candidate_exposure_mr <- suppressMessages(mr(output$candidate_exposure_dat, method_list=c("mr_wald_ratio", mr_method)))
	}

#####

	message("Looking up exposure instruments for ", length(unique(out2$id.outcome)), " candidate traits")
	output$exposure_candidate <- extract_outcome_data(unique(output$dat$SNP), unique(out2$id.outcome), proxies=use_proxies)
	if(is.null(output$exposure_candidate))
	{
		message("None of the candidate trait instruments available for ", dat$exposure[1])
		return(output)
	}
	message(nrow(output$candidate_exposure), " instruments extracted")
	
	temp <- subset(output$dat, select=c(SNP, beta.exposure, se.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, id.exposure, exposure, units.exposure))

	if(!include_outliers)
	{
		message("Removing outlier SNPs from candidate trait outcome lists")
		n1 <- nrow(output$exposure_candidate)
		output$exposure_candidate <- group_by(output$exposure_candidate, id.outcome) %>%
			do({
				x <- .
				y <- subset(out2, id.outcome == x$id.outcome[1])
				x <- subset(x, !SNP %in% y$SNP)
				x
			})
		message("Removed ", n1 - nrow(output$exposure_candidate), " outlier SNPs")
	}


	output$exposure_candidate_dat <- suppressMessages(harmonise_data(temp, output$exposure_candidate))
	output$exposure_candidate_dat <- subset(output$exposure_candidate_dat, mr_keep)
	if(nrow(output$exposure_candidate_dat) == 0)
	{
		message("None of the candidate traits have the ", dat$exposure[1], " instruments after harmonising")
		return(output)
	}

	message("Performing MR of ", dat$exposure[1], " against ", length(unique(output$exposure_candidate_dat$id.outcome)), " candidate traits")
	if(mr_method == "strategy1")
	{
		temp <- strategy1(output$exposure_candidate_dat)
		output$exposure_candidate_mr <- temp$res
		output$exposure_candidate_mr_full <- temp
	} else {
		output$exposure_candidate_mr <- suppressMessages(mr(output$exposure_candidate_dat, method_list=c("mr_wald_ratio", mr_method)))
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
strategy1 <- function(dat, het_threshold=0.05, ivw_max_snp=1)
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


#' Identify putatively significant associations in the outlier scan
#' 
#' @param mr_threshold_method This is the argument to be passed to \code{p.adjust}. Default is "fdr". If no p-value adjustment is to be applied then specify "unadjusted"
#' @param mr_threshold Threshold to declare significance
#' @export
#' @return Same as outlier_scan but the candidate_exposure_mr and candidate_outcome_mr objects have an extra pval_adj and sig column each
outlier_sig <- function(outlierscan, mr_threshold_method = "fdr", mr_threshold = 0.05)
{
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
		outlierscan$exposure_candidate_mr$pval_adj <- p.adjust(outlierscan$exposure_candidate_mr$pval, mr_threshold_method)
		outlierscan$exposure_candidate_mr$sig <- outlierscan$exposure_candidate_mr$pval_adj < mr_threshold
	} else {
		outlierscan$candidate_outcome_mr$sig <- outlierscan$candidate_outcome_mr$pval < mr_threshold
		outlierscan$candidate_exposure_mr$sig <- outlierscan$candidate_exposure_mr$pval < mr_threshold
		outlierscan$exposure_candidate_mr$sig <- outlierscan$exposure_candidate_mr$pval < mr_threshold
	}

	message("Number of candidate - outcome associations: ", sum(outlierscan$candidate_outcome_mr$sig))
	message("Number of candidate - exposure associations: ", sum(outlierscan$candidate_exposure_mr$sig))
	message("Number of exposure - candidate associations: ", sum(outlierscan$exposure_candidate_mr$sig))
	return(outlierscan)
}

#' Plot results from outlier_scan in a network
#' 
#' Creates a simple network depicting the connections between outlier instruments, the original exposure and outcome traits, and the detected candidate associations
#' 
#' @param outlierscan Output from outlier_scan function
#' @export
#' @return Prints plot, and returns dataframe of the connections
outlier_network <- function(outlierscan)
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

	if(!"sig" %in% names(outlierscan$candidate_outcome_mr) | !"sig" %in% names(outlierscan$candidate_exposure_mr))
	{
		stop("Significant hits have not been specified. Please run outlier_sig.")
	}


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
		data.frame(
			from=ao$trait[ao$id == outlierscan$dat$id.exposure[1]],
			to=ao$trait[ao$id %in% subset(outlierscan$exposure_candidate_mr, sig)$id.exposure],
			what="Exposure to candidate traits"
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
}





simulate <- function(nid)
{
	# Simulate 3 outliers
	# scenario 1 - pleiotropy g -> u -> y
	# scenario 2 - confounder g -> u; u -> x; u -> y
	# (ignore for now) scenario 3 - pathway g -> u; x -> u; u -> y

	out <- list()

	library(simulateGP)

	gx <- make_geno(nid, 10, 0.5)
	gu1 <- cbind(gx[,1], make_geno(nid, 6, 0.5))
	gu2 <- cbind(gx[,2], make_geno(nid, 6, 0.5))
	gu3 <- cbind(make_geno(nid, 7, 0.5))

	u1 <- gu1 %*% c(1, runif(6)) + rnorm(nid, sd=6)
	u2 <- gu2 %*% c(1, runif(6)) + rnorm(nid, sd=6)
	x <- gx %*% c(0.5, runif(9)) + u1 * -3 + u2 * -3 + rnorm(nid, sd=16)
	# u3 <- gu %*% c(0.5, runif(6)) + x * -3 + rnorm(nid, sd=6)
	u3 <- gu %*% c(0.5, runif(6)) + rnorm(nid, sd=6)
	y <- x * 4 + u1 * -3 + u2 * -3 + u3 * -3 + rnorm(nid, sd=10)

	out$dat <- get_effs(x, y, gx, "X", "Y")
	# radial <- RadialMR::RadialMR(out$dat$beta.exposure, out$dat$beta.outcome, out$dat$se.exposure, out$dat$se.outcome, out$dat$SNP, "IVW", "YES", "NO", 0.05/nrow(out$dat), "NO")

	# radialmr(out$dat, 1:2)

	outliers <- as.character(1:2)
	# output$radialmr <- radial
	nout <- length(outliers)
	outliers <- subset(out$dat, SNP %in% outliers)$SNP
	out$outliers <- outliers
	out$id_list <- c("U1", "U2", "U3")

	temp1 <- gwas(u1, gx[,1:2])
	temp1$SNP <- 1:2
	temp1$outcome <- "U1"
	temp2 <- gwas(u2, gx[,1:2])
	temp2$SNP <- 1:2
	temp2$outcome <- "U2"
	temp <- rbind(temp1, temp2)
	out$search <- data.frame(SNP=temp$SNP, beta.outcome=temp$bhat, se.outcome=temp$se, pval.outcome=temp$pval, id.outcome=temp$outcome, outcome=temp$outcome, mr_keep=TRUE, effect_allele.outcome="A", other_allele.outcome="G", eaf.outcome=0.5)
	out$search$sig <- out$search$pval.outcome < 5e-8


	## get effects
	u1x <- get_effs(u1, x, gu1, "U1", "X")
	u1y <- get_effs(u1, y, gu1, "U1", "Y")
	u1x$SNP <- u1y$SNP <- c(1, paste0("u1_",1:6))
	u1x$effect_allele.exposure <- u1y$effect_allele.exposure <- "A"
	u1x$other_allele.exposure <- u1y$other_allele.exposure <- "G"
	u1x$effect_allele.outcome <- u1y$effect_allele.outcome <- "A"
	u1x$other_allele.outcome <- u1y$other_allele.outcome <- "G"
	u1x$eaf.exposure <- u1y$eaf.exposure <- 0.5
	u1x$eaf.outcome <- u1y$eaf.outcome <- 0.5

	u2x <- get_effs(u2, x, gu2, "U2", "X")
	u2y <- get_effs(u2, y, gu2, "U2", "Y")
	u2x$SNP <- u2y$SNP <- c(2, paste0("u2_",1:6))
	u2x$effect_allele.exposure <- u2y$effect_allele.exposure <- "A"
	u2x$other_allele.exposure <- u2y$other_allele.exposure <- "G"
	u2x$effect_allele.outcome <- u2y$effect_allele.outcome <- "A"
	u2x$other_allele.outcome <- u2y$other_allele.outcome <- "G"
	u2x$eaf.exposure <- u2y$eaf.exposure <- 0.5
	u2x$eaf.outcome <- u2y$eaf.outcome <- 0.5


	out$candidate_instruments <- subset(rbind(u1x, u2x), select=c(
		SNP, beta.exposure, se.exposure, id.exposure, exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, pval.exposure
	))
	out2 <- subset(out$search, sig)
	out$candidate_instruments <- group_by(out$candidate_instruments, id.exposure) %>%
		do({
			x <- .
			y <- subset(out2, id.outcome == x$id.exposure[1])
			x <- subset(x, !SNP %in% y$SNP)
			x
		})

	out$candidate_outcome <- subset(rbind(u1y, u2y), select=c(
		SNP, beta.outcome, se.outcome, id.outcome, outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, pval.outcome
	))

	out$candidate_outcome_dat <- suppressMessages(harmonise_data(out$candidate_instruments, out$candidate_outcome))
	out$candidate_outcome_mr <- suppressMessages(mr(out$candidate_outcome_dat, method_list="mr_ivw"))


	out$candidate_exposure <- subset(rbind(u1x, u2x), select=c(
		SNP, beta.outcome, se.outcome, id.outcome, outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, pval.outcome
	))

	out$candidate_exposure_dat <- suppressMessages(harmonise_data(out$candidate_instruments, out$candidate_exposure))
	out$candidate_exposure_mr <- suppressMessages(mr(out$candidate_exposure_dat, method_list="mr_ivw"))


	return(out)
}

