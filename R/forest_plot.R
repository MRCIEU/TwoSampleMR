### Functions to generate forest plots or scatterplots
### Version: 0.52
### Date: 20160315
### Authors: cl14022, ph14916


### This script does:
### 1. Read in a set of meta-analysis data, for moderate numbers of SNPs and exposures
### 2. Count the numbers of exposures (SNPs) and of outcomes
### 3. Detect the presence of summary data
### 4. Determine the appropriate plot according to the following table
###     Scatterplot: if the number of exposures and outcomes together exceeds 50 ( or go straight to summaries if it's 1:many)
###     Grouped forest plot: if there are multiple outcomes (regardless of the number of exposures)
###     Forest plot with summaries: if there are multiple exposures and a single outcome ( or multiple outcomes and one exposure)
### 5. Produces the appropriate plot in the R environment
### 6. If a filename is given, saves the appropriate plot as a png


### To do this, it uses several sub-functions
### 1. one_exp_one_out(...,summaryOnly=FALSE): produces a forest plot with multiple SNP results and multiple summaries, or optionally only the summaries
### 2. one_exp_many_out(): produces a forest plot  of summary results based on one exposure and several outcomes
### 3. many_exp_one_out(): produces a forest plot  of summary results based on several exposures and one outcome
### 4. many_exp_many_out(...,snpsOrsummary="snps"): [NOT IMPLEMENTED YET]--produces a scatterplot of effect estimates, either SNPs or summaries
### 5. output_forestplot(): produces output files and R plots




##### Begin one_exp_many_out()
one_exp_many_out <- function( inpt_df = NULL ,  eff_Col = "effect.outcome", exposure_Name="snps.test", outcome_Name="effect.estimate", instrument_Size = 'n.snps',forest_Title = 'Effect', outfile_Name = 'annot_FP.pdf', left_Col_Names=c("outcome", "snps.test", "chr.names"), left_Col_Titles = c("Outcome", "SNP", "Chr"), right_Col_Names = NULL, right_Col_Titles = NULL, log_ES = FALSE, exp_ES = FALSE, se_Col = "se.outcome", ub_Col =  NULL, lb_Col = NULL, summary_list = c("Inverse variance weighted")) {
    library(ggplot2)
    library(plyr)
    library(gtable)
    library(reshape)
    library(scales)

    # inpt_df = (character) filename of the data
    # eff_Col = (character) character giving the MR effect column
    # exposure_Name = (character) name of the SNP Name/Summary Name column
    # outcome_Name = (character) name of the column giving the outcome
    # instrument_Size = (character) name of the column giving the number of SNPs in each instrument
    # forest_Title = (character) title of the forest plot
    # outfile_Name = (character) name of the file to save the png of the forest plot to
    # left_Col_Names = (character) vector of the column names of the LHS annotations to be used from the original data file
    # left_Col_Titles = (character) vector of headings for the forest plot columns of the LHS
    # right_Col_Names = (character) vector of the column names of the RHS annotations to be used from the original data file
    # right_Col_Titles = (character) vector of headings for the forest plot columns of the RHS
    # log_ES = (logical) log-transform the effect size and confidence interval limits Y/N?
    # exp_ES = (logical) exponential-transform the effect size and confidence interval limits Y/N?
    # decrease = (logical) sort effect sizes in decreasing order Y/N?
    # se_Col = (character) column name of the standard error column in the original data
    # ub_Col =  (character) column name of the upper bound of CI column in the original data
    # lb_Col = (character) column name of the lower bound of CI column in the original data

    options(warn=-1)
    keep_rows <- NULL
    for (i in summary_list){
        keep_rows <- c(which(i == inpt_df[,"method"]), keep_rows)
    }
    number_Exp <- length(unique(inpt_df[, exposure_Name])) # Number of exposures
    number_Out <- length(unique(inpt_df[, outcome_Name])) # Number of outcomes
    total_Rows <- nrow(inpt_df)

    inpt_df <- inpt_df[keep_rows,]

    if(length(lb_Col)  == 0){
        inpt_df$lb <- inpt_df[,eff_Col] - qnorm(p = 0.975) * inpt_df[,se_Col]
        inpt_df$ub <- inpt_df[,eff_Col] + qnorm(p = 0.975) * inpt_df[,se_Col]
        lb_Col <- "lb"
        ub_Col <- "ub"
        }


    spacer <- function( eff_col, outcome, Data_Fm) {

        ## This function fills in a 'spacing vector', to determine placement of rows in a study plot, an attribute list of textual attributes
        ## which ensures that subgroups have bold face headers, a content list, which gives the content for headers, and a row list, of the rows to be filled in with this header data


        # eff_col = (character), column name of the column giving the effect size of interest
        # outcome = (character), column name giving phenotype
        # Data_Fm = (name) name of the data frame containing the meta-analytic cat
        # decrease = (logical) order the SNPs by decreasing effect size y/n?
        # summary_rows = (numeric), column listing rows with summaries in them

        N_study <- nrow(Data_Fm)  ## includes summary data

        # Count the number of categories of outcomes; in this case there is only 1
        N_sub <- 1
        ## one for each SNP, and then one for each summary

        # The layout of rows is one per SNP, then, buffering room added for each outcome and also for summaries
        N_space <- N_study + 2 * N_sub

        spacing_vec <- vector(mode = "numeric", length = N_space)
        content_list <- list()
        attr_list <- list()
        row_list <- list()
        idx <- 1
        for (snp in unique(Data_Fm[,outcome])) {
            # loop over snp ids to fill the 'spacing vector,' which determines placement of the rows in the forest plot
            snp_idxs <- which(Data_Fm[,outcome] == snp)
            n_snp <- length(snp_idxs)
            spacing_vec[[idx]] <- idx
            spacing_vec[[(idx + 1 + n_snp)]] <- (idx + 1 + n_snp)

            spacing_vec[((idx + 1):(idx + n_snp))] <- idx + 1:n_snp

            # In the same loop, fill a list with study names but also with headers to serve as the primary annotation of the forest plot
            content_list[[idx]] <- snp
            content_list[[(idx + 1 + n_snp)]] <- ""
            content_list[((idx + 1):(idx + n_snp))] <- Data_Fm[snp_idxs, outcome]

            # A list of attributes for the typeface so that headers are in bold face
            attr_list[[idx]] <- "bold"
            attr_list[[(idx + 1 + n_snp)]] <- "plain"
            attr_list[((idx + 1):(idx + n_snp))] <- "plain"

            # The most important list, contains mapping between spaces on the forest plot and rows in the meta-analytic data frame
            row_list[[idx]] <- NA
            row_list[[(idx + 1 + n_snp)]] <- NA
            row_list[((idx + 1):(idx + n_snp))] <- snp_idxs

            idx <- idx + 2 + n_snp
        }

        ## repeat for summary_rows
        # returns data frame giving spacings, primary annotation (content_list), typeface attributes for the primary annotation (attr_list), and  mapping between rows in forest plot and rows in meta-analytic data frame
        return( data.frame( spacing_vec = (1 + length(spacing_vec) - spacing_vec), content_list = unlist(content_list), row_list = unlist(row_list), attr_list = unlist(attr_list) ) )
    }


    space_Out <- function(data_Fm, space_Fm) {
        # a function to expand the initial data frame and add the spacing vector in a proper way
        # data_Fm = (name) data frame containing the meta-analytic data
        # space_Fm = (name) data frame containing the spacing vector and the mapping between spaces in the forest plot and rows in the meta-analytic data
        exp_data_Fm  <- data.frame(space_Fm)
        for (jj in colnames(data_Fm)) {
            exp_data_Fm[,jj] <- data_Fm[space_Fm$row_list,jj]
        }
        # Returns a data frame containing spacing, primary annotations, typeface attributes for the primary annotation, and correctly-mapped rows in the meta-analytic data frame
        return(exp_data_Fm)
    }

    ### 4. Generate a properly-spaced forest plot (without a summary effect) for all studies
    ggforest <- function(data_Fm, space_col, eff_col, lb_col, ub_col,se_col, title_text = '', log_ES = FALSE, exp_ES = FALSE,  instrument_size) {
        # data_Fm = (name) spaced data frame (must be from output from space_Out()) containing meta-analytic data and appropriate annotation, spacing cols etc
        # space_col = (character) column name giving the column listing spacing of effects in the forest plot
        # eff_col = (character) column name of effect-size column
        # lb_col = (character) column name of the confidence interval lower bound column
        # ub_col = (character) column name of the confidence interval upper bound column
        # title_text = (character) title text for main part of forest plot (will use column name if no text is supplied)
        # log_ES = (logical) convert effect sizes TO natural log scale y/n?
        # exp_ES = (logical) convert effect sizes FROM natural log scale y/n?
        # summary_rows = (integer) vector listing rows with summaries in them
        # instrument_size = (integer) vector listing the number of SNPs to have gone into each instrument

        if (nchar(title_text) <= 1) {
            title_text <- eff_col
        }
        ## reassignment of column names to avoid ggplot scope problems
        data_Fm$space_col <- data_Fm[,space_col]

        data_Fm$lb_col <- as.numeric(data_Fm[,lb_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,lb_col]))*(exp_ES)

        data_Fm$ub_col <- as.numeric(data_Fm[,ub_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,ub_col]))*(exp_ES)

        data_Fm$eff_col <- as.numeric(data_Fm[,eff_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,eff_col]))*(exp_ES)
        data_Fm$se_col <- as.numeric(data_Fm[,se_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,se_col]))*(exp_ES)

        data_Fm$nsnps <- data_Fm[,instrument_size]


        if(log_ES == TRUE){
            stopifnot(all(as.numeric(data_Fm[,lb_col])) >= 0) # throw error if there are negative limits

            data_Fm$lb_col <-  log(as.numeric(data_Fm[,lb_col]))

            data_Fm$ub_col <-  log(as.numeric(data_Fm[,ub_col]))

            data_Fm$eff_col <- log(as.numeric(data_Fm[,eff_col]))

            data_Fm$se_col <- log(as.numeric(data_Fm[,se_col]))
        }


        # ggplot code to generate the forest plot using geom_segments and geom_points and to make a relatively minimal theme
        raw_forest  <- ggplot(data = data_Fm, aes( y = space_col, yend = space_col, x = as.numeric(lb_col), xend = as.numeric(ub_col) )) + geom_segment()  + geom_point(aes( y = space_col,  x = as.numeric(eff_col), size = (se_col)^(-2)), shape = 15) + scale_size_continuous(range = c(1.25, 8.5))  + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),axis.line=element_line(), axis.line.y=element_blank(), title = element_text(size = 30), legend.position = 'none', axis.line.x = element_line(size = 1)) + expand_limits(y = c(data_Fm[,space_col] - 1, data_Fm[,space_col] + 2),x = c(min(as.numeric(lb_col),na.rm = TRUE), max(as.numeric(ub_col), na.rm = TRUE))) + labs(title = title_text) # returns ggplot2 object with the (un-annotated) forest plot
        if(exp_ES == TRUE){
            raw_forest <- raw_forest  + scale_x_continuous(trans = log2_trans()) + geom_vline(xintercept = 1)
        } else {
            raw_forest <- raw_forest + geom_vline(xintercept = 0)
        }
        return(raw_forest)
    }



    ### 5. Generate a properly-spaced set of annotations for the left-hand side (lhs) of the figure
    ### and
    ### 6. Generate a properly-spaced set of annotations for the right-hand side (rhs) of the figure
    ## 5.1. Function to generate a single column of annotations
    anot_col <- function(data_Fm, text_col, space_col, title_text = '', notitles) {
        # data_Fm = (name) space_Out data frame
        # text_col = (character) name of column containing annotations
        # title_text = (character) title of annotation column, defaults to the column name if nothing is entered

        data_Fm$space_col <- data_Fm[,space_col]
        data_Fm$text_col <- data_Fm[,text_col]
        data_Fm$text_col[!is.na(data_Fm$text_col)] <-  format(data_Fm$text_col[!is.na(data_Fm$text_col)],digits = 3, width = 16)

        # A hard rule to set the width of the annotation column, which sometimes truncates very wide columns (complex disease names, numbers with 16 digits, etc)
        text_widths <- c(-1, max(10,0.5 * max(sapply( as.character(data_Fm[,text_col]),nchar ))))

        # GGplot rendering of the annotation column
        lefttext  <- ggplot(data = data_Fm, aes( y = space_col, x = 0, label = text_col, fontface = attr_list )) + geom_text(hjust = "inward") + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_text(colour = "white"),axis.ticks.x = element_line(colour = "white"), axis.title = element_blank(), rect = element_blank(), panel.grid = element_blank(), title = element_text(size = 30) ) + expand_limits(x = text_widths,  y = c(data_Fm[,space_col] - 1, data_Fm[,space_col] + 2)) + labs(title = title_text, size = 30)  # returns two-item list with left_text, the GGplot annotations, and text_widths, the x-axis limits of the plot
        if(notitles == TRUE){
            lefttext <- lefttext + theme(title = element_blank())
        }
        return(list(left_text = lefttext, text_widths = text_widths))
    }

    ## 5.2. Function to aggregate all of the left-hand (all of the right-hand) notation columns into a single subplot, which works by applying the anot_col function to each item in the list of annotation columns
    anot_side <- function(data_Fm, space_col, col_names, title_list = '', notitles = FALSE) {
        # data_Fm = (name) space_Out data frame containing effect sizes, all annotations, and spacings
        # space_col = (character) name of the column containing the spacing vector
        # col_names = (character vector) vector of the column names for the columns to be the left-hand-side or right-hand-side annotations
        # title_list = (character vector) vector of the titles to be given to each of the annotation columns, defaults to the column names

        relative_widths <- vector(mode = "numeric", length = length(col_names))
        output <- vector(mode = "list")

        if (length(title_list) <= 1) {
            title_list <- col_names
        }

        for (i in 1:length(col_names)) {
            # loop to get the widths of each annotation column and to group the annotation objects together
            col <- anot_col( data_Fm = data_Fm, text_col = col_names[i], space_col = space_col, title_text = title_list[[i]], notitles = notitles)
            relative_widths[i] <- col$text_widths[2] - col$text_widths[1]
            output[[col_names[i]]] <- ggplotGrob(col$left_text)
        }
        # calculate the widths of the annotation columns, relative to each other, and then multiply it by a third so that it fits onto one side of the plot
        output$relative_widths <- relative_widths / (sum(relative_widths)) * 0.29
        # returns an output object containing: each annotation plot, listed under its own name, and the relative_widths of the columns
        return(output)
    }




    ### 7. Combine the three sub-figures into a single figure
    ### 7.1.Generate a list of spacings for annotation columns and for


    ## order and structure effect sizes and CIs for forest plot
    space1 <- spacer( outcome = outcome_Name, eff_col = eff_Col, Data_Fm = inpt_df)

    expand_data  <- space_Out(data_Fm = inpt_df, space_Fm = space1)


    ## Make the forest plot
    fo1  <- ggforest( data_Fm = expand_data, space_col = 'spacing_vec', eff_col = eff_Col, lb_col = lb_Col, ub_col = ub_Col, log_ES = log_ES, title_text = forest_Title, exp_ES = exp_ES,instrument_size = instrument_Size, se_col = se_Col )

    ## Construct left-hand-side annotations
    left <- anot_side( data_Fm = expand_data, space_col = 'spacing_vec', col_names = left_Col_Names, title_list = left_Col_Titles, notitles = FALSE )

    ## Construct right-hand-side annotations
    right <- NULL
    if (length(right_Col_Names) > 0){
      right <- anot_side(data_Fm = expand_data, space_col = 'spacing_vec', col_names = right_Col_Names, title_list = right_Col_Titles, notitles = FALSE )
      }

    options(warn=0)
    output_forestplot(forest = fo1, left =  left, right = right, outfile_Name = outfile_Name, nrows = nrow(expand_data), ncols = length(right_Col_Names) + length(left_Col_Names), plotType = "oEmO")
    return(list(forest = fo1, left = left, right = right))

}
##### END one_exp_many_out()
##### BEGIN many_exp_one_out()
many_exp_one_out <- function( inpt_df = NULL ,  eff_Col = "effect.outcome", exposure_Name="snps.test", outcome_Name="effect.estimate", instrument_Size = 'n.snps',forest_Title = 'Effect', outfile_Name = 'annot_FP.pdf', left_Col_Names=c("outcome", "snps.test", "chr.names"), left_Col_Titles = c("Outcome", "SNP", "Chr"), right_Col_Names = NULL, right_Col_Titles = NULL, log_ES = FALSE, exp_ES = FALSE,  se_Col = "se.outcome", ub_Col =  NULL, lb_Col = NULL, summary_list = c("Inverse variance weighted")) {
    library(ggplot2)
    library(plyr)
    library(gtable)
    library(reshape)
    library(scales)

    # inpt_df = (character) filename of the data
    # eff_Col = (character) character giving the MR effect column
    # exposure_Name = (character) name of the SNP Name/Summary Name column
    # outcome_Name = (character) name of the column giving the outcome
    # instrument_Size = (character) name of the column giving the number of SNPs in each instrument
    # forest_Title = (character) title of the forest plot
    # outfile_Name = (character) name of the file to save the png of the forest plot to
    # left_Col_Names = (character) vector of the column names of the LHS annotations to be used from the original data file
    # left_Col_Titles = (character) vector of headings for the forest plot columns of the LHS
    # right_Col_Names = (character) vector of the column names of the RHS annotations to be used from the original data file
    # right_Col_Titles = (character) vector of headings for the forest plot columns of the RHS
    # log_ES = (logical) log-transform the effect size and confidence interval limits Y/N?
    # exp_ES = (logical) exponential-transform the effect size and confidence interval limits Y/N?
    # decrease = (logical) sort effect sizes in decreasing order Y/N?
    # se_Col = (character) column name of the standard error column in the original data
    # ub_Col =  (character) column name of the upper bound of CI column in the original data
    # lb_Col = (character) column name of the lower bound of CI column in the original data

    options(warn=-1)

    inpt_df <- inpt_df[inpt_df[,instrument_Size] > 1,]

    keep_rows <- NULL
    for (i in summary_list){
        keep_rows <- c(which(i == inpt_df[,"method"]), keep_rows)
    }
    inpt_df <- inpt_df[keep_rows,]
    number_Exp <- length(unique(inpt_df[, exposure_Name])) # Number of exposures
    number_Out <- length(unique(inpt_df[, outcome_Name])) # Number of outcomes
    total_Rows <- nrow(inpt_df)


    if(length(lb_Col)  == 0){
        inpt_df$lb <- inpt_df[,eff_Col] - qnorm(p = 0.975) * inpt_df[,se_Col]
        inpt_df$ub <- inpt_df[,eff_Col] + qnorm(p = 0.975) * inpt_df[,se_Col]
        lb_Col <- "lb"
        ub_Col <- "ub"
    }


    spacer <- function(exposure, eff_col, outcome, Data_Fm) {

        ## This function fills in a 'spacing vector', to determine placement of rows in a study plot, an attribute list of textual attributes
        ## which ensures that subgroups have bold face headers, a content list, which gives the content for headers, and a row list, of the rows to be filled in with this header data

        # exposure = (character), column name of the column containing the name of the particular exposure (SNP for it)
        # eff_col = (character), column name of the column giving the effect size of interest
        # outcome = (character), column name giving phenotype
        # Data_Fm = (name) name of the data frame containing the meta-analytic cat
        # decrease = (logical) order the SNPs by decreasing effect size y/n?
        # summary_rows = (numeric), column listing rows with summaries in them

        N_study <- nrow(Data_Fm)  ## includes summary data

        # Count the number of categories of exposures
        N_sub <- 1

        # The layout of rows is one per SNP, then, buffering room added for each exposure and also for summaries
        N_space <- N_study + 2 * N_sub

        spacing_vec <- vector(mode = "numeric", length = N_space)
        content_list <- list()
        attr_list <- list()
        row_list <- list()
        idx <- 1
        for (snp in unique(Data_Fm[,exposure])) {
            # loop over snp ids to fill the 'spacing vector,' which determines placement of the rows in the forest plot
            snp_idxs <- which(Data_Fm[,exposure] == snp)
            n_snp <- length(snp_idxs)
            spacing_vec[[idx]] <- idx
            spacing_vec[[(idx + 1 + n_snp)]] <- (idx + 1 + n_snp)

            spacing_vec[((idx + 1):(idx + n_snp))] <- idx + 1:n_snp

            # In the same loop, fill a list with study names but also with headers to serve as the primary annotation of the forest plot
            content_list[[idx]] <- snp
            content_list[[(idx + 1 + n_snp)]] <- ""
            content_list[((idx + 1):(idx + n_snp))] <- Data_Fm[snp_idxs, outcome]

            # A list of attributes for the typeface so that headers are in bold face
            attr_list[[idx]] <- "bold"
            attr_list[[(idx + 1 + n_snp)]] <- "plain"
            attr_list[((idx + 1):(idx + n_snp))] <- "plain"

            # The most important list, contains mapping between spaces on the forest plot and rows in the meta-analytic data frame
            row_list[[idx]] <- NA
            row_list[[(idx + 1 + n_snp)]] <- NA
            row_list[((idx + 1):(idx + n_snp))] <- snp_idxs

            idx <- idx + 2 + n_snp
        }

        ## repeat for summary_rows
        # returns data frame giving spacings, primary annotation (content_list), typeface attributes for the primary annotation (attr_list), and  mapping between rows in forest plot and rows in meta-analytic data frame
        return( data.frame( spacing_vec = (1 + length(spacing_vec) - spacing_vec), content_list = unlist(content_list), row_list = unlist(row_list), attr_list = unlist(attr_list) ) )
    }


    space_Out <- function(data_Fm, space_Fm) {
        # a function to expand the initial data frame and add the spacing vector in a proper way
        # data_Fm = (name) data frame containing the meta-analytic data
        # space_Fm = (name) data frame containing the spacing vector and the mapping between spaces in the forest plot and rows in the meta-analytic data
        exp_data_Fm  <- data.frame(space_Fm)
        for (jj in colnames(data_Fm)) {
            exp_data_Fm[,jj] <- data_Fm[space_Fm$row_list,jj]
        }
        # Returns a data frame containing spacing, primary annotations, typeface attributes for the primary annotation, and correctly-mapped rows in the meta-analytic data frame
        return(exp_data_Fm)
    }

    ### 4. Generate a properly-spaced forest plot (without a summary effect) for all studies
    ggforest <- function(data_Fm, space_col, eff_col, lb_col, ub_col,se_col, title_text = '', log_ES = FALSE, exp_ES = FALSE,  instrument_size) {
        # data_Fm = (name) spaced data frame (must be from output from space_Out()) containing meta-analytic data and appropriate annotation, spacing cols etc
        # space_col = (character) column name giving the column listing spacing of effects in the forest plot
        # eff_col = (character) column name of effect-size column
        # lb_col = (character) column name of the confidence interval lower bound column
        # ub_col = (character) column name of the confidence interval upper bound column
        # title_text = (character) title text for main part of forest plot (will use column name if no text is supplied)
        # log_ES = (logical) convert effect sizes TO natural log scale y/n?
        # exp_ES = (logical) convert effect sizes FROM natural log scale y/n?
        # summary_rows = (integer) vector listing rows with summaries in them
        # instrument_size = (integer) vector listing the number of SNPs to have gone into each instrument

        if (nchar(title_text) <= 1) {
            title_text <- eff_col
        }
        ## reassignment of column names to avoid ggplot scope problems
        data_Fm$space_col <- data_Fm[,space_col]

        data_Fm$lb_col <- as.numeric(data_Fm[,lb_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,lb_col]))*(exp_ES)

        data_Fm$ub_col <- as.numeric(data_Fm[,ub_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,ub_col]))*(exp_ES)

        data_Fm$eff_col <- as.numeric(data_Fm[,eff_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,eff_col]))*(exp_ES)
        data_Fm$se_col <- as.numeric(data_Fm[,se_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,se_col]))*(exp_ES)

        data_Fm$nsnps <- data_Fm[,instrument_size]


        if(log_ES == TRUE){
            stopifnot(all(as.numeric(data_Fm[,lb_col])) >= 0) # throw error if there are negative limits

            data_Fm$lb_col <-  log(as.numeric(data_Fm[,lb_col]))

            data_Fm$ub_col <-  log(as.numeric(data_Fm[,ub_col]))

            data_Fm$eff_col <- log(as.numeric(data_Fm[,eff_col]))

            data_Fm$se_col <- log(as.numeric(data_Fm[,se_col]))
        }


        # ggplot code to generate the forest plot using geom_segments and geom_points and to make a relatively minimal theme
        raw_forest  <- ggplot(data = data_Fm, aes( y = space_col, yend = space_col, x = as.numeric(lb_col), xend = as.numeric(ub_col) )) + geom_segment()  + geom_point(aes( y = space_col,  x = as.numeric(eff_col), size = (se_col)^(-2)), shape = 15) + scale_size_continuous(range = c(1.25, 8.5))  + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),axis.line=element_line(), axis.line.y=element_blank(), title = element_text(size = 30), legend.position = 'none', axis.line.x = element_line(size = 1)) + expand_limits(y = c(data_Fm[,space_col] - 1, data_Fm[,space_col] + 2),x = c(min(as.numeric(lb_col),na.rm = TRUE), max(as.numeric(ub_col), na.rm = TRUE))) + labs(title = title_text) # returns ggplot2 object with the (un-annotated) forest plot
        if(exp_ES == TRUE){
            raw_forest <- raw_forest  + scale_x_continuous(trans = log2_trans()) + geom_vline(xintercept = 1)
        } else {
            raw_forest <- raw_forest + geom_vline(xintercept = 0)
        }
        return(raw_forest)
    }



    ### 5. Generate a properly-spaced set of annotations for the left-hand side (lhs) of the figure
    ### and
    ### 6. Generate a properly-spaced set of annotations for the right-hand side (rhs) of the figure
    ## 5.1. Function to generate a single column of annotations
    anot_col <- function(data_Fm, text_col, space_col, title_text = '', notitles) {
        # data_Fm = (name) space_Out data frame
        # text_col = (character) name of column containing annotations
        # title_text = (character) title of annotation column, defaults to the column name if nothing is entered

        data_Fm$space_col <- data_Fm[,space_col]
        data_Fm$text_col <- data_Fm[,text_col]
        data_Fm$text_col[!is.na(data_Fm$text_col)] <-  format(data_Fm$text_col[!is.na(data_Fm$text_col)],digits = 3, width = 16)

        # A hard rule to set the width of the annotation column, which sometimes truncates very wide columns (complex disease names, numbers with 16 digits, etc)
        text_widths <- c(-1, max(10,0.5 * max(sapply( as.character(data_Fm[,text_col]),nchar ))))

        # GGplot rendering of the annotation column
        lefttext  <- ggplot(data = data_Fm, aes( y = space_col, x = 0, label = text_col, fontface = attr_list )) + geom_text(hjust = "inward") + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_text(colour = "white"),axis.ticks.x = element_line(colour = "white"), axis.title = element_blank(), rect = element_blank(), panel.grid = element_blank(), title = element_text(size = 30) ) + expand_limits(x = text_widths,  y = c(data_Fm[,space_col] - 1, data_Fm[,space_col] + 2)) + labs(title = title_text, size = 30)  # returns two-item list with left_text, the GGplot annotations, and text_widths, the x-axis limits of the plot
        if(notitles == TRUE){
            lefttext <- lefttext + theme(title = element_blank())
        }
        return(list(left_text = lefttext, text_widths = text_widths))
    }

    ## 5.2. Function to aggregate all of the left-hand (all of the right-hand) notation columns into a single subplot, which works by applying the anot_col function to each item in the list of annotation columns
    anot_side <- function(data_Fm, space_col, col_names, title_list = '', notitles = FALSE) {
        # data_Fm = (name) space_Out data frame containing effect sizes, all annotations, and spacings
        # space_col = (character) name of the column containing the spacing vector
        # col_names = (character vector) vector of the column names for the columns to be the left-hand-side or right-hand-side annotations
        # title_list = (character vector) vector of the titles to be given to each of the annotation columns, defaults to the column names

        relative_widths <- vector(mode = "numeric", length = length(col_names))
        output <- vector(mode = "list")

        if (length(title_list) <= 1) {
            title_list <- col_names
        }

        for (i in 1:length(col_names)) {
            # loop to get the widths of each annotation column and to group the annotation objects together
            col <- anot_col( data_Fm = data_Fm, text_col = col_names[i], space_col = space_col, title_text = title_list[[i]], notitles = notitles)
            relative_widths[i] <- col$text_widths[2] - col$text_widths[1]
            output[[col_names[i]]] <- ggplotGrob(col$left_text)
        }
        # calculate the widths of the annotation columns, relative to each other, and then multiply it by a third so that it fits onto one side of the plot
        output$relative_widths <- relative_widths / (sum(relative_widths)) * 0.29
        # returns an output object containing: each annotation plot, listed under its own name, and the relative_widths of the columns
        return(output)
    }




    ### 7. Combine the three sub-figures into a single figure
    ### 7.1.Generate a list of spacings for annotation columns and for


    ## order and structure effect sizes and CIs for forest plot
    space1 <- spacer( exposure = exposure_Name, outcome = outcome_Name, eff_col = eff_Col, Data_Fm = inpt_df)

    expand_data  <- space_Out(data_Fm = inpt_df, space_Fm = space1)


    ## Make the forest plot
    fo1  <- ggforest( data_Fm = expand_data, space_col = 'spacing_vec', eff_col = eff_Col, lb_col = lb_Col, ub_col = ub_Col, log_ES = log_ES, title_text = forest_Title, exp_ES = exp_ES,instrument_size = instrument_Size, se_col = se_Col )

    ## Construct left-hand-side annotations
    left <- anot_side( data_Fm = expand_data, space_col = 'spacing_vec', col_names = left_Col_Names, title_list = left_Col_Titles, notitles = FALSE )

    ## Construct right-hand-side annotations
    right <- NULL
    if (length(right_Col_Names) > 0){
        right <- anot_side(data_Fm = expand_data, space_col = 'spacing_vec', col_names = right_Col_Names, title_list = right_Col_Titles, notitles = FALSE )
    }

    options(warn=0)
    output_forestplot(forest = fo1, left =  left, right = right, outfile_Name = outfile_Name, nrows = nrow(expand_data), ncols = length(right_Col_Names) + length(left_Col_Names), plotType = "mEoO")
    return(list(forest = fo1, left = left, right = right))

}
##### END many_exp_one_out()
##### BEGIN one_exp_one_out(...,summaryOnly=FALSE)
one_exp_one_out <- function( inpt_df = NULL ,  eff_Col = "effect.outcome", exposure_Name="snps.test", outcome_Name="effect.estimate", instrument_Size = 'n.snps',forest_Title = 'Effect', outfile_Name = 'annot_FP.pdf', left_Col_Names=c("outcome", "snps.test", "chr.names"), left_Col_Titles = c("Outcome", "SNP", "Chr"), right_Col_Names = NULL, right_Col_Titles = NULL, log_ES = FALSE, exp_ES = FALSE, decrease = TRUE, se_Col = "se.outcome", ub_Col =  NULL, lb_Col = NULL, summaryOnly=FALSE, summary_list = c("Inverse variance weighted")) {
    library(ggplot2)
    library(plyr)
    library(gtable)
    library(reshape)
    library(scales)

    # inpt_df = (character) filename of the data
    # eff_Col = (character) character giving the MR effect column
    # exposure_Name = (character) name of the SNP Name/Summary Name column
    # outcome_Name = (character) name of the column giving the outcome
    # instrument_Size = (character) name of the column giving the number of SNPs in each instrument
    # forest_Title = (character) title of the forest plot
    # outfile_Name = (character) name of the file to save the png of the forest plot to
    # left_Col_Names = (character) vector of the column names of the LHS annotations to be used from the original data file
    # left_Col_Titles = (character) vector of headings for the forest plot columns of the LHS
    # right_Col_Names = (character) vector of the column names of the RHS annotations to be used from the original data file
    # right_Col_Titles = (character) vector of headings for the forest plot columns of the RHS
    # log_ES = (logical) log-transform the effect size and confidence interval limits Y/N?
    # exp_ES = (logical) exponential-transform the effect size and confidence interval limits Y/N?
    # decrease = (logical) sort effect sizes in decreasing order Y/N?
    # se_Col = (character) column name of the standard error column in the original data
    # ub_Col =  (character) column name of the upper bound of CI column in the original data
    # lb_Col = (character) column name of the lower bound of CI column in the original data

    options(warn=-1)


    summary_Rows <- which(inpt_df[, instrument_Size] > 1)
    if (summaryOnly == FALSE){
        keep_rows <- NULL
        for (i in summary_list){
            keep_rows <- c(which(i == inpt_df[,"method"]), keep_rows)
        }

        inpt_df <- rbind(inpt_df[-summary_Rows,], inpt_df[intersect(summary_Rows, keep_rows),]) ## put summaries at the bottom of the plot
        summary_Rows <- (nrow(inpt_df) - length(intersect(summary_Rows, keep_rows)) + 1) : nrow(inpt_df)
    }
    number_Exp <- length(unique(inpt_df[, exposure_Name])) # Number of exposures
    number_Out <- length(unique(inpt_df[, outcome_Name])) # Number of outcomes
    total_Rows <- nrow(inpt_df)

    if(length(lb_Col)  == 0){
        inpt_df$lb <- inpt_df[,eff_Col] - qnorm(p = 0.975) * inpt_df[,se_Col]
        inpt_df$ub <- inpt_df[,eff_Col] + qnorm(p = 0.975) * inpt_df[,se_Col]
        lb_Col <- "lb"
        ub_Col <- "ub"
    }


    if(length(summary_Rows) > 0) {


    }

    if(summaryOnly == TRUE){
        summary_Rows <- NULL
    }

    spacer <- function(exposure, eff_col, outcome, Data_Fm, decrease = TRUE, summary_rows) {

        ## This function fills in a 'spacing vector', to determine placement of rows in a study plot, an attribute list of textual attributes
        ## which ensures that subgroups have bold face headers, a content list, which gives the content for headers, and a row list, of the rows to be filled in with this header data

        # exposure = (character), column name of the column containing the name of the particular exposure (SNP for it)
        # eff_col = (character), column name of the column giving the effect size of interest
        # outcome = (character), column name giving phenotype
        # Data_Fm = (name) name of the data frame containing the meta-analytic cat
        # decrease = (logical) order the SNPs by decreasing effect size y/n?
        # summary_rows = (numeric), column listing rows with summaries in them

        N_study <- nrow(Data_Fm)  ## includes summary data


        # The layout of rows is one per SNP, then, buffering room added for each exposure and also for summaries
        N_space <- N_study

        spacing_vec <- vector(mode = "numeric", length = N_space)
        content_list <- list()
        attr_list <- list()
        row_list <- list()
        idx <- 1
        for (snp in unique(Data_Fm[,exposure])) {
            # loop over snp ids to fill the 'spacing vector,' which determines placement of the rows in the forest plot
            snp_idxs <- which(Data_Fm[,exposure] == snp)
            n_snp <- length(snp_idxs)
            spacing_vec[[idx]] <- idx
            spacing_vec[[(idx + 1 + n_snp)]] <- (idx + 1 + n_snp)
            snp_sort <- c(order(as.numeric(Data_Fm[snp_idxs[setdiff(snp_idxs,summary_rows)], eff_col]), decreasing = decrease),(length(setdiff(snp_idxs,summary_rows)) + order(as.numeric(Data_Fm[snp_idxs[summary_rows], eff_col]), decreasing = decrease )) )
            spacing_vec[((idx + 1):(idx + n_snp))] <- idx + 1:n_snp

            # In the same loop, fill a list with study names but also with headers to serve as the primary annotation of the forest plot
            content_list[[idx]] <- snp
            content_list[[(idx + 1 + n_snp)]] <- ""
            content_list[((idx + 1):(idx + n_snp))] <- Data_Fm[snp_idxs[snp_sort], outcome]

            # A list of attributes for the typeface so that headers are in bold face
            attr_list[[idx]] <- "bold"
            attr_list[[(idx + 1 + n_snp)]] <- "plain"
            attr_list[((idx + 1):(idx + n_snp))] <- "plain"

            # The most important list, contains mapping between spaces on the forest plot and rows in the meta-analytic data frame
            row_list[[idx]] <- NA
            row_list[[(idx + 1 + n_snp)]] <- NA
            row_list[((idx + 1):(idx + n_snp))] <- snp_idxs[snp_sort]

            idx <- idx + 2 + n_snp
        }

        ## repeat for summary_rows
        # returns data frame giving spacings, primary annotation (content_list), typeface attributes for the primary annotation (attr_list), and  mapping between rows in forest plot and rows in meta-analytic data frame
        return( data.frame( spacing_vec = (1 + length(spacing_vec) - spacing_vec), content_list = unlist(content_list), row_list = unlist(row_list), attr_list = unlist(attr_list) ) )
    }


    space_Out <- function(data_Fm, space_Fm) {
        # a function to expand the initial data frame and add the spacing vector in a proper way
        # data_Fm = (name) data frame containing the meta-analytic data
        # space_Fm = (name) data frame containing the spacing vector and the mapping between spaces in the forest plot and rows in the meta-analytic data
        exp_data_Fm  <- data.frame(space_Fm)
        for (jj in colnames(data_Fm)) {
            exp_data_Fm[,jj] <- data_Fm[space_Fm$row_list,jj]
        }
        # Returns a data frame containing spacing, primary annotations, typeface attributes for the primary annotation, and correctly-mapped rows in the meta-analytic data frame
        return(exp_data_Fm)
    }

    ### 4. Generate a properly-spaced forest plot (without a summary effect) for all studies
    ggforest <- function(data_Fm, space_col, eff_col, lb_col, ub_col,se_col, title_text = '', log_ES = FALSE, exp_ES = FALSE, summary_rows, instrument_size) {
        # data_Fm = (name) spaced data frame (must be from output from space_Out()) containing meta-analytic data and appropriate annotation, spacing cols etc
        # space_col = (character) column name giving the column listing spacing of effects in the forest plot
        # eff_col = (character) column name of effect-size column
        # lb_col = (character) column name of the confidence interval lower bound column
        # ub_col = (character) column name of the confidence interval upper bound column
        # title_text = (character) title text for main part of forest plot (will use column name if no text is supplied)
        # log_ES = (logical) convert effect sizes TO natural log scale y/n?
        # exp_ES = (logical) convert effect sizes FROM natural log scale y/n?
        # summary_rows = (integer) vector listing rows with summaries in them
        # instrument_size = (integer) vector listing the number of SNPs to have gone into each instrument

        if (nchar(title_text) <= 1) {
            title_text <- eff_col
        }
        ## reassignment of column names to avoid ggplot scope problems
        data_Fm$space_col <- data_Fm[,space_col]

        data_Fm$lb_col <- as.numeric(data_Fm[,lb_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,lb_col]))*(exp_ES)

        data_Fm$ub_col <- as.numeric(data_Fm[,ub_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,ub_col]))*(exp_ES)

        data_Fm$eff_col <- as.numeric(data_Fm[,eff_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,eff_col]))*(exp_ES)
        data_Fm$se_col <- as.numeric(data_Fm[,se_col])*(1 - exp_ES) + exp(as.numeric(data_Fm[,se_col]))*(exp_ES)

        data_Fm$nsnps <- data_Fm[,instrument_size]

        data_Fm$is_summ <- 0
        if(length(summary_rows) > 0 ){ # variable showing whether something is a summary

            data_Fm$is_summ <- data_Fm$is_summ + data_Fm$nsnps > 1
        }

        data_Fm$is_summ <- factor(data_Fm$is_summ)
        if(log_ES == TRUE){
            stopifnot(all(as.numeric(data_Fm[,lb_col])) >= 0) # throw error if there are negative limits

            data_Fm$lb_col <-  log(as.numeric(data_Fm[,lb_col]))

            data_Fm$ub_col <-  log(as.numeric(data_Fm[,ub_col]))

            data_Fm$eff_col <- log(as.numeric(data_Fm[,eff_col]))

            data_Fm$se_col <- log(as.numeric(data_Fm[,se_col]))
        }


        # ggplot code to generate the forest plot using geom_segments and geom_points and to make a relatively minimal theme
        raw_forest  <- ggplot(data = data_Fm, aes( y = space_col, yend = space_col, x = as.numeric(lb_col), xend = as.numeric(ub_col) )) + geom_segment(aes(colour = is_summ))  + geom_point(aes( y = space_col,  x = as.numeric(eff_col), size = (se_col)^(-2), shape = is_summ, fill = is_summ)) + scale_shape_manual(values = c(15,18)) + scale_size_continuous(range = c(1.25, 8.5)) + scale_colour_manual(values = c('gray40', 'black'))  + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),axis.line=element_line(), axis.line.y=element_blank(), title = element_text(size = 30), legend.position = 'none', axis.line.x = element_line(size = 1)) + expand_limits(y = c(data_Fm[,space_col] - 1, data_Fm[,space_col] + 2),x = c(min(as.numeric(lb_col),na.rm = TRUE), max(as.numeric(ub_col), na.rm = TRUE))) + labs(title = title_text) # returns ggplot2 object with the (un-annotated) forest plot
        if(exp_ES == TRUE){
            raw_forest <- raw_forest  + scale_x_continuous(trans = log2_trans()) + geom_vline(xintercept = 1)
        } else {
            raw_forest <- raw_forest + geom_vline(xintercept = 0)
        }
        return(raw_forest)
    }



    ### 5. Generate a properly-spaced set of annotations for the left-hand side (lhs) of the figure
    ### and
    ### 6. Generate a properly-spaced set of annotations for the right-hand side (rhs) of the figure
    ## 5.1. Function to generate a single column of annotations
    anot_col <- function(data_Fm, text_col, space_col, title_text = '', notitles) {
        # data_Fm = (name) space_Out data frame
        # text_col = (character) name of column containing annotations
        # title_text = (character) title of annotation column, defaults to the column name if nothing is entered

        data_Fm$space_col <- data_Fm[,space_col]
        data_Fm$text_col <- data_Fm[,text_col]
        data_Fm$text_col[!is.na(data_Fm$text_col)] <-  format(data_Fm$text_col[!is.na(data_Fm$text_col)],digits = 3, width = 16)

        # A hard rule to set the width of the annotation column, which sometimes truncates very wide columns (complex disease names, numbers with 16 digits, etc)
        text_widths <- c(-1, max(10,0.5 * max(sapply( as.character(data_Fm[,text_col]),nchar ))))

        # GGplot rendering of the annotation column
        lefttext  <- ggplot(data = data_Fm, aes( y = space_col, x = 0, label = text_col, fontface = attr_list )) + geom_text(hjust = "inward") + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_text(colour = "white"),axis.ticks.x = element_line(colour = "white"), axis.title = element_blank(), rect = element_blank(), panel.grid = element_blank(), title = element_text(size = 30) ) + expand_limits(x = text_widths,  y = c(data_Fm[,space_col] - 1, data_Fm[,space_col] + 2)) + labs(title = title_text, size = 30)  # returns two-item list with left_text, the GGplot annotations, and text_widths, the x-axis limits of the plot
        if(notitles == TRUE){
            lefttext <- lefttext + theme(title = element_blank())
        }
        return(list(left_text = lefttext, text_widths = text_widths))
    }

    ## 5.2. Function to aggregate all of the left-hand (all of the right-hand) notation columns into a single subplot, which works by applying the anot_col function to each item in the list of annotation columns
    anot_side <- function(data_Fm, space_col, col_names, title_list = '', notitles = FALSE) {
        # data_Fm = (name) space_Out data frame containing effect sizes, all annotations, and spacings
        # space_col = (character) name of the column containing the spacing vector
        # col_names = (character vector) vector of the column names for the columns to be the left-hand-side or right-hand-side annotations
        # title_list = (character vector) vector of the titles to be given to each of the annotation columns, defaults to the column names

        relative_widths <- vector(mode = "numeric", length = length(col_names))
        output <- vector(mode = "list")

        if (length(title_list) <= 1) {
            title_list <- col_names
        }

        for (i in 1:length(col_names)) {
            # loop to get the widths of each annotation column and to group the annotation objects together
            col <- anot_col( data_Fm = data_Fm, text_col = col_names[i], space_col = space_col, title_text = title_list[[i]], notitles = notitles)
            relative_widths[i] <- col$text_widths[2] - col$text_widths[1]
            output[[col_names[i]]] <- ggplotGrob(col$left_text)
        }
        # calculate the widths of the annotation columns, relative to each other, and then multiply it by a third so that it fits onto one side of the plot
        output$relative_widths <- relative_widths / (sum(relative_widths)) * 0.29
        # returns an output object containing: each annotation plot, listed under its own name, and the relative_widths of the columns
        return(output)
    }




    ### 7. Combine the three sub-figures into a single figure
    ### 7.1.Generate a list of spacings for annotation columns and for


    ## order and structure effect sizes and CIs for forest plot
    space1 <- spacer( exposure = exposure_Name, eff_col = eff_Col, outcome = outcome_Name, Data_Fm = inpt_df, decrease = decrease, summary_rows = summary_Rows)

    expand_data  <- space_Out(data_Fm = inpt_df, space_Fm = space1)


    ## Make the forest plot
    fo1  <- ggforest( data_Fm = expand_data, space_col = 'spacing_vec', eff_col = eff_Col, lb_col = lb_Col, ub_col = ub_Col, log_ES = log_ES, title_text = forest_Title, exp_ES = exp_ES,summary_rows = summary_Rows,instrument_size = instrument_Size, se_col = se_Col )

    ## Construct left-hand-side annotations
    left <- anot_side( data_Fm = expand_data, space_col = 'spacing_vec', col_names = left_Col_Names, title_list = left_Col_Titles, notitles = FALSE )

    ## Construct right-hand-side annotations
    right <- NULL
    if (length(right_Col_Names) > 0){
        right <- anot_side(data_Fm = expand_data, space_col = 'spacing_vec', col_names = right_Col_Names, title_list = right_Col_Titles, notitles = FALSE )
    }

    options(warn=0)
    output_forestplot(forest = fo1, left =  left, right = right, outfile_Name = outfile_Name, nrows = nrow(expand_data), ncols = length(right_Col_Names) + length(left_Col_Names), plotType = ifelse(summaryOnly,"oEoOo","oEoO"))
    return(list(forest = fo1, left = left, right = right))

}
##### END one_exp_one_out(...,summaryOnly=FALSE)

##### BEGIN output_forestplot()
output_forestplot <- function(forest, left, right, outfile_Name, nrows, ncols, plotType){

    group_top_Plots <- function(forst_Pt, left_Hs, right_Hs = NULL) {
        # forst_Pt = (name) the name of the forest plot object
        # left_Hs = (name) the name of the aggregate annotation object for the left hand side of the plot
        # right_Hs  = (name) the name of the aggregate annotation object for the right hand side of the plot

        # Aggregate all of the plots into a single grid object for plotting
        grob_Bag <- vector(mode = 'list')

        left_RW <- left_Hs$relative_widths
        left_Grobs <- left_Hs
        left_Grobs$relative_widths <- NULL

        for (i in 1:length(left_Grobs)) {
            grob_Bag[paste('l',names(left_Grobs)[i],sep = '')] <- left_Grobs[i]
        }

        grob_Bag$m_forest <- ggplotGrob(forst_Pt)
        right_RW <- NULL
        if(!is.null(right_Hs)){
            right_RW <- right_Hs$relative_widths
            right_Grobs <- right_Hs
            right_Grobs$relative_widths <- NULL
            for (i in 1:length(right_Grobs)) {
                grob_Bag[paste('r',names(right_Grobs)[i], sep = '')] <- right_Grobs[i]
            }
        }


        width_vec <- c(left_RW,0.42,right_RW)
        width_vec <- width_vec / sum(width_vec)
        # convert the grid objects (now grouped) into a table of grid objects that can be plotted using grid.draw
        grp_FP <- gtable_matrix( name = "groupplot", grobs = matrix(grob_Bag, nrow = 1), widths = unit(width_vec, "npc"), heights = unit(1,"npc") )

        # return the grid object table, to be plotted
        return(grp_FP)
    }

    group <- group_top_Plots(forst_Pt = forest, left_Hs = left, right_Hs = right)
    ## draw and export the annotated forest plot
    grid.newpage()
    grid.draw(group)

    typeSpec <- data.frame(c('oEmO', 'mEoO', 'oEoO', 'oEoOo'),1:4)
    colnames(typeSpec) <- c('type', 'spec')
    widthSpec <- switch(typeSpec$spec[match(plotType,typeSpec$type)], 2500 + 180*(ncols%%30), 6000 + 680*(ncols%%30), 4000 + 350*(ncols%%30), 2500 + 180*(ncols%%30))
    heightSpec <- switch(typeSpec$spec[match(plotType,typeSpec$type)],3000+ 250*(nrows%%13), 3000+ 250*(nrows%%13), 6000+ 500*(nrows%%13), 2000+ 125*(nrows%%13))
    if(length(outfile_Name) > 0 ){
        png(file = outfile_Name,width = widthSpec, height = heightSpec, res = 300)
        grid.newpage()
        grid.draw(group)
        dev.off()
    }
    return(group)
}
##### END output_forestplot()
