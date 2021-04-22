### Functions to generate grouped forest plots
### Version: 0.4
### Date: 20150914
### Authors: cl14022


### This script does:
### 1. Read in a set of meta-analysis data that  arise from mr_singlesnp
### 2. Order sub-analyses alphabetically
### 2.1. Order the studies within each sub-analysis by effect size
### 3. Generate a list of spacings so that sub-analyses can have headers and be spaced apart for easy reading
### 4. Generate a properly-spaced forest plot (without a summary effect) for all studies
### 5. Generate a properly-spaced set of annotations for the left-hand side (lhs) of the figure
### 6. Generate a properly-spaced set of annotations for the right-hand side (rhs) of the figure
### 7. Combine the three sub-figures into a single figure
### 7.1.Generate a list of spacings for annotation columns and for the forest-plot column
### 8. Export results as a pdf/windows meta file and to screen


# library(ggplot2)
# library(plyr)
# library(gtable)



## 8.0.1. Unifiying all into a single function that calls each of the above methods


#' Grouped forest plot
#'
#' @param name (character) name of the delimited file containing all of the results on the first sheet (needs to have headers), or of the r object.
#' @param eff_Col (character) name of the column in the delimited file that contains the effect sizes.
#' @param exposure_Name (character) name of the column in the delimited file containing the *types* of studies.
#' @param outcome_Name (character) name of the column in the delimited file containing the names of each study.
#' @param outfile_Name (character) name to be used for output file (*.pdf) or (*.wmf).
#' @param forest_Title (character) the title to be used for a forest plot.
#' @param left_Col_Names (character vector) vector containing the names of the left-hand-side annotation columns in the delimited file.
#' @param left_Col_Titles (character vector) vector containing the titles for each left-hand-side annotation column.
#' @param right_Col_Names (character vector) vector containing the names of the right-hand-side annotation columns in the delimited file.
#' @param right_Col_Titles (character vector) vector containing the titles for each right-hand-side annotation column.
#' @param debug (logical) show warnings `TRUE`/`FALSE`?
#' @param log_ES (logical) perform natural log transform of effect sizes and confidence bounds `TRUE`/`FALSE`?
#' @param decrease (logical) sort the studies by decreasing effect sizes `TRUE`/`FALSE`?
#' @param se_Col (character) name of the column giving the standard error of the effect sizes.
#' @param returnRobj (logical) return the graph as an internal R object `TRUE`/`FALSE`?
#' 
#' @keywords internal
#' @return grid object giving the forest plot (or plot as pdf)
#' @importFrom grDevices dev.off pdf
#' @importFrom stats qnorm
mr_forest_plot_grouped <-
    function(name, eff_Col = "b", exposure_Name="exposure", outcome_Name="outcome", forest_Title = '', outfile_Name = 'annot_FP.pdf', left_Col_Names=c("Exposure", "Outcome"), left_Col_Titles = NULL, right_Col_Names = c("p", "Outcome.n.case", "Outcome.n.control", "Outcome.sample.size"), right_Col_Titles =
                 NULL, debug = FALSE,  log_ES = FALSE, decrease = TRUE,  returnRobj = TRUE, se_Col = "se") {
        requireNamespace("gtable", quietly = TRUE)
        # name = (character) name of the r object from mr_singlesnp()
        # inRobj = (logical) is the data to be used an internal R object y/n?
        # eff_Col = (character) name of the column in the delimited file that contains the effect sizes

        # exposure_Name = (character) name of the column in the delimited file containing the exposure to use
        # outcome_Name = (character) name of the column in the delimited file containing outcome to be analyzed
        # outfile_Name = (character) name to be used for output file (*.pdf) or (*.wmf)
        # forest_Title = (character) the title to be used for a forest plot
        # left_Col_Names = (character vector) vector containing the names of the left-hand-side annotation columns in the delimited file
        # left_Col_Titles = (character vector) vector containing the titles for each left-hand-side annotation column
        # right_Col_Names = (character vector) vector containing the names of the right-hand-side annotation columns in the delimited file
        # right_Col_Titles = (character vector) vector containing the titles for each right-hand-side annotation column
        # debug = (logical) show warnings y/n?

        # log_ES = (logical) perform natural log transform of effect sizes and confidence bounds y/n?
        # decrease = (logical) sort the studies by decreasing effect sizes y/n?

        # returnRobj = (logical) return the graph as an internal R object y/n?
        name$lb <- name[,eff_Col] - qnorm(p = 0.95) * name[,se_Col]
        name$ub <- name[,eff_Col] + qnorm(p = 0.95) * name[,se_Col]
        data <- name
        spacer <- function(exposure, eff_col, outcome, Data_Fm, decrease = TRUE) {
                # exposure = (character), column name of the column containing the name of the particular exposure
                # eff_col = (character), column name of the column giving the effect size of interest
                # outcome = (character), column name giving study names/phenotype names
                # Data_Fm = (name) name of the data frame containing the meta-analytic cata
                # decrease = (logical) order the studies by decreasing effect size y/n?

                N_study <- nrow(Data_Fm)

                # Count the number of categories of studies (e.g. by type of illness )
                N_sub <- length(unique(Data_Fm[,exposure]))

                # The layout of rows is one per study, then, buffering room added for the subgroups
                N_space <- N_study + 2 * N_sub

                spacing_vec <- vector(mode = "numeric", length = N_space)
                content_list <- list()
                attr_list <- list()
                row_list <- list()
                idx <- 1
                for (subgp in unique(Data_Fm[,exposure])) {
                    # loop over substudy ids to fill the 'spacing vector,' which determines placement of the rows in the forest plot
                    sg_idxs <- which(Data_Fm[,exposure] == subgp)
                    n_std <- length(sg_idxs)
                    spacing_vec[[idx]] <- idx
                    spacing_vec[[(idx + 1 + n_std)]] <- (idx + 1 + n_std)
                    std_sort <- order(as.numeric(Data_Fm[sg_idxs, eff_col]), decreasing = decrease)
                    spacing_vec[((idx + 1):(idx + n_std))] <- idx + 1:n_std

                    # In the same loop, fill a list with study names but also with headers to serve as the primary annotation of the forest plot
                    content_list[[idx]] <- subgp
                    content_list[[(idx + 1 + n_std)]] <- ""
                    content_list[((idx + 1):(idx + n_std))] <- Data_Fm[sg_idxs[std_sort], outcome]

                    # A list of attributes for the typeface so that headers are in bold face
                    attr_list[[idx]] <- "bold"
                    attr_list[[(idx + 1 + n_std)]] <- "plain"
                    attr_list[((idx + 1):(idx + n_std))] <- "plain"

                    # The most important list, contains mapping between spaces on the forest plot and rows in the meta-analytic data frame
                    row_list[[idx]] <- NA
                    row_list[[(idx + 1 + n_std)]] <- NA
                    row_list[((idx + 1):(idx + n_std))] <- sg_idxs[std_sort]

                    idx <- idx + 2 + n_std
                }
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
        ggforest <- function(data_Fm, space_col, eff_col, lb_col = "lb", ub_col = "ub", title_text = '', log_ES = FALSE) {
                # data_Fm = (name) spaced data frame (must be from output from space_Out()) containing meta-analytic data and appropriate annotation, spacing cols etc
                # space_col = (character) column name giving the column listing spacing of effects in the forest plot
                # eff_col = (character) column name of effect-size column
                # lb_col = (character) column name of the confidence interval lower bound column
                # ub_col = (character) column name of the confidence interval upper bound column
                # title_text = (character) title text for main part of forest plot (will use column name if no text is supplied)
                # log_ES = (logical) convert effect sizes to natural log scale y/n?

                if (nchar(title_text) <= 1) {
                    title_text <- eff_col
                }
                ## reassignment of column names to avoid ggplot scope problems
                data_Fm$space_col <- data_Fm[,space_col]
                data_Fm$lb_col <- data_Fm[,lb_col]
                if (log_ES) {
                    data_Fm$lb_col <- log(as.numeric(data_Fm[,lb_col]))
                }
                data_Fm$ub_col <- data_Fm[, ub_col]
                if (log_ES) {
                    data_Fm$ub_col <- log(as.numeric(data_Fm[,ub_col]))
                }
                data_Fm$eff_col <- data_Fm[, eff_col]
                if (log_ES) {
                    data_Fm$eff_col <- log(as.numeric(data_Fm[,eff_col]))
                }
                # ggplot code to generate the forest plot using geom_segments and geom_points and to make a relatively minimal theme
                raw_forest  <- ggplot(data = data_Fm, aes( y = space_col, yend = space_col, x = as.numeric(lb_col), xend = as.numeric(ub_col) )) + geom_segment() + geom_point(aes( y = space_col,  x = as.numeric(eff_col), size = 4 )) + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), rect = element_blank(), title = element_text(size = 23), legend.position = 'none' ) + expand_limits(y = c(data_Fm[,space_col] - 1, data_Fm[,space_col] + 2)) + labs(title = title_text) # returns ggplot2 object with the (un-annotated) forest plot
                return(raw_forest)
            }


        ### 5. Generate a properly-spaced set of annotations for the left-hand side (lhs) of the figure
        ### and
        ### 6. Generate a properly-spaced set of annotations for the right-hand side (rhs) of the figure
        ## 5.1. Function to generate a single column of annotations
        anot_col <- function(data_Fm, text_col, space_col, title_text = '') {
                # data_Fm = (name) space_Out data frame
                # text_col = (character) name of column containing annotations
                # title_text = (character) title of annotation column, defaults to the column name if nothing is entered

                data_Fm$space_col <- data_Fm[,space_col]
                data_Fm$text_col <- data_Fm[,text_col]

                # A hard rule to set the width of the annotation column, which sometimes truncates very wide columns (complex disease names, numbers with 16 digits, etc)
                text_widths <- c(-1, max(10,0.5 * max(sapply( as.character(data_Fm[,text_col]),nchar ))))

                # GGplot rendering of the annotation column
                lefttext  <- ggplot(data = data_Fm, aes( y = space_col, x = 0, label = text_col, fontface = attr_list )) + geom_text(hjust = 0) + theme_bw() + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_text(colour = "white"),axis.ticks.x = element_line(colour = "white"), axis.title = element_blank(), rect = element_blank(), panel.grid = element_blank(), title = element_text(size = 23) ) + expand_limits(x = text_widths,  y = c(data_Fm[,space_col] - 1, data_Fm[,space_col] + 2)) + labs(title = title_text, size = 40)  # returns two-item list with left_text, the GGplot annotations, and text_widths, the x-axis limits of the plot
                return(list(left_text = lefttext, text_widths = text_widths))
            }

        ## 5.2. Function to aggregate all of the left-hand (all of the right-hand) notation columns into a single subplot, which works by applying the anot_col function to each item in the list of annotation columns
        anot_side <- function(data_Fm, space_col, col_names, title_list = '') {
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
                    col <- anot_col( data_Fm = data_Fm, text_col = col_names[i], space_col = space_col, title_text = title_list[[i]] )
                    relative_widths[i] <- col$text_widths[2] - col$text_widths[1]
                    output[[col_names[i]]] <- ggplotGrob(col$left_text)
                }
                # calculate the widths of the annotation columns, relative to each other, and then multiply it by a third so that it fits onto one side of the plot
                output$relative_widths <- relative_widths / (sum(relative_widths)) * 0.33
                # returns an output object containing: each annotation plot, listed under its own name, and the relative_widths of the columns
                return(output)
            }




        ### 7. Combine the three sub-figures into a single figure
        ### 7.1.Generate a list of spacings for annotation columns and for
        group_Plots <- function(forst_Pt, left_Hs, right_Hs) {
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

            right_RW <- right_Hs$relative_widths
            right_Grobs <- right_Hs
            right_Grobs$relative_widths <- NULL
            for (i in 1:length(right_Grobs)) {
                grob_Bag[paste('r',names(right_Grobs)[i], sep = '')] <- right_Grobs[i]
            }

            width_vec <- c(left_RW,0.34,right_RW)
            width_vec <- width_vec / sum(width_vec)
            # convert the grid objects (now grouped) into a table of grid objects that can be plotted using grid.draw
            grp_FP <- gtable::gtable_matrix( name = "groupplot", grobs = matrix(grob_Bag, nrow = 1), widths = unit(width_vec, "npc"), heights = unit(1,"npc") )

            # return the grid object table, to be plotted
            return(grp_FP)
        }


        ## act on debug flag
        if (debug == FALSE) {
            options(warn = -1)
        }

        ## read in data: beginning of the OVERALL FUNCTION



        ## order and structure effect sizes and CIs for forest plot
        space1 <- spacer( exposure = exposure_Name, eff_col = eff_Col, outcome = outcome_Name, Data_Fm = data, decrease = decrease)

        expand_data  <- space_Out(data_Fm = data, space_Fm = space1)

        ## Make the forest plot
        fo1  <- ggforest( data_Fm = expand_data, space_col = 'spacing_vec', eff_col = eff_Col, lb_col = "lb" ,ub_col = "ub", log_ES = log_ES, title_text = forest_Title )

        ## Construct left-hand-side annotations
        left <- anot_side( data_Fm = expand_data, space_col = 'spacing_vec', col_names = left_Col_Names, title_list = left_Col_Titles )

        ## Construct right-hand-side annotations
        right <- anot_side( data_Fm = expand_data, space_col = 'spacing_vec', col_names = right_Col_Names, title_list = right_Col_Titles )

        ## group all plots together
        group <- group_Plots(forst_Pt = fo1, left_Hs = left, right_Hs = right)
        ## draw and export the annotated forest plot
        grid.newpage()
        grid.draw(group)
        if (returnRobj == FALSE) {
            pdf(outfile_Name,width = 23.4, height = 16.5)
            grid.draw(group)
            dev.off()
        } else {
            return(group)
        }
        ## make sure not to leave warnings turned off
        options(warn = 0)
}
