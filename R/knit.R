#' Knit report using working environment
#'
#' Warning: It is quite likely that this will be called within an Rmd file
#' implying a recursive call to \code{knit()}. This will generate "duplicate label"
#' errors for unlabelled chunks. To avoid this, all code chunks
#' in your Rmd file should be named.
#' Supposedly this error can also be avoided by setting the following option:
#' \code{options(knitr.duplicate.label = 'allow')}.
#' I tried this but it didn't seem to help.
#'
#' @param input_filename Rmd file.
#' @param output_filename Markdown or HTML output file.  An HTML file
#' is specified using the .htm, .html, .HTM or .HTML file extension.
#' When html is specified, a similarly named markdown file is also
#' generated.
#' All output files including cache and figures will appear in the
#' same folder as \code{output_filename}.
#' @param  ... Arguments to be passed to \code{\link[knitr:knit]{knitr::knit}}
#' @return NULL
#' @keywords internal
knit_report <- function(input_filename, output_filename, ...)
{
	requireNamespace("knitr", quietly = TRUE)
    requireNamespace("markdown", quietly = TRUE)
    output_filename <- normalizePath(output_filename)

    output_dir <- dirname(output_filename)
    if (!file.exists(output_dir))
        dir.create(output_dir)

    current_dir <- getwd()
    on.exit(setwd(current_dir))
    setwd(output_dir)

    name <- gsub("\\.[^.]+$", "", basename(output_filename))
    suffix <- gsub(".*\\.([^.]+)$", "\\1", output_filename)

    is.html <- tolower(suffix) %in% c("htm","html")
    is.pdf <- tolower(suffix) == "pdf"
    is.docx <- tolower(suffix) %in% c("doc", "docx", "word")
    is.md <- tolower(suffix) %in% c("md", "markdown")

    if (is.html)
        return(knitr::knit2html(input_filename, output=paste0(name, ".html"), envir=parent.frame(), ...))
    else if (is.md)
        return(knitr::knit(input_filename, output=paste0(name, ".md"), envir=parent.frame(), ...))
    else if (is.pdf)
    {        
        requireNamespace("rmarkdown", quietly = TRUE)
        return(rmarkdown::render(input_filename, rmarkdown::pdf_document(), intermediates_dir=getwd(), output_dir=getwd(), output_file=paste0(name, ".pdf"), clean = TRUE, envir=parent.frame(), ...))
    }
    else if (is.docx)
    {        
        requireNamespace("rmarkdown", quietly = TRUE)
        return(rmarkdown::render(input_filename, rmarkdown::word_document(), intermediates_dir=getwd(), output_dir=getwd(), output_file=paste0(name, ".docx"), clean = TRUE, envir=parent.frame(), ...))
    }
    else
        stop("Please choose a filename with pdf, html, docx or md suffix")
}


#' Generate MR report
#'
#' Using the output from the \code{mr} function this report will generate a report containing tables and graphs summarising the results.
#' A separate report is produced for each exposure - outcome pair that was analysed.
#'
#' @md
#' @param dat Output from \code{\link{harmonise_data}}
#' @param output_path Directory in which reports should be saved.
#' @param output_type Choose `"html"` or `"md"`. Default is `"html"`.
#' All output files including cache and figures will appear in the
#' folder specified in \code{output_path}.
#' @param author Author name.
#' @param study Study title.
#' @param path The filepath to the report template.
#' @param ... Extra options to be passed to \code{\link[knitr:knit]{knitr::knit}}.
#'
#' @export
#' @return NULL
mr_report <- function(dat, output_path = ".", output_type = "html", author = "Analyst", study = "Two Sample MR", path=system.file("reports", package="TwoSampleMR"), ...)
{
    message("Writing report as ", output_type, " file to ", output_path)

    message("Performing analysis")
    m <- list(
        mr = mr(dat),
        enrichment = enrichment(dat),
        directionality_test = directionality_test(dat),
        mr_heterogeneity = mr_heterogeneity(dat),
        mr_pleiotropy_test = mr_pleiotropy_test(dat),
        mr_singlesnp = mr_singlesnp(dat),
        mr_leaveoneout = mr_leaveoneout(dat)
    )

    message("Generating graphs")
    p <- list(
        mr_scatter_plot = mr_scatter_plot(m$mr, dat),
        mr_forest_plot = mr_forest_plot(m$mr_singlesnp),
        mr_funnel_plot = mr_funnel_plot(m$mr_singlesnp),
        mr_leaveoneout_plot = mr_leaveoneout_plot(m$mr_leaveoneout)
    )

    combinations <- plyr::ddply(dat, c("id.exposure", "id.outcome"), plyr::summarise, n=length(exposure), exposure=exposure[1], outcome=outcome[1])

    output_file <- array("", nrow(combinations))
    for(i in 1:nrow(combinations))
    {
        title <- paste(combinations$exposure[i], "against", combinations$outcome[i])
        tablist <- lapply(m[c("mr", "enrichment", "directionality_test", "mr_heterogeneity", "mr_pleiotropy_test")], function(x)
            {
                if(is.null(x))
                {
                    return(NULL)
                } else {
                   subset(x, id.exposure == combinations$id.exposure[i] & id.outcome == combinations$id.outcome[i], select=-c(id.exposure, id.outcome, exposure, outcome))
                }
            }
        )

        plotlist <- lapply(p, function(x) {
            d <- attributes(x)$split_labels
            index <- which(d$id.exposure == combinations$id.exposure[i] & d$id.outcome == combinations$id.outcome[i])
            if(length(index) < 1)
            {
                return(blank_plot("Insufficient number of SNPs"))
            } else {
                return(x[[index]])
            }
        })

        output_file[i] <- file.path(output_path, paste("TwoSampleMR", sanitise_string(title), output_type, sep="."))
        output_file[i] <- knit_report(file.path(path, "mr_report.Rmd"), output_file[i], ...)
    }
    return(output_file)
}

sanitise_string <- function(x)
{
    gsub(" ", "_", gsub("[^[:alnum:] ]", "", x))
}