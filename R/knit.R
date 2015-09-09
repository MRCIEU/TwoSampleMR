#' Knit report using working environment
#'
#' Warning: It is quite likely that this will be called within an RMD file
#' implying a recursive call to knit(). This will generate "duplicate label"
#' errors for unlabelled chunks. To avoid this, all code chunks
#' in \code{rmd.filename} should be named.
#' Supposedly this error can also be avoided by setting the following option:
#'      options(knitr.duplicate.label='allow')
#' I tried this but it didn't seem to help.
#'
#' @param input_filename Rmd file.
#' @param output_filename Markdown or HTML output file.  An HTML file
#' is specified using the .htm, .html, .HTM or .HTML file extension.
#' When html is specified, a similarly named markdown file is also
#' generated.
#' All output files including cache and figures will appear in the
#' same folder as \code{output_filename}.
#' @param  ... Arguments to be passed to \code{\link{knitr::knit}}
#' @return NULL
knit_report <- function(input_filename, output_filename, ...)
{
	require(knitr)
    require(markdown)
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

    if (is.html)
        md.filename <- paste(name, "md", sep=".")
    else
        md.filename <- basename(output_filename)
    
    knit(input_filename, output=md.filename, envir=parent.frame(), ...)

    if (is.html)
        markdownToHTML(md.filename, basename(output_filename))
}


#' Generate MR report
#'
#' Using the output from the \code{mr} function this report will generate a report containing tables and graphs summarising the results.
#' A separate report is produced for each exposure - outcome pair that was analysed
#'
#' @param mr_results Output from \code{mr}
#' @param dat Output from \code{harmonise_exposure_outcome}
#' @param output_path Directory in which reports should be saved
#' @param output_type Choose "html" or "md". Default is "html".
#' All output files including cache and figures will appear in the
#' folder specified in \code{output_path}.
#' @param author Author name
#' @param study Study title
#' @param ... Extra options to be passed to knitr
#'
#' @export
#' @return NULL
mr_report <- function(mr_results, dat, output_path = ".", output_type = "html", author = "Analyst", study = "Two Sample MR", path=system.file("reports", package="meffil"), ...)
{
    message("Writing report as html file to ", output_path)

    combinations <- unique(paste(dat$exposure, dat$outcome, sep="@@@@@@"))
    combinations <- as.data.frame(do.call(rbind, strsplit(combinations, split="@@@@@@")))
    names(combinations) <- c("exposure", "outcome")

    plots <- mr_scatter_plot(mr_results, dat)
    combinations <- expand.grid(exposure=unique(dat$exposure), outcome=unique(dat$outcome))

    maintab <- lapply(1:nrow(combinations), function(i) {
        mrres <- subset(mr_results$mr, Exposure==combinations$exposure[i] & Outcome==combinations$outcome[i], select=c(Test, b, se, pval))
        names(mrres) <- c("Test", "Effect", "SE", "p-value")
        return(mrres)
    })

    eggertab <- lapply(1:nrow(combinations), function(i) {

        egger <- subset(mr_results$extra, Exposure==combinations$exposure[i] & Outcome==combinations$outcome[i], select=c(b, se, pval))
        names(egger) <- c("Effect", "SE", "p-value")
        return(egger)
    })

    # return(list(maintab, eggertab, plots))

    for(i in 1:nrow(combinations))
    {
        title <- paste(combinations$exposure[i], "against", combinations$outcome[i])
        p <- plots[[i]]
        mt <- maintab[[i]]
        et <- eggertab[[i]]
        output_file <- paste("TwoSampleMR", gsub(" ", "_", title), output_type, sep=".")
        knit_report(file.path(path, "mr_report.Rmd"), output_file, ...)
    }


}