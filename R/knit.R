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
#' Using the output from the \code{mr} function this report will generate a report containing tables and graphs summarising the results
#'
#' @param mr_results Output from \code{mr}
#' @param output_file Markdown or HTML output file.  An HTML file
#' is specified using the .htm, .html, .HTM or .HTML file extension.
#' When html is specified, a similarly named markdown file is also
#' generated.
#' All output files including cache and figures will appear in the
#' same folder as \code{output_filename}.
#' @param author Author name
#' @param study Study title
#' @param ... Extra options to be passed to knitr
#'
#' @export
#' @return NULL
mr_report <- function(mr_results, output_file = "mr_report.md", author = "Analyst", study = "Two Sample MR", ...)
{
    message("Writing report as html file to", output_file)
    path <- system.file("reports", package="meffil")
    knit_report(file.path(path, "mr_report.Rmd"), output_file, ...)

}