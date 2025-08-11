.onAttach <- function(libname, pkgname) {

	packageStartupMessage(paste("TwoSampleMR version", utils::packageVersion("TwoSampleMR"), "\n"))

	b <- suppressWarnings(try(jsonlite::read_json("https://raw.githubusercontent.com/MRCIEU/opengwas/main/messages-twosamplemr.json"), silent=TRUE))
	if (!inherits(b, 'try-error')) {
		if (length(b) > 0) {
			o <- lapply(b, function(x) {
				# packageStartupMessage(" Message date: ", x[["date"]])
				sapply(x[["message"]], function(j) packageStartupMessage(paste(" ", j)))
			})
		}
	}

	a <- suppressWarnings(try(readLines("https://raw.githubusercontent.com/MRCIEU/TwoSampleMR/master/DESCRIPTION"), silent=TRUE))

	if (!inherits(a, 'try-error')) {
		latest <- gsub("Version: ", "", a[grep("Version", a)])
		current = utils::packageDescription('TwoSampleMR')

		test <- utils::compareVersion(latest, current$Version)
		if (test == 1) {
			packageStartupMessage("\nWarning:\nYou are running an old version of the TwoSampleMR package.\n",
				"This version:   ", current$Version, "\n",
				"Latest version: ", latest, "\n",
				"Please consider updating using remotes::install_github('MRCIEU/TwoSampleMR')")
		}
	}
}
