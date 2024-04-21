.onAttach <- function(libname, pkgname) {

        packageStartupMessage(
                paste("TwoSampleMR version", utils::packageVersion("TwoSampleMR"), "\n"),
                "[>] New: Improved API performance and stability\n",
				"[>] New: Authentication required for all OpenGWAS queries from 1st May 2024\n",
				"[>] For guidance see https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication\n"
        )

	a <- suppressWarnings(try(readLines("https://raw.githubusercontent.com/MRCIEU/TwoSampleMR/master/DESCRIPTION"), silent=TRUE))

	if(!inherits(a, 'try-error'))
	{
		latest <- gsub("Version: ", "", a[grep("Version", a)])
		current = utils::packageDescription('TwoSampleMR')

		test <- utils::compareVersion(latest, current$Version)
		if(test == 1)
		{
			packageStartupMessage("\nWarning:\nYou are running an old version of the TwoSampleMR package.\n",
				"This version:   ", current$Version, "\n",
				"Latest version: ", latest, "\n",
				"Please consider updating using remotes::install_github('MRCIEU/TwoSampleMR')")
		}
	}
}
