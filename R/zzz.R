.onLoad <- function(libname, pkgname) {

	a <- suppressWarnings(try(readLines("https://raw.githubusercontent.com/MRCIEU/TwoSampleMR/master/DESCRIPTION"), silent=TRUE))

	if(!class(a) == 'try-error')
	{
		latest <- gsub("Version: ", "", a[grep("Version", a)])
		current = utils::packageDescription('TwoSampleMR')

		test <- utils::compareVersion(latest, current$Version)
		if(test == 1)
		{
			packageStartupMessage("\nWarning:\nYou are running an old version of the TwoSampleMR package.\n",
				"This version:   ", current$Version, "\n",
				"Latest version: ", latest, "\n",
				"Please consider updating using devtools::install_github('MRCIEU/TwoSampleMR')")
		}
	}

	invisible()

}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "TwoSampleMR 0.5.1\n",
    "[>] IMPORTANT: Some datasets have been updated, and some are disabled while we update them.\n",
    "    Apologies for this inconvenience, they will be back up as soon as possible.\n",
    "    See news(package='TwoSampleMR') and https://gwas.mrcieu.ac.uk for more information\n",
    "[>] To temporarily revert to the previous database, see here:\n",
    "    https://mrcieu.github.io/TwoSampleMR/articles/gwas2020.html\n"
  )
}
