.onLoad <- function(libname, pkgname) {

        packageStartupMessage(
                "New version 0.5.0 of TwoSampleMR\n",
                "[>] Full documentation: https://mrcieu.github.io/TwoSampleMR\n",
                "[>] List of changes and information about backgwards compatibility here:\n",
                "    https://mrcieu.github.io/TwoSampleMR/articles/gwas2020.html\n",
                "[>] Browse available data at https://gwas.mrcieu.ac.uk\n"
        )

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
