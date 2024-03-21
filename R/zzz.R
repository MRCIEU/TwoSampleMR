.onAttach <- function(libname, pkgname) {

        packageStartupMessage(
                paste("TwoSampleMR version", utils::packageVersion("TwoSampleMR"), "\n"),
                "[>] New: Option to use non-European LD reference panels for clumping etc\n",
                "[>] Some studies temporarily quarantined to verify effect allele\n",
                "[>] See news(package = 'TwoSampleMR') and https://gwas.mrcieu.ac.uk for further details\n"
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
