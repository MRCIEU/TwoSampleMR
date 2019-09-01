.onLoad <- function(libname, pkgname) {

	packageStartupMessage(
		"Welcome to TwoSampleMR.\n",
		"[>] Full documentation: https://mrcieu.github.io/TwoSampleMR\n",
		"[>] Check the MR-Base database status using db_status()",
		"[>] Check news(package='TwoSampleMR') for bug fixes and updates.\n" 
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


	op <- options()
	op.googleAuthR <- list(
		googleAuthR.httr_oauth_cache = "mrbase.oauth",
		googleAuthR.verbose = 3,
		# googleAuthR.client_id = "906514199468-1jpkqgngur8emoqfg9j460s47fdo2euo.apps.googleusercontent.com",
		# googleAuthR.client_secret = "I7Gqp83Ku4KJxL9zHWYxG_gD",
		googleAuthR.webapp.client_id = "906514199468-1jpkqgngur8emoqfg9j460s47fdo2euo.apps.googleusercontent.com",
		googleAuthR.webapp.client_secret = "I7Gqp83Ku4KJxL9zHWYxG_gD",

		googleAuthR.client_id = "906514199468-m9thhcept50gu26ng494376iipt125d6.apps.googleusercontent.com",
		googleAuthR.client_secret = "zkihPnJnNRlHTinpzI0NUs4R",


		googleAuthR.webapp.port = 4018,
		googleAuthR.jsonlite.simplifyVector = TRUE,
		googleAuthR.scopes.selected = c("https://www.googleapis.com/auth/userinfo.profile",
										"https://www.googleapis.com/auth/userinfo.email"),
		googleAuthR.ok_content_types=c("application/json; charset=UTF-8", ("text/html; charset=UTF-8")),
		googleAuthR.securitycode = 
			paste0(sample(c(1:9, LETTERS, letters), 20, replace = T), collapse=''),
		googleAuthR.tryAttempts = 5
	)
	# toset <- !(names(op.googleAuthR) %in% names(op))
	# if(any(toset)) options(op.googleAuthR[toset])
	options(op.googleAuthR)

	# options(mrbaseapi="http://scmv-webapps.epi.bris.ac.uk:5000/")
	# options(mrbaseapi="http://apitest.mrbase.org/")
	toggle_api("dev")

	invisible()

}
