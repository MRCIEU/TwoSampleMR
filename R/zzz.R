.onLoad <- function(libname, pkgname) {

	# packageStartupMessage(
	# 	"Welcome to TwoSampleMR.\n",
	# 	"[>] Check for updates: https://github.com/MRCIEU/TwoSampleMR\n",
	# 	"[>] Full documentation: https://mrcieu.github.io/TwoSampleMR"
	# )
	op <- options()
	op.googleAuthR <- list(
		googleAuthR.httr_oauth_cache = TRUE,
		googleAuthR.verbose = 3,
		googleAuthR.webapp.client_id = "906514199468-1jpkqgngur8emoqfg9j460s47fdo2euo.apps.googleusercontent.com",
		googleAuthR.webapp.client_secret = "I7Gqp83Ku4KJxL9zHWYxG_gD",
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
	options(mrbaseapi="http://api.mrbase.org/")
	
	invisible()

}