
#' Check password for MR Base API
#'
#' @param password String
#' @export
#' @return Boolean whether password is accepted or not.
check_password <- function(password)
{
	url <- paste0("http://scmv-webapps.epi.bris.ac.uk:5000/check_password?password=", password)
	a <- fromJSON(url)
	if(a == 0)
	{
		message("Password accepted")
		return(TRUE)
	} else {
		message("Wrong password. You will not have access to some restricted studies.")
		return(FALSE)
	}
}


#' Set password as a global option
#'
#' You can set the password as a global option at the beginning of a session to avoid having to write it repeatedly in analysis
#'
#' @param password String
#' @export
#' @return Nothing, just prints message about whether password is correct or not, and sets global option
#' @alias
#' @examples \dontrun{
#'
#'}
set_password <- function(password)
{
	if(check_password(password)) options(password=password)
}


#' Get password from either global options or argument
#'
#' @param password String
#' @return Password string, or empty string if nothing supplied
get_password <- function(password)
{
	if(is.null(password))
	{
		if(is.null(options()$password))
		{
			return("")
		} else {
			check_password(options()$password)
			return(options()$password)
		}
	} else {
		check_password(password)
		return(password)
	}
}
