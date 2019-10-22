ids_new_to_old <- function(id)
{
	gsub("IEU-a:", "", id)
}

ids_new_to_old2 <- function(id)
{
	gsub("IEU-a-", "", id)
}


ids_old_to_new <- function(id)
{
	ieu_index <- ! grepl("[A-Z]+-[a-z]+:", id)
	id[ieu_index] <- paste0("IEU-a:", id[ieu_index])
	return(id)
}

ids_old_to_new2 <- function(id)
{
	id <- gsub(":", "-", id)
	ieu_index <- ! grepl("[A-Z]+-[a-z]+-", id)
	id[ieu_index] <- paste0("IEU-a-", id[ieu_index])
	return(id)
}

