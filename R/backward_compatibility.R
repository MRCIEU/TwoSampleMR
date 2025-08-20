ids_new_to_old <- function(id) {
  gsub("IEU-a:", "", id)
}

ids_new_to_old2 <- function(id) {
  id <- gsub("IEU-a-", "", id)
  id <- gsub("-a-", "-a:", id)
  id <- gsub("-b-", "-b:", id)
  id <- gsub("-c-", "-c:", id)
  id <- gsub("-d-", "-d:", id)
  return(id)
}


ids_old_to_new <- function(id) {
  ieu_index <- !grepl("[A-Z]+-[a-z]+:", id)
  id[ieu_index] <- paste0("IEU-a:", id[ieu_index])
  return(id)
}

ids_old_to_new2 <- function(id) {
  id <- gsub(":", "-", id)
  ieu_index <- !grepl("[A-Z]+-[a-z]+-", id)
  id[ieu_index] <- paste0("IEU-a-", id[ieu_index])
  id <- gsub(":", "-", id)
  return(id)
}
