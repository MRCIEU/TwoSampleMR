library(TwoSampleMR)

a <- extract_instruments(2)
b <- clump_data(a)
message(nrow(b))
c <- ld_matrix(a$SNP)
message(nrow(c))
d <- extract_outcome_data(a$SNP, 7)
e <- available_outcomes()
message(nrow(e))






load("~/repo/mr-eve/results/03/exposure_dat.rdata")
snplist <- unique(exposure_dat$SNP)

t1 <- Sys.time()
options(mrbaseapi="http://ieu-db-interface.epi.bris.ac.uk:8080/")
out1 <- extract_outcome_data(snplist, 1001, proxies=FALSE)
Sys.time()-t1


t1 <- Sys.time()
toggle_api("release")
out2 <- extract_outcome_data(snplist, 1001, proxies=FALSE)
Sys.time()-t1


t1 <- Sys.time()
options(mrbaseapi="http://ieu-db-interface.epi.bris.ac.uk:8080/")
pout1 <- extract_outcome_data(snplist[1:5000], 1001, proxies=TRUE)
Sys.time()-t1


t1 <- Sys.time()
toggle_api("release")
pout2 <- extract_outcome_data(snplist[1:5000], 1001, proxies=TRUE)
Sys.time()-t1



library(dplyr)

t1 <- Sys.time()
options(mrbaseapi="http://ieu-db-interface.epi.bris.ac.uk:8080/")
pout1 <- lapply(split(snplist[1:5000], 1:5), function(x) extract_outcome_data(x, 1001, proxies=TRUE)) %>% bind_rows()
Sys.time()-t1
# 48sec

t1 <- Sys.time()
options(mrbaseapi="http://ieu-db-interface.epi.bris.ac.uk:8080/")
pout2 <- extract_outcome_data(snplist[1:5000], 1001, proxies=TRUE)
Sys.time()-t1
# 39sec

t1 <- Sys.time()
toggle_api("release")
pout2 <- lapply(split(snplist[1:5000], 1:5), function(x) extract_outcome_data(x, 1001, proxies=TRUE)) %>% bind_rows()
Sys.time()-t1
# 53sec

t1 <- Sys.time()
toggle_api("release")
pout2 <- extract_outcome_data(snplist[1:5000], 1001, proxies=TRUE)
Sys.time()-t1
# 40sec

####

t1 <- Sys.time()
options(mrbaseapi="http://ieu-db-interface.epi.bris.ac.uk:8080/")
pout1 <- lapply(split(snplist, 1:4), function(x) extract_outcome_data(x, 1001, proxies=TRUE)) %>% bind_rows()
Sys.time()-t1
# 3,7m

t1 <- Sys.time()
options(mrbaseapi="http://ieu-db-interface.epi.bris.ac.uk:8080/")
pout2 <- extract_outcome_data(snplist, 1001, proxies=TRUE, proxy_splitsize=50)
Sys.time()-t1
# 2.6m - 500


t1 <- Sys.time()
options(mrbaseapi="http://ieu-db-interface.epi.bris.ac.uk:8080/")
pout2 <- extract_outcome_data(snplist, 1001, proxies=TRUE, proxy_splitsize=100)
Sys.time()-t1
# 2.3m - 500



t1 <- Sys.time()
options(mrbaseapi="http://ieu-db-interface.epi.bris.ac.uk:8080/")
pout2 <- extract_outcome_data(snplist, 1001, proxies=TRUE)
Sys.time()-t1


t1 <- Sys.time()
toggle_api("release")
pout2 <- extract_outcome_data(snplist, 1001, proxies=TRUE)
Sys.time()-t1



t1 <- Sys.time()
toggle_api("release")
pout2 <- lapply(split(snplist[1:5000], 1:5), function(x) extract_outcome_data(x, 1001, proxies=TRUE)) %>% bind_rows()
Sys.time()-t1
# 53sec

t1 <- Sys.time()
toggle_api("release")
pout2 <- extract_outcome_data(snplist[1:5000], 1001, proxies=TRUE)
Sys.time()-t1
# 40sec
