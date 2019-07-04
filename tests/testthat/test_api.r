library(devtools)
load_all()
toggle_api("dev")
api_query('status')
api_query('gwasinfo/2')
api_query('gwasinfo/2,1001')
api_query('gwasinfo') %>% nrow
api_query('gwasinfo', access_token=NULL) %>% nrow
api_query('gwasinfo', query=list(id=c(2,1001)))
api_query('gwasinfo', query=list(id=c(2,1001,987)))
api_query('gwasinfo', query=list(id=c(2,1001,987)), access_token=NULL)
api_query('associations/2/rs234')
api_query('associations', query=list(id=c(2,1001,987), rsid=c('rs234', 'rs123')))

o <- api_query('tophits', query=list(id=2))
api_query('tophits', query=list(id=987))
api_query('tophits', query=list(id=987), access_token=NULL)

p <- api_query('ldmatrix', query=list(rsid=o$name))
q <- api_query('clump', query=list(rsid=o$name, pval=o$p))

gwasinfo() %>% dim
gwasinfo(access_token=NULL) %>% dim
gwasinfo(id=c(2,987)) %>% dim
gwasinfo(id=c(123123123)) %>% dim
gwasinfo(id=c(2,987), access_token=NULL) %>% dim


a <- extract_instruments(2)
b <- extract_outcome_data(a$SNP, 7)


a1 <- extract_instruments(2, force_server=TRUE)
dim(a1)
dim(a)

a1 <- extract_instruments(c(300,301), force_server=TRUE)
table(a1$exposure)

a <-  extract_instruments(c(300,301))
table(a$exposure)
a$SNP[!a$SNP %in% a1$SNP]
a1$SNP[!a1$SNP %in% a$SNP]
a1$SNP[duplicated(a1$SNP)]
a$SNP[duplicated(a$SNP)]
