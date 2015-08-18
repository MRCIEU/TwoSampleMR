


mysql -u epxjz -h epi-franklin.epi.bris.ac.uk -p
password is: wv-92n_YjB





The database is called epxjz_2sampleMR.

use epxjz_2sampleMR;

mydb <- dbConnect(MySQL(), user='epxjz', password='wv-92n_YjB', dbname='epxjz_2sampleMR', host='gh13047@epi-franklin.epi.bris.ac.uk')


dbListTables(mydb)

This will return a list of the tables in our connection. 

dbListFields(mydb, 'mrbase_assoc_info')



rs <- dbSendQuery(mydb, "select * from mrbase_file_info")
d <- dbFetch(rs, n=-1)


load("inst/data/tel.RData")
snps <- as.character(unique(tel$SNP))
te <- paste("+", paste(snps, collapse=" +"), sep="")

query <- paste("SELECT * FROM mrbase_assoc_info WHERE MATCH (snp_id) AGAINST ('", te, "' IN BOOLEAN MODE);", sep="")

query <- "SELECT * FROM mrbase_assoc_info WHERE snp_id='rs10936599'"
rs <- dbSendQuery(mydb, query)
d <- fetch(rs, n=10)
dim(d)


SELECT COUNT(*) FROM mrbase_assoc_info;
SELECT COUNT(*) FROM mrbase_assoc_info where snp_id='rs10900000';

SELECT COUNT(*) FROM mrbase_GWAScat_assoc_info


SELECT * FROM mrbase_assoc_info WHERE snp_id='rs10900000' LIMIT 10;
SELECT * FROM mrbase_file_info WHERE file_name LIKE '%data%' LIMIT 10;

SELECT * FROM mrbase_assoc_info WHERE snp_id LIKE 'rs109%' LIMIT 10;

rsid_study