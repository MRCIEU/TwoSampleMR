# mysql -u epxjz -h epi-franklin.epi.bris.ac.uk -p
# wv-92n_YjB

# mysql -u mruser -h epi-franklin.epi.bris.ac.uk -p
# TMG_F1WnTL

# mysql -u gh13047 -h epi-franklin.epi.bris.ac.uk -p
# ri.K-2Gbvd

mydb <- dbConnect(MySQL(), user='epxjz', password='wv-92n_YjB', dbname='mrbase', host='epi-franklin.epi.bris.ac.uk')
dbListTables(mydb)
dbListFields(mydb, 'assoc')


rs <- dbSendQuery(mydb, "select * from study")
d <- dbFetch(rs, n=-1)


load("~/repo/mr_base/inst/data/tel.RData")
snps <- as.character(unique(tel$SNP))
te <- paste("+", paste(snps, collapse=" +"), sep="")

query <- paste("SELECT * FROM mrbase_assoc_info WHERE MATCH (snp_id) AGAINST ('", te, "' IN BOOLEAN MODE);", sep="")

query <- "SELECT * FROM assoc WHERE name='rs10936599'"
rs <- dbSendQuery(mydb, query)
d <- fetch(rs, n=10)
dim(d)


ssh -L 3306:localhost:3306 gh13047@epi-franklin.epi.bris.ac.uk
mysql -u gh13047 -h 127.0.0.1 -P 3306 -p
ri.K-2Gbvd

mysql -u mruser -h 127.0.0.1 -P 3306 -p
TMG_F1WnTL

use mrbase;

describe assoc;
describe snps;
describe study;

SELECT COUNT(*) FROM study;
SELECT COUNT(*) FROM snps;
# SELECT COUNT(*) FROM assoc;
# 1.7 billion rows

SELECT * FROM study limit 10;
SELECT * FROM snps WHERE name='rs13078807';
SELECT * FROM assoc WHERE snp=207707;

SELECT * FROM assoc limit 10;
SELECT * FROM assoc WHERE snp=2223704;

SELECT a.*, b.*, c.*
FROM assoc a, snps b, study c
WHERE a.snp=b.id AND a.study=c.id
AND (b.name='rs10900000' OR b.name='rs10000010' OR b.name='rs10000092')
AND (c.filename='cardiogramplusc4d_180814_update_data.txt.uniform.af.txt' OR c.filename='All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz.uniform.af.txt')
ORDER BY filename;


SELECT a.*, b.*, c.*
FROM assoc a, snps b, study c
WHERE a.snp=b.id AND a.study=c.id
AND b.name IN ('rs10900000', 'rs10000010', 'rs10000092')
AND c.filename IN ('cardiogramplusc4d_180814_update_data.txt.uniform.af.txt', 'All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz.uniform.af.txt', 'MAGIC_INSULIN_SECRETION_DI_for_release_HMrel27.txt.uniform.af.txt')
ORDER BY filename;
