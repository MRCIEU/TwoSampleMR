mysql -u epxjz -h epi-franklin.epi.bris.ac.uk -p
wv-92n_YjB


mydb <- dbConnect(MySQL(), user='mruser', password='TMG_F1WnTL', dbname='mrbase', host='epi-franklin.epi.bris.ac.uk')
dbListTables(mydb)
dbListFields(mydb, 'assoc')


rs <- dbSendQuery(mydb, "select * from snps")
d <- dbFetch(rs, n=10)


load("~/repo/mr_base/inst/data/tel.RData")
snps <- as.character(unique(tel$SNP))
te <- paste("+", paste(snps, collapse=" +"), sep="")

query <- paste("SELECT * FROM mrbase_assoc_info WHERE MATCH (snp_id) AGAINST ('", te, "' IN BOOLEAN MODE);", sep="")

query <- "SELECT * FROM assoc WHERE name='rs10936599'"
rs <- dbSendQuery(mydb, query)
d <- fetch(rs, n=10)
dim(d)



use mrbase;

describe assoc;
describe snps;
describe study;

SELECT COUNT(*) FROM study;
SELECT COUNT(*) FROM snps;
# SELECT COUNT(*) FROM assoc;
# 1.7 billion rows

SELECT * FROM study limit 10;
SELECT * FROM snps WHERE name='rs10900000';
SELECT * FROM assoc WHERE snp=207707;


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
AND c.filename IN ('cardiogramplusc4d_180814_update_data.txt.uniform.af.txt', 'All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz.uniform.af.txt')
ORDER BY filename;
