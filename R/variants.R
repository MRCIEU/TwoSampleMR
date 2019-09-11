variants_gene <- function(gene, radius=0)
{
	l <- list()
	print(gene)
	for(i in 1:length(gene))
	{
		l[[gene[i]]] <- api_query(paste0('variants/gene/', gene[i], "?radius=", radius))
	}
	return(l)
}


variants_rsid <- function(rsid, radius)
{

}


variants_chrpos <- function(chrpos, radius)
{

}

