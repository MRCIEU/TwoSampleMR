#' Obtain variants around a gene
#'
#' Provide a gene identified, either Ensembl or Entrez, e.g. ENSG00000123374 or 1017
#'
#' @param gene Vector of genes, either Ensembl or Entrez, e.g. ENSG00000123374 or 1017
#' @param radius=0 radius around the gene region to include
#'
#' @export
#' @return data frame 
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


#' Obtain information about rsid
#'
#'
#' @param rsid Vector of rsids
#'
#' @export
#' @return data frame
variants_rsid <- function(rsid)
{
	api_query("variants/rsid", list(rsid = rsid))
}


#' Obtain information about chr pos and surrounding region
#'
#' For a list of chromosome and positions, finds all variants within a given radius
#'
#' @param chrpos list of <chr>:<pos> in build 37, e.g. c("3:46414943", "3:122991235")
#' @param radius Radius around each chrpos, default = 0
#'
#' @export
#' @return Data frame
variants_chrpos <- function(chrpos, radius=0)
{
	api_query("variants/chrpos", list(chrpos = chrpos, radius=radius)) %>% bind_rows()
}


