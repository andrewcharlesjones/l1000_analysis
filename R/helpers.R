



#' Load gene set data
#'
#' @param gsc
#' @param gene_set_name
#' @keywords
#' @export
#' @examples
#' load_gsc_data()
load_gsc_data <- function() {
  
  gsc_data <- readr::read_rds(download.raw.from.taiga(data.name='msigdb-gene-set-collections-8453', data.version=1))
  return(gsc_data)
  
}



#' Fetch the gene names belonging to a certain gene set
#'
#' @param gsc
#' @param gene_set_name
#' @keywords
#' @export
#' @examples
#' get_gene_set_genes()
get_gene_set_genes <- function(gsc, collection_name, gene_set_name) {
  
  gene_set_data <- gsc[[collection_name]]
  
  ii <- 1
  for (curr_gs in gene_set_data@.Data) {
    if (curr_gs@setName == gene_set_name) {
      return(curr_gs@geneIds)
    }
    
    ii <- ii + 1
  }
  
}





