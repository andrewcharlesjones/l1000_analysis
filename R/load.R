library(taigr)
library(magrittr)



#' Load L1000 compound data from Taiga (includes matched viability data)
#'
#' @param
#' @keywords
#' @export
#' @examples
#' load_l1000_cpd_data()
load_l1000_cpd_data <- function(normalize_sensitivity = TRUE) {
  l1k_data <- load.from.taiga(data.name='l1000-data-full-matrix-b0c3', data.version=2, data.file='l1000_drug_mat')
  gene_names <- colnames(l1k_data)[1:12328]
  
  sens_columns <- c("auc_prism", "auc_gdsc", "auc_ctrp")
  
  rows_to_keep <- which(rowSums(is.na(l1k_data[,sens_columns])) != 3)
  l1k_data <- l1k_data[rows_to_keep,]
  
  if (normalize_sensitivity) {
    
    l1k_data_collapsed <- l1k_data %>% 
      group_by(pert_iname, cell_id) %>% 
      summarise_all(mean) %>% 
      as.data.frame()
    
    sens_data <- l1k_data_collapsed[,sens_columns]
    sens_data_normalized <- l1000.analysis::quantile_normalize_viability_data(sens_data) %>% 
      as.data.frame() %>% 
      set_colnames(sens_columns)
    
    sens_data_normalized <- cbind(sens_data_normalized, l1k_data_collapsed[,c("pert_iname", "cell_id")])
    
    l1k_data_new_sens <- merge(l1k_data[,!(colnames(l1k_data) %in% sens_columns)], sens_data_normalized, 
                               by = c("pert_iname", "cell_id"))
    
    
    l1k_data <- l1k_data_new_sens
  }
  
  repurposing_metadata <- l1000.analysis::load_repurposing_cpd_annotations()
  repurposing_metadata <- repurposing_metadata[!duplicated(repurposing_metadata$name),]
  l1k_data <- merge(l1k_data, 
        repurposing_metadata[,c("name", "moa", "target")], 
        by.x = "pert_iname", 
        by.y = "name", 
        all.x = T)
  
  sig_metrics <- load.from.taiga(data.name='l1000-metadata-8489', data.version=1, data.file='GSE92742_Broad_LINCS_sig_metrics')
  l1k_data <- merge(l1k_data, sig_metrics[,c("sig_id", "tas")], by = "sig_id")
  
  
  tas_l2 <- sqrt(rowSums(l1k_data[,gene_names])^2)
  tas_l2 <- (tas_l2 - min(tas_l2)) / max(tas_l2 - min(tas_l2))
  l1k_data['tas_l2'] <- tas_l2

  
  list_to_return <- list(
    expression = l1k_data[,gene_names],
    sensitivity = 1 - l1k_data[,sens_columns],
    metadata = l1k_data[,c("sig_id", "pert_iname", "cell_id", "tas", "tas_l2", "pert_time", "pert_dose")]
  )
  
  list_to_return$sensitivity['auc_avg'] <- rowMeans(list_to_return$sensitivity[,sens_columns], na.rm = T)
  
  return(list_to_return)
}



#' Load L1000 shRNA data from Taiga (includes matched viability data)
#'
#' @param
#' @keywords
#' @export
#' @examples
#' load_l1000_shrna_data()
load_l1000_shrna_data <- function() {
  
  l1k_shrna_data <- load.from.taiga(data.name='l1000-shrna-data-full-matrix-798f', data.file='l1000_shrna_mat')
  
  gene_names <- colnames(l1k_shrna_data)[1:12328]
  
  rows_to_keep <- l1k_shrna_data[!is.na(l1k_shrna_data$shrna_viability),]
  
  list_to_return <- list(
    expression = l1k_shrna_data[,gene_names],
    sensitivity = -1 * l1k_shrna_data[,"shrna_viability", drop = F],
    metadata = l1k_shrna_data[,c("sig_id", "pert_iname", "cell_id", "pert_dose", "pert_time", "tas")]
  )
  
  return(list_to_return)
  
}


#' Load L1000 CRISPR data from Taiga (includes matched viability data)
#'
#' @param
#' @keywords
#' @export
#' @examples
#' load_l1000_crispr_data()
load_l1000_crispr_data <- function() {
  
  l1k_crispr_data <- load.from.taiga(data.name='l1000-crispr-data-full-matrix-8ed9', data.version=1, data.file='l1000_crispr_mat')
  
  gene_names <- colnames(l1k_crispr_data)[1:12328]
  
  rows_to_keep <- l1k_crispr_data[!is.na(l1k_crispr_data$shrna_viability),]
  
  list_to_return <- list(
    expression = l1k_crispr_data[,gene_names],
    sensitivity = 1 - l1k_crispr_data[,"crispr_viability", drop = F],
    metadata = l1k_crispr_data[,c("sig_id", "pert_iname", "cell_id")]
  )
  
  return(list_to_return)
  
}



#' Load L1000 gene metadata
#'
#' @param
#' @keywords
#' @export
#' @examples
#' load_gene_info()
load_gene_info <- function() {
  load.from.taiga(data.name='l1000-metadata-8489', data.version=1, data.file='gene_info')
}


#' Load compound annotations (MOA, target, etc.)
#'
#' @param
#' @keywords
#' @export
#' @examples
#' load_repurposing_cpd_annotations()
load_repurposing_cpd_annotations <- function() {
  repurposing.export.processed <- load.from.taiga(data.name='repurposing-compound-annotations-e8de', data.version=9, data.file='repurposing_export_processed')
  
  # load.from.taiga(data.name='repurposing-compound-annotations-0431', data.version=1, data.file='repurposing_export_10_4')
}





