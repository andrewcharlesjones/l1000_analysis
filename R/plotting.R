library(ggplot2)
library(ggrepel)



#' Plot viability signature output
#'
#' @param via_sig_output
#' @keywords
#' @export
#' @examples
#' make_via_sig_volcano()
make_via_sig_volcano <- function(via_sig_output, cpd_name) {
  
  # fdr_is_signif <- via_sig_output$viability_related$adj.P.Val < 0.1
  # via_sig_output$viability_related['fdr_is_signif'] <- fdr_is_signif
  
  p1 <- via_sig_output$viability_related %>%
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point() + 
    cdsr::theme_Publication() + 
    geom_text_repel(data = unique(rbind(
      via_sig_output$viability_related %>% arrange(P.Value) %>% head(5),
      via_sig_output$viability_related %>% arrange(-logFC) %>% head(7),
      via_sig_output$viability_related %>% arrange(logFC) %>% head(7))),
                    aes(label = Gene),
      size = 10) +
    labs(x = "Viability-related signature") + 
    ggtitle(cpd_name) +
    theme_scatterplot() + 
    geom_vline(xintercept = 0, linetype = "dashed")
  
  p2 <- via_sig_output$viability_unrelated %>%
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point() + 
    cdsr::theme_Publication() + 
    geom_text_repel(data = unique(rbind(
      via_sig_output$viability_unrelated %>% arrange(P.Value) %>% head(5),
      via_sig_output$viability_unrelated %>% arrange(-logFC) %>% head(8),
      via_sig_output$viability_unrelated %>% arrange(logFC) %>% head(8))),
      aes(label = Gene),
      size = 10) +
    labs(x = "Viability-unrelated signature") + 
    ggtitle(cpd_name) + 
    theme_scatterplot() + 
    geom_vline(xintercept = 0, linetype = "dashed")
  
  cowplot::plot_grid(p1, p2)

}



#' Plot viability-related and -unrelated signatures against each other
#'
#' @param via_sig_output
#' @keywords
#' @export
#' @examples
#' plot_via_sig_comparison()
plot_via_sig_comparison <- function(via_sig_output, gene_class_df = NULL, cpd_to_plot = NULL) {
  
  merged_sigs <- merge(via_sig_output$viability_related, via_sig_output$viability_unrelated,
                       by = "Gene",
                       suffixes = c("_slope", "_intercept"))
  
  if (!is.null(gene_class_df)) {
    merged_sigs <- merge(merged_sigs, gene_class_df, by = "Gene")
    
    merged_sigs %>% ggplot(aes(logFC_slope, logFC_intercept, color = gene_class)) + 
      geom_point() + 
      cdsr::theme_Publication() + 
      geom_text_repel(data = merged_sigs[merged_sigs$gene_class != "",],
        aes(label = Gene),
        size = 10) + 
      labs(x = "Viability-related signature", y = "Viability-unrelated signature") + 
      scale_color_manual(values=c("gray", "red", "blue")) + 
      l1000.analysis::theme_scatterplot() + 
      theme(legend.title=element_blank()) + 
      ggtitle(cpd_to_plot) + 
      theme(legend.text=element_text(size=30))
    
  } else {
    merged_sigs %>% ggplot(aes(logFC_slope, logFC_intercept)) + 
      geom_point() + 
      cdsr::theme_Publication() + 
      geom_text_repel(data = unique(rbind(
        merged_sigs %>% arrange(-abs(logFC_slope)) %>% head(10),
        merged_sigs %>% arrange(-abs(logFC_intercept)) %>% head(10))),
        aes(label = Gene),
        size = 10) + 
      labs(x = "Viability-related signature", y = "Viability-unrelated signature") + 
      l1000.analysis::theme_scatterplot()
  }
  
}


#' Change aesthetics of scatterplot to a standard style
#'
#' @param ggplot_obj
#' @keywords
#' @export
#' @examples
#' theme_scatterplot()
theme_scatterplot <- function() {
  
  theme(text = element_text(size=40))
  
}



#' Save plots in a standardized way
#'
#' @param ggplot_obj
#' @keywords
#' @export
#' @examples
#' save_publication_plot()
save_publication_plot <- function(ggplot_obj, filename) {
  
  ggsave(filename = filename, ggplot_obj, height = 10, width = 10)
  
}





#' Plot an expression-sensitivity relationship for one gene
#'
#' @param cpd
#' @param via_sig
#' @param data_list
#' @param gene
#' @param sens_column
#' @keywords
#' @export
#' @examples
#' plot_gene_auc_relationship()
plot_gene_auc_relationship <- function(cpd, via_sig, data_list, gene, sens_column = "auc_avg", is_plot_abline = T, is_plot_intercept_dot = T) {
  
  one_cpd_idx <- which(data_list$metadata$pert_iname == cpd)
  one_cpd_data <- cbind(data_list$expression[one_cpd_idx,],
                        data_list$sensitivity[one_cpd_idx,,drop = F],
                        data_list$metadata[one_cpd_idx,]) %>% as.data.frame()
  one_cpd_data_collapsed <- one_cpd_data %>% 
    group_by(cell_id) %>% 
    summarise_all(mean)
  
  p <- one_cpd_data_collapsed %>% ggplot(aes_string(sens_column, gene)) + 
    geom_point(size = 3) + 
    cdsr::theme_Publication() + 
    labs(x = "Sensitivity", y = paste(gene, "DE", sep = " ")) + 
    l1000.analysis::theme_scatterplot() + 
    geom_smooth(method = "lm") + 
    geom_vline(xintercept=0, linetype="dashed", size = 3)
  
  if (is_plot_abline) {
    p <- p + geom_abline(slope = via_sig$viability_related$logFC[via_sig$viability_related$Gene == gene], 
                    intercept = via_sig$viability_unrelated$logFC[via_sig$viability_unrelated$Gene == gene])
  }
  
  if (is_plot_intercept_dot) {
    p <- p + geom_point(aes(x=0, y=via_sig$viability_unrelated$logFC[via_sig$viability_unrelated$Gene == gene]), color="red", size = 8)
  }
  
  return(p)
  
}




#' Plot comparison of two compounds' viability signatures
#'
#' @param all_via_sig_df
#' @param cpd1
#' @param cpd2
#' @keywords
#' @export
#' @examples
#' plot_two_via_sigs()
plot_two_via_sigs <- function(all_via_sig_df, cpd1, cpd2) {
  
  comparison_df <- data.frame(
    cpd1 = all_via_sig_df[[cpd1]]$via_sig_coeff,
    cpd2 = all_via_sig_df[[cpd2]]$via_sig_coeff,
    gene = all_via_sig_df[[cpd1]]$Gene
  )

  comparison_df %>% ggplot(aes(cpd1, cpd2)) + 
    geom_point() + 
    cdsr::theme_Publication() + 
    geom_text_repel(data = unique(rbind(
      comparison_df %>% arrange(-abs(cpd1)) %>% head(10),
      comparison_df %>% arrange(-abs(cpd2)) %>% head(10))),
    aes(label = gene)) + 
    labs(x = cpd1, y = cpd2) + 
    l1000.analysis::theme_scatterplot()
}



#' Plot comparison of one cpd 
#'
#' @param all_via_sig_df
#' @param cpd1
#' @param cpd2
#' @keywords
#' @export
#' @examples
#' plot_cpd_and_shrna_sigs()
plot_cpd_and_shrna_sigs <- function(cpd_via_sig_df, shrna_via_sig_df, cpd_name, shrna_name, shared_pathway_genes = NULL) {
  
  comparison_df <- merge(cpd_via_sig_df[[cpd_name]][['viability_related']],
                         shrna_via_sig_df[[shrna_name]][['viability_related']],
                         by = "Gene",
                         suffixes = c("_cpd", "_shrna"))
  
  comparison_df['shared_pathway'] <- llply(comparison_df$Gene, function(x) {
    if (x %in% shared_pathway_genes) {
      return("Shared pathway")
    } else {
      return("")
    }
  }) %>% as.character()
  
  if (!is.null(shared_pathway_genes)) {
    
    comparison_df %>% ggplot(aes(logFC_cpd, logFC_shrna, color = shared_pathway)) + 
      geom_point() + 
      cdsr::theme_Publication() + 
      geom_text_repel(data =comparison_df[comparison_df$Gene %in% shared_pathway_genes,],
        aes(label = Gene),
        size = 10) + 
      labs(x = cpd_name, y = paste(shrna_name, "KD", sep = " ")) + 
      scale_color_manual(values=c("gray", "red")) + 
      l1000.analysis::theme_scatterplot() + 
      theme(legend.title=element_blank()) + 
      theme(legend.text=element_text(size=30))
    
    
  } else {
    
    comparison_df %>% ggplot(aes(logFC_cpd, logFC_shrna)) + 
      geom_point() + 
      cdsr::theme_Publication() + 
      geom_text_repel(data = unique(rbind(
        comparison_df %>% arrange(-abs(logFC_cpd)) %>% head(10),
        comparison_df %>% arrange(-abs(logFC_shrna)) %>% head(10))),
        aes(label = Gene),
        size = 10) + 
      labs(x = cpd_name, y = shrna_name) + 
      l1000.analysis::theme_scatterplot()
  }
  
  
}


#' Plot comparison of two compounds' viability signatures
#'
#' @param data_list
#' @param via_sig_list
#' @param cpd_name
#' @keywords
#' @export
#' @examples
#' plot_pert_data_heatmap()
plot_pert_data_heatmap <- function(data_list, via_sig_list, cpd_name, num_genes = 30, additional_metadata = NULL, sens_columns = c("auc_avg", "auc_prism", "auc_gdsc", "auc_ctrp"), sort_column = "auc_avg", remove_via_nas = F, ...) {
  
  full_df <- cbind(data_list$expression, data_list$metadata, data_list$sensitivity)
  
  cpd_df <- full_df[full_df$pert_iname == cpd_name,]
  cpd_df_collapsed <- cpd_df %>% 
    group_by(pert_iname, cell_id) %>% 
    summarise_all(mean) %>% 
    as.data.frame()
  
  cpd_via_sig <- via_sig_list[[cpd_name]][['viability_related']]
  
  top_via_sigs <- cpd_via_sig %>% 
    arrange(-abs(t)) %>% 
    head(num_genes)
  top_via_sig_genes <- top_via_sigs$Gene
  
  cpd_df_collapsed %<>% 
    set_rownames(cpd_df_collapsed$cell_id)
  sorted_idx <- order(cpd_df_collapsed[,sort_column])
  annot_df <- cpd_df_collapsed[sorted_idx, sens_columns, drop = F]
  annot_df <- annot_df[,colSums(is.na(annot_df)) != nrow(annot_df), drop = F]
  
  if (!is.null(additional_metadata)) {
    annot_df <- merge(annot_df, additional_metadata, by.x = "row.names", by.y = "cell_id", all.x = T)
    annot_df <- annot_df %>% column_to_rownames("Row.names")
  }
  
  heatmap_data_df <- cpd_df_collapsed[sorted_idx, top_via_sig_genes]
  
  if (remove_via_nas) {
    na_idx <- which(is.na(annot_df[,sort_column]))
    heatmap_data_df <- heatmap_data_df[which(!is.na(annot_df[,sort_column])),]
    annot_df <- annot_df[which(!is.na(annot_df[,sort_column])),,drop = F]
  }
  
  pheatmap::pheatmap(heatmap_data_df %>% scale() %>% t(), annotation_col = annot_df, cluster_cols = F, ...)
  
}






