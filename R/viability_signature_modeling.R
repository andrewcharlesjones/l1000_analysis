library(limma)
library(tibble)



#' Compute viability signature
#'
#' @param data_list
#' @param compound
#' @param global
#' @param sens_column
#' 
#' @keywords
#' @export
#' @examples
#' compute_viability_signature()
compute_viability_signature <- function(data_list, compound = NA, global = F, sens_column = "auc_avg") {
  
  if (global == T & !is.na(compound)) {
    stop("Cannot specify global and a specific compound. Use one or the other, idiot.")
  }
  
  if (global == T) {
    cpd_idx <- 1:nrow(data_list$expression)
  } else if (!is.na(compound)) {
    cpd_idx <- which(data_list$metadata$pert_iname == compound)
  } else {
    stop("Need to either specify global or specific compound.")
  }

  cpd_data_with_meta <- cbind(data_list$expression[cpd_idx,], 
                              data_list$sensitivity[cpd_idx,,drop = F], 
                              data_list$metadata[cpd_idx,])
  
  cpd_data_collapsed <- cpd_data_with_meta %>% 
    group_by(pert_iname, cell_id) %>% 
    summarise_all(mean)
  
  via_sig_list <- l1000.analysis::run_limma_via_sig(cpd_data_collapsed, gene_names = colnames(data_list$expression), sens_column = sens_column)
  
  return(via_sig_list)
  
}



#' Run limma linear model to compute viability signature
#'
#' @param data_df
#' @param sens_column
#' 
#' @keywords
#' @export
#' @examples
#' run_limma_via_sig()
run_limma_via_sig <- function(data_df, gene_names, sens_column = "auc_avg", covariates = "") {
  
  
  
  cpd_data_collapsed <- data_df[!is.na(data_df[,sens_column]),]
  
  if (covariates != "") {
    sample_info <- cpd_data_collapsed[,c("pert_iname", "cell_id", sens_column, covariates)]
  } else {
    sample_info <- cpd_data_collapsed[,c("pert_iname", "cell_id", sens_column)]
  }
  
  
  if (covariates != "") {
    design <- model.matrix(as.formula(
      paste(
        paste(" ~ ", sens_column, sep = ""),
        covariates, sep = "+")),
      data = sample_info)
  } else {
    design <- model.matrix(as.formula(paste(" ~ ", sens_column)),
      data = sample_info)
  }
  
  
  fit <- lmFit(cpd_data_collapsed[,gene_names] %>% t() %>% as.data.frame(), design = design)
  fit <- eBayes(fit)
  
  viability_unrelated <- topTable(fit, number = Inf, coef = 1) %>% rownames_to_column("Gene")
  viability_related <- topTable(fit, number = Inf, coef = 2) %>% rownames_to_column("Gene")
  
  return(
    list(
      viability_unrelated = viability_unrelated,
      viability_related = viability_related
    )
  )
}




#' Compute a unique viability signature for each drug
#'
#' @param data_list
#' @param sens_column
#' 
#' @keywords
#' @export
#' @examples
#' compute_viability_signature_all_cpds()
compute_viability_signature_all_cpds <- function(data_list, sens_column = "auc_avg", selectivity_threshold = 0.0, covariates = "") {
  
  full_df <- cbind(data_list$expression, data_list$metadata, data_list$sensitivity)
  
  full_df_collapsed <- full_df %>% 
    group_by(pert_iname, cell_id) %>% 
    summarise_all(mean)
  
  pert_counts <- full_df_collapsed$pert_iname %>% table() %>% as.data.frame() %>% set_colnames(c("pert_iname", "count"))
  usable_perts <- pert_counts$pert_iname[pert_counts$count >= 3]
  
  full_df_collapsed <- full_df_collapsed[full_df_collapsed$pert_iname %in% usable_perts,]
  
  
  all_cpd_via_sigs <- dlply(full_df_collapsed, ~ pert_iname, .progress = "text", function(x) {
    
    via_sig_list <- l1000.analysis::run_limma_via_sig(x, gene_names = colnames(data_list$expression), sens_column = sens_column, covariates = covariates)
    
    return(via_sig_list)
    
  })
  
  compute_sensitivity_selectivities(data_list = data_list, sens_column = sens_column)
  
  return(all_cpd_via_sigs)
  
}

#' Compute a unique viability signature for each drug, controlling for dose
#'
#' @param data_list
#' @param sens_column
#' 
#' @keywords
#' @export
#' @examples
#' compute_viability_signature_all_cpds_control_dose_and_time()
compute_viability_signature_all_cpds_control_dose_and_time <- function(data_list, sens_column = "auc_avg", selectivity_threshold = 0.0, covariates = "") {
  
  full_df <- cbind(data_list$expression, data_list$metadata, data_list$sensitivity)
  
  full_df_collapsed <- full_df %>% 
    group_by(pert_iname, cell_id, pert_dose, pert_time) %>% 
    summarise_all(mean)
  
  pert_counts <- full_df_collapsed$pert_iname %>% table() %>% as.data.frame() %>% set_colnames(c("pert_iname", "count"))
  usable_perts <- pert_counts$pert_iname[pert_counts$count >= 3]
  
  full_df_collapsed <- full_df_collapsed[full_df_collapsed$pert_iname %in% usable_perts,]
  
  
  all_cpd_via_sigs <- dlply(full_df_collapsed, ~ pert_iname, .progress = "text", function(x) {
    
    
    is_single_dose <- length(x$pert_dose %>% unique()) == 1
    is_single_time <- length(x$pert_time %>% unique()) == 1
    if (is_single_dose) {
      if (is_single_time) {
        covariates <- c("pert_dose", "pert_time")
      } else {
        covariates <- c("pert_dose")
      }
    }
    if (is_single_time) {
      covariates <- c("pert_time")
    } else {
      covariates <- c("")
    }
    
    via_sig_list <- l1000.analysis::run_limma_via_sig(x, gene_names = colnames(data_list$expression), sens_column = sens_column, covariates = covariates)
    
    return(via_sig_list)
    
  })
  
  compute_sensitivity_selectivities(data_list = data_list, sens_column = sens_column)
  
  return(all_cpd_via_sigs)
  
}


#' Get a matrix of viability sigantures given a list of them
#'
#' @param via_sig_list
#' @param which_sig
#' 
#' @keywords
#' @export
#' @examples
#' get_cpd_viability_signature_matrix()
get_cpd_viability_signature_matrix <- function(via_sig_list, which_sig = "viability_related") {
  
  via_sig_mat_list <- list()
  for (curr_cpd in names(via_sig_list)) {
    curr_sig <- via_sig_list[[curr_cpd]][[which_sig]]
    via_sig_mat_list <- rbind(via_sig_mat_list, curr_sig)
  }
  via_sig_df <- via_sig_mat_list %>% as.matrix() %>% as.data.frame() %>% 
    set_colnames(names(via_sig_list)) %>% 
    set_rownames(curr_sig$Gene)
  
  return(via_sig_df)
  
}



#' Compute a unique viability signature for each drug using Pearson correlation (other function uses limma)
#'
#' @param data_list
#' @param sens_column
#' 
#' @keywords
#' @export
#' @examples
#' compute_viability_signature_pearson_all_cpds()
compute_viability_signature_pearson_all_cpds <- function(data_list, sens_column = "auc_avg") {
  
  full_df <- cbind(data_list$expression, data_list$metadata, data_list$sensitivity)
  
  full_df_collapsed <- full_df %>% 
    group_by(pert_iname, cell_id) %>% 
    summarise_all(mean)
  
  all_cpd_via_sigs <- dlply(full_df_collapsed, ~ pert_iname, .progress = "text", function(x) {
    
    curr_via_sig <- cor(x[,colnames(data_list$expression)], x[,sens_column]) %>% 
      as.data.frame() %>% 
      set_colnames("via_sig_coeff") %>% 
      mutate(pert_iname = rep(x$pert_iname[1]), ncol(data_list$expression)) %>% 
      rownames_to_column("Gene")
    
    return(curr_via_sig)
    
  })
  
  return(all_cpd_via_sigs)
  
}



#' Compute correlations between TAS and sensitivity
#'
#' @param data_list
#' @param sens_column
#' 
#' @keywords
#' @export
#' @examples
#' compute_tas_viability_cors()
compute_tas_viability_cors <- function(data_list, sens_column = "auc_avg") {
  
  full_df <- cbind(data_list$expression, data_list$metadata, data_list$sensitivity)
  
  full_df_collapsed <- full_df %>% 
    group_by(pert_iname, cell_id) %>% 
    summarise_all(mean)
  
  all_cpd_via_sigs <- dlply(full_df_collapsed, ~ pert_iname, .progress = "text", function(x) {
    
    curr_tas_via_cor <- cor(x$tas_l2, x[,sens_column])
    
    return(curr_tas_via_cor)
    
  })
  
  return(all_cpd_via_sigs)
  
}


#' Compute the strength of each compound's viability signature
#'
#' @param via_sig_list
#' @param which_sig
#' 
#' @keywords
#' @export
#' @examples
#' compute_viability_signature_strengths()
compute_viability_signature_strengths <- function(via_sig_list, which_sig = "viability_related") {
  
  sig_strength_list <- ldply(via_sig_list, function(x) {
    
    # curr_strength <- sqrt(sum(x[['viability_related']][['logFC']] ^ 2))
    curr_strength <- x[['viability_related']][['t']] %>% 
      abs() %>% 
      mean()
    return(curr_strength)
    
  }) %>% 
    set_colnames(c("pert_iname", "via_sig_strength"))
  
  return(sig_strength_list)
  
}


#' Compute selectivity of the killing of each compound
#'
#' @param data_list
#' @param sens_column
#' 
#' @keywords
#' @export
#' @examples
#' compute_sensitivity_selectivities()
compute_sensitivity_selectivities <- function(data_list, sens_column = "auc_avg") {
  
  full_df <- cbind(data_list$expression, data_list$metadata, data_list$sensitivity)
  
  full_df_collapsed <- full_df[,c("pert_iname", "cell_id", sens_column)] %>% 
    group_by(pert_iname, cell_id) %>% 
    summarise_all(mean)
  
  selectivities <- full_df_collapsed %>% 
    group_by(pert_iname) %>% 
    summarise_all(sd, na.rm = T) %>% 
    as.data.frame() %>% 
    dplyr::select(pert_iname, sens_column) %>% 
    set_colnames(c("pert_iname", "selectivity"))
  
  return(selectivities)
  
}

#' Compute average sensitivity of each compound
#'
#' @param data_list
#' @param sens_column
#' 
#' @keywords
#' @export
#' @examples
#' compute_avg_sensitivities()
compute_avg_sensitivities <- function(data_list, sens_column = "auc_avg") {
  
  full_df <- cbind(data_list$expression, data_list$metadata, data_list$sensitivity)
  
  full_df_collapsed <- full_df[,c("pert_iname", "cell_id", "auc_avg")] %>% 
    group_by(pert_iname, cell_id) %>% 
    summarise_all(mean)
  
  avg_sens <- full_df_collapsed %>% 
    group_by(pert_iname) %>% 
    summarise_all(mean) %>% 
    as.data.frame() %>% 
    dplyr::select(pert_iname, auc_avg) %>% 
    set_colnames(c("pert_iname", "mean_sensitivity"))
  
  return(avg_sens)
  
}







