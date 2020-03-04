library(magrittr)


#' Quantile normalize a dataset so that each column has the same distribution (typically used for normalizing across sensitivity datasets)
#'
#' @param data_df dataframe to be normalized
#' 
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' 
#' 
#' @keywords
#' @export
#' @examples
#' quantile_normalize_viability_data()
quantile_normalize_viability_data <- function(data_df) {
  
  # get samples with full data across all features (columns)
  na_counts <- rowSums(!is.na(data_df))
  all_nonna_idx <- which(na_counts == ncol(data_df))
  
  # fit a normalization distribution to this dense matrix
  normalization_dist <- preprocessCore::normalize.quantiles.determine.target(x = data_df[all_nonna_idx,] %>% as.matrix())
  
  # normalize the entire dataset by this distribution
  normalized_df <- preprocessCore::normalize.quantiles.use.target(x = data_df %>% as.matrix(), target = normalization_dist)
  
  return(normalized_df)
}
