% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_pert_data_heatmap}
\alias{plot_pert_data_heatmap}
\title{Plot comparison of two compounds' viability signatures}
\usage{
plot_pert_data_heatmap(data_list, via_sig_list, cpd_name, num_genes = 30,
  additional_metadata = NULL, sens_columns = c("auc_avg", "auc_prism",
  "auc_gdsc", "auc_ctrp"), sort_column = "auc_avg", remove_via_nas = F,
  ...)
}
\description{
Plot comparison of two compounds' viability signatures
}
\examples{
plot_pert_data_heatmap()
}
