% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viability_signature_modeling.R
\name{compute_viability_signature_all_cpds_control_dose_and_time}
\alias{compute_viability_signature_all_cpds_control_dose_and_time}
\title{Compute a unique viability signature for each drug, controlling for dose}
\usage{
compute_viability_signature_all_cpds_control_dose_and_time(data_list,
  sens_column = "auc_avg", selectivity_threshold = 0,
  covariates = "")
}
\arguments{
\item{sens_column}{}
}
\description{
Compute a unique viability signature for each drug, controlling for dose
}
\examples{
compute_viability_signature_all_cpds_control_dose_and_time()
}
