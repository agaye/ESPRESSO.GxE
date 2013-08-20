\name{pheno.params}
\alias{pheno.params}
\docType{data}
\title{Parameters for the outcome variable}
\description{A table with 6 columns: scenario.id, pheno.model, disease.prev, pheno.sensitivity, pheno.specificity, pheno.reliability}
\usage{data(pheno.params)}
\format{
  A table with 6 columns and 1 row for each scenario
  \describe{
    \item{\code{scenario.id}}{The ID of the scenario}
    \item{\code{pheno.model}}{Type of outcome; 0 for binary and 1 for continuous}
    \item{\code{disease.prev}}{Prevalence of the binary outcome}
    \item{\code{pheno.sensitivity}}{Sensitivity of the assessment of binary outcome}
    \item{\code{pheno.specificity}}{Specificity of the assessment of binary outcome}
    \item{\code{pheno.reliability}}{Reliability of the assessment of continuous outcome}
  }
}

\examples{data(pheno.params)}