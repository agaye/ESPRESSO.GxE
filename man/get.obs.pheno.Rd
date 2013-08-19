\name{get.obs.pheno}
\alias{get.obs.pheno}
\title{Generates some observed environmental exposure data}
\usage{
  get.obs.pheno(phenotype = NULL, pheno.model,
    pheno.error.1.0 = 0.05, pheno.error.0.1 = 0.05,
    pheno.reliability = 0.9)
}
\arguments{
  \item{phenotype}{phenotype}

  \item{pheno.model}{Model of the phenotype: binary=0
  ,quantitative-normal=1 or quantitative-uniform=2}

  \item{pheno.error.1.0}{1 to 0 misclassification rate}

  \item{pheno.error.0.1}{0 to 1 misclassification rate}

  \item{pheno.reliability}{Reliability of the assessment of
  quantitative phenotype}
}
\value{
  A dataframe containing: \code{observed.phenotype}
  Observed data
}
\description{
  Adds a set level of error to simulated binary or
  quantitative data(the true data) to obtain data with a
  larger variance (the observed data).
}
\author{
  Amadou Gaye
}

