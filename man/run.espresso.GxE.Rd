\name{run.espresso.GxE}
\alias{run.espresso.GxE}
\title{Runs a full ESPRESSO analysis}
\usage{
  run.espresso.GxE(simulation.params = NULL,
    pheno.params = NULL, geno.params = NULL,
    env.params = NULL, scenarios2run = 1)
}
\arguments{
  \item{simulation.params}{general parameters for the
  scenario(s) to analyse}

  \item{pheno.params}{paramaters for the outcome variables}

  \item{geno.params}{parameters for the genetic
  determinant}

  \item{env.params}{parameters for the environmental
  determinant}

  \item{scenarios2run}{the indices of the scenarios one
  wish to analyse}
}
\value{
  a summary table that contains both the input parameters
  and the results of the analysis
}
\description{
  This function calls the functions required to run a full
  ESPRESSO analysis where the model consit of an outcome
  (binary or continuous) determinant by a binary or
  continuous environmental determinant.
}
\author{
  Amadou Gaye
}

