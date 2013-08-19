\name{get.observed.data}
\alias{get.observed.data}
\title{Generates exposure data with some error}
\usage{
  get.observed.data(true.data = NULL,
    geno.error = c(0.05, 0.05), geno.model = 0, MAF = 0.1,
    env.error = c(0.15, 0.15), env.model = 0,
    env.prev = 0.1, env.sd = 1, env.reliability = 0.8,
    pheno.model = 0, pheno.error = c(0.1, 0.1),
    pheno.reliability = 0.9)
}
\arguments{
  \item{true.data}{Input table of simulated data considered
  as true data}

  \item{geno.error}{Misclassification rates in genetic
  assessment: 1-sensitivity and 1-specificity}

  \item{geno.model}{Genetic model; 0 for binary and 1 for
  continuous}

  \item{MAF}{Minor allele frequency}

  \item{env.error}{Misclassification rates in environmental
  exposures assessment: 1-sensitivity and 1-specificity}

  \item{env.model}{Model of the exposure: binary=0,
  quantitative-normal=1 or quantitative-uniform=2}

  \item{env.prev}{Prevalence of environmental exposure}

  \item{env.sd}{Standard Deviation}

  \item{env.reliability}{Reliability of the assessment of
  quantitative exposure}

  \item{pheno.model}{Model of the phenotype: binary=0
  ,quantitative-normal=1 or quantitative-uniform=2}

  \item{pheno.error}{Misclassification rate: 1 to 0, 0 to
  1}

  \item{pheno.reliability}{Reliability of the assessment of
  quantitative phenotype}
}
\value{
  A matrix
}
\description{
  Uses functions make.obs.geno and make.obs.env to generate
  effect data with a set level of error
}
\author{
  Amadou Gaye
}

