\name{geno.params}
\alias{geno.params}
\docType{data}
\title{Parameters to simulate the genetic determinant}
\description{A table with 7 columns: scenario.id, geno.model, MAF, geno.OR, geno.efkt, geno.sensitivity, geno.specificity.}
\usage{data(geno.params)}
\format{
  A table with 7 columns and 1 row for each scenario
  \describe{
    \item{\code{scenario.id}}{The ID of the scenario}
    \item{\code{geno.model}}{Genetic model; 0 for binary and 1 for continuous}
    \item{\code{MAF}}{Minor allele frequency}
    \item{\code{geno.OR}}{Odds ratios of the genetic determinants}
    \item{\code{geno.efkt}}{Effects of the genetic variants}
    \item{\code{geno.sensitivity}}{Sensitivity of the assessment to the genetic determinants}
    \item{\code{geno.specificity}}{Specificity of the assessment to the genetic determinants}
  }
}
\examples{data(geno.params)}