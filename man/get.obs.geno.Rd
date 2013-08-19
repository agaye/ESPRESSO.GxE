\name{get.obs.geno}
\alias{get.obs.geno}
\title{Adds some error to genotype data}
\usage{
  get.obs.geno(allele.A = NULL, allele.B = NULL,
    geno.model = 0, MAF = 0.1, geno.error.0.1 = 0.05,
    geno.error.1.0 = 0.05)
}
\arguments{
  \item{allele.A}{Allele A}

  \item{allele.B}{Allele B}

  \item{geno.model}{Genetic model; 0 for binary and 1 for
  continuous}

  \item{MAF}{Minor Allele Frequency}

  \item{geno.error.0.1}{0 to 1 misclassififcation rate}

  \item{geno.error.1.0}{1 to 0 misclassififcation rate}
}
\value{
  A dataframe containing: \code{observed.genotype} Observed
  genotypes \code{observed.allele.A} Observed A alleles
  \code{observed.allele.B} Observed B alleles
}
\description{
  Simulates errors and adds it to the true data to obtain
  observed data. The alleles simulated by the function
  sim.geno.data are randomly misclassified and used to form
  new genotypes that represent the observed genotypes.
}
\author{
  Amadou Gaye
}

