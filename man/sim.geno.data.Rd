\name{sim.geno.data}
\alias{sim.geno.data}
\title{Simulates genotypes for a genetic variant}
\usage{
  sim.geno.data(num.obs = 20000, MAF = 0.1, geno.model = 0)
}
\arguments{
  \item{num.obs}{Number of observations to simulate}

  \item{MAF}{Minor allele frequency of the variant}

  \item{geno.model}{Genetic model; 0 for binary and 1 for
  continuous}
}
\value{
  A dataframe containing the following variables:
  \code{allele.A} Major allele \code{allele.B} Minor allele
  \code{genotype} Genotype
}
\description{
  Generates two alleles and combines them to form the
  genotype of a SNP under a binary or additive genetic
  model
}

