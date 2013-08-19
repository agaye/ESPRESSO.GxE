\name{sim.pheno.bin}
\alias{sim.pheno.bin}
\title{Generates phenotype status}
\usage{
  sim.pheno.bin(num.obs = 20000, disease.prev = 0.1,
    genotype = NULL, environment = NULL,
    interaction = NULL, subject.effect.data = NULL,
    geno.OR = 1.5, env.OR = 1.5, int.OR = 1.5)
}
\arguments{
  \item{num.obs}{Number of observations to generate per
  iteration}

  \item{disease.prev}{Prevalence of the binary outcome}

  \item{genotype}{Exposure data for genetic determinate}

  \item{environment}{Exposure data for environment}

  \item{interaction}{<-----------------------------------------------------------------------}

  \item{subject.effect.data}{Subject effect data, reflects
  the heterogenity in baseline disease risk}

  \item{geno.OR}{Odds ratios of the two genetic
  determinates}

  \item{env.OR}{Odds ratios of the two environments}

  \item{int.OR}{<-----------------------------------------------------------------------}
}
\value{
  A dataframe of phenotype
}
\description{
  Generates affected and non-affected subjects
}
\author{
  Amadou Gaye
}

