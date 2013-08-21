\name{sim.pheno.qtl}
\alias{sim.pheno.qtl}
\title{Simulates continuous outcome data}
\usage{
  sim.pheno.qtl(numsubjects = 10000, genotype = NULL,
    geno.efkt = 0.25, environment = NULL, env.efkt = 0.5,
    interaction = NULL, int.efkt = 1.5)
}
\arguments{
  \item{num.subjects}{Number of subjects to simulate}

  \item{genotype}{Genotype}

  \item{geno.efkt}{Effects if the genetic variants}

  \item{environment}{Exposure data for environment}

  \item{env.efkt}{Effects of the environmental determiants}

  \item{interaction}{Effect model: main effects=0,
  Gene-Environment interaction=1, Gene-Gene interaction=2
  and Environment-Enviroment interaction=3}

  \item{int.efkt}{Interaction effect}
}
\value{
  A dataframe of phenotype
}
\description{
  Uses the effects data of the determinants to construct a
  linear predictor(LP). The outcome is normally distributed
  variable generated with a mean equal to LP and a standard
  deviation of 1. Some error is then added to the simulated
  outcome to obtained the observed outcome.
}
\author{
  Amadou Gaye
}

