\name{sim.CC.data}
\alias{sim.CC.data}
\title{Simulates case and controls}
\usage{
  sim.CC.data(num.obs = 20000, numcases = 2000,
    numcontrols = 8000, allowed.sample.size = 2e+07,
    disease.prev = 0.1, MAF = 0.1, geno.model = 0,
    geno.OR = 1.5, env.model = 0, env.prev = 0.1,
    env.mean = 0, env.sd = 1, env.low.lim = 0,
    env.up.lim = 1, env.OR = 1.5, int.OR = 1.5,
    baseline.OR = 12.36)
}
\arguments{
  \item{num.obs}{Number of observations to generate per
  iteration}

  \item{numcases}{Number of cases to simulate}

  \item{numcontrols}{Number of controls to simulate}

  \item{allowed.sample.size}{Maximum number of observations
  allowed}

  \item{disease.prev}{Prevalence of the binary outcome}

  \item{MAF}{Minor allele frequency}

  \item{geno.model}{Genetic model; 0 for binary and 1 for
  continuous}

  \item{geno.OR}{Odds ratios of the genetic determinants}

  \item{env.model}{Model of the environmental exposure}

  \item{env.prev}{Prevelance of the environmental
  determinates}

  \item{env.mean}{Mean under quantitative-normal model}

  \item{env.sd}{Standard deviation under
  quantitative-normal model}

  \item{env.low.lim}{Lower limit under quantitative-uniform
  model}

  \item{env.up.lim}{Upper limit under quantitative-uniform
  model}

  \item{env.OR}{Odds ratios of the environmental
  determinants}

  \item{int.OR}{Odds ration of the interaction}

  \item{baseline.OR}{Baseline odds ratio for subject on 95
  percent population centile versus 5 percentile. This
  parameter reflects the heterogeneity in disease risk
  arising from determinates that have not been measured or
  have not been included in the model}
}
\value{
  A matrix
}
\description{
  Generates affected and non-affected subjects until the
  set sample size is achieved
}
\author{
  Amadou Gaye
}

