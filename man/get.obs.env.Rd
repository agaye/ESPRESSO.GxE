\name{get.obs.env}
\alias{get.obs.env}
\title{Generates some observed environmental exposure data}
\usage{
  get.obs.env(env.data = NULL, env.model = 0,
    env.prev = 0.1, env.sd = 1, env.error = c(0.15, 0.15),
    env.reliability = 0.8)
}
\arguments{
  \item{env.data}{A vector of environmental measures that
  represents the true data}

  \item{env.model}{Model of the exposure: binary=0
  ,quantitative-normal=1 or quantitative-uniform=2}

  \item{env.prev}{Prevalence of the environmental
  determinat}

  \item{env.sd}{Standard deviation}

  \item{env.error}{misclassification rates: 1-sensitivity
  and 1-specificity}

  \item{env.reliability}{Reliability of the assessment of
  quantitative exposure}
}
\value{
  A dataframe of two coloumns: \code{true.environment}
  Input data (true data) \code{observed.environment}
  Observed data
}
\description{
  Adds a set level of error to simulated binary or
  quantitative data(the true data) to obtain data with a
  larger variance (the observed data). The level of error
  is determined by the misclassification rates in binary
  data and by the set level of variance in the quantitative
  data.
}
\author{
  Amadou Gaye
}

