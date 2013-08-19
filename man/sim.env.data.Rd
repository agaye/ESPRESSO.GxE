\name{sim.env.data}
\alias{sim.env.data}
\title{Simulates cases and controls}
\usage{
  sim.env.data(num.obs = 20000, env.model = 0,
    env.prev = 0.1, env.mean = 0, env.sd = 1,
    env.low.lim = 0, env.up.lim = 1)
}
\arguments{
  \item{num.obs}{Number of observations to simulate}

  \item{env.model}{Model of the exposure: binary=0,
  quantitative-normal=1, quantitative-uniform=2}

  \item{env.prev}{Prevalence of the environmental exposure}

  \item{env.mean}{Mean under quantitative-normal model}

  \item{env.sd}{Standard deviation under
  quantitative-normal model}

  \item{env.low.lim}{Lower limit under quantitative-uniform
  model}

  \item{env.up.lim}{Upper limit under quantitative-uniform
  model}
}
\value{
  A vector of continuous or binary values
}
\description{
  Generates data for a binary, quantitative-normal, or
  quantitative-uniform environmental determinate
}
\author{
  Amadou Gaye
}

