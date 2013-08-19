\name{power.calc}
\alias{power.calc}
\title{Calculates the empirical and theoretical power}
\usage{
  power.calc(pval = 1e-04, z.values = NULL,
    mean.model.z = NULL)
}
\arguments{
  \item{pval}{Cut-off p-value defining statistical
  significance}

  \item{z.values}{z-statistic of the determinant}

  \item{mean.model.z}{Mean z-statistic of the environmental
  determinant}
}
\value{
  A table containing empirical power and theoretical power
}
\description{
  Determines the empirical and theoretical power. The
  empirical power is the proportion of simulations in which
  the z-statistic for the parameter of interest exceeds the
  z-statistic for the desured level if statistical
  significance. Theoretical power is the power of the
  study.
}
\author{
  Amadou Gaye
}

