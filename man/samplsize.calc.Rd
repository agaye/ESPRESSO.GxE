\name{samplsize.calc}
\alias{samplsize.calc}
\title{Calculates the sample size required to achieve the desired power}
\usage{
  samplsize.calc(numcases = 2000, numcontrols = 8000,
    num.subjects = 500, pheno.model = 0, pval = 1e-04,
    power = 0.8, mean.model.z = NULL)
}
\arguments{
  \item{numcases}{Number of cases when outcome is binary}

  \item{numcontrols}{Number of controls when outcome is
  binary}

  \item{num.subjects}{Number of subjects when outcome is
  continuous}

  \item{pheno.model}{Outcome type, 0 for binary and 1 for
  continuous}

  \item{pval}{Cut-off p-value defining statistical
  significance}

  \item{power}{Desired power}

  \item{mean.model.z}{Ratio of mean beta estimate over mean
  se estimate}
}
\value{
  A table containing: \code{numcases.required} Number of
  cases required to achieve the desired power under binary
  outcome model \code{numcontrols.required} Number of
  controls required to achieve the desired power under
  binary outcome model \code{numsubjects.required} Number
  of subjects required to achieve the desired power under a
  quantatative outcome model
}
\description{
  Estimates by how much the simulated study size needs to
  be inflated or shrank in order to obtain the specified
  level of power. The ratio of z-statistic required for
  desired power to mean model z-statistic obtained
  indicates the relative changes in standard error
  required. This cirresponds to relative change on scale
  and square root of sample size.
}
\author{
  Amadou Gaye
}

