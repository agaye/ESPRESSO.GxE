\name{env.params}
\alias{env.params}
\docType{data}
\title{Parameters to simulate environmental exposures data}
\description{A table with 12 columns: scenario.id, env.model, env.prevalence, env.OR, env.efkt, env.mean, env.sd, env.low.lim, env.up.lim, env.sensitivity, env.specificity, env.reliability}
\usage{data(env.params)}
\format{
   A table with 12 columns and 1 row for each scenario
  \describe{
    \item{\code{scenario.id}}{The ID of the scenario}
    \item{\code{env.model}}{Model of the environmental exposure; 0 for binary and 1 for continuous}
    \item{\code{env.prevalence}}{Prevelance of environment}
    \item{\code{env.OR}}{Odds ratio of environmental exposure}
    \item{\code{env.efkt}}{Effect size of environment}
    \item{\code{env.mean}}{The mean measure for the environment under quantitative-normal model}
    \item{\code{env.sd}}{Standard deviation under quantitative-normal model}
    \item{\code{env.low.lim}}{Lower limit under quantitative-uniform model}
    \item{\code{env.up.lim}}{Upper limit under quantitative-uniform model}
    \item{\code{env.sensitivity}}{Sensitivity of the assessment to the environment}
    \item{\code{env.specificity}}{Specificity of the assessment to the environment}
    \item{\code{env.reliability}}{Reliability of the assessment of quantative exposure}
  }
}

\examples{data(env.params)}