\name{sim.QTL.data.GxE}
\alias{sim.QTL.data.GxE}
\title{Simulates subjects for continuous outcome}
\usage{
sim.QTL.data.GxE(numsubjects = 500, ph.mean = 0, ph.sd = 1, MAF = 0.1,
  geno.model = 0, geno.efkt = 0.25, env.model = 0, env.efkt = 0.5,
  env.prev = 0.1, env.mean = 0, env.sd = 1, env.low.lim = 0,
  env.up.lim = 1, int.efkt = 1.5, pheno.reliability = 0.9)
}
\arguments{
  \item{num.subjects}{Number of subjects to simulate}

  \item{ph.mean}{statistical mean}

  \item{ph.sd}{standard deviation}

  \item{MAF}{Minor allele frequencies of the two genetic
  variants}

  \item{geno.model}{Genetic model; 0 for binary and 1 for
  continuous}

  \item{geno.efkt}{Effects of the genetic variants}

  \item{env.model}{Model of the environmental exposure}

  \item{env.efkt}{Effects of the environment determinats}

  \item{env.prev}{Prevalences of the environmental
  determinats}

  \item{env.mean}{Mean under quantitative-normal model}

  \item{env.sd}{Standard deviation under
  quantitative-normal model}

  \item{env.low.lim}{lower limit under quantitative-uniform
  model}

  \item{env.up.lim}{upper limit under quantitative-uniform
  model}

  \item{int.efkt}{Interaction effect}

  \item{pheno.reliability}{reliability of the assessment
  for a quantitative outcome.}
}
\value{
A matrix
}
\description{
Generates the specified number of subjects
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}
