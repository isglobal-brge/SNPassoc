\name{dscore}
\alias{dscore}
\alias{dscore.default}
\alias{dscore.character}
\alias{plot.dscore}
\alias{pscore}


\title{Exact distribution of a genetic score}
\description{
  Computes the exact distribution of a genetic score on a set of SNPs
}
\usage{
dscore(x, \dots)

\method{dscore}{default}(x, \dots)

\method{dscore}{character}(x, \dots)

pscore(q, x, \dots)

}
\arguments{
  \item{x}{a numeric vector containing either the minor allele frequencies (MAF) of the SNPs used to created the genetic score or the names of the SNPs.}
  \item{q}{quantile}
  \item{\dots}{further arguments to be passed to useMart function from biomaRt (not necessary).}
}
\references{
I Subirana, G Lucas, C Lluis-Granella, R Elosua, JR Gonzalez. Predicting individual disease probability using the exact distribution of a genetic risk score at the population level. Submitted   
}
\value{
  A vector with the probabilities of having 0, 1, 2, 3, ... 2*K alleles of a genetic score (K indicates the number of SNPs). When name of the SNPs is provided, an attribute called 'MAFs' is also returned 
}

\examples{
gs <- dscore(c(0.23, 0.16, 0.12, 0.02))
gs 
plot(gs)

}
\keyword{utilities}