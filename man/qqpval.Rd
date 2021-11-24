\name{qqpval}
\alias{qqpval}
\title{Functions for inspecting population substructure}
\description{
  This function plots ranked observed p values against the corresponding expected p values in -log scale. 
}
\usage{
qqpval(p, pch=16, col=4, \dots)
}

\arguments{
  \item{p}{a vector of p values}
  \item{pch}{symbol to use for points}
  \item{col}{color for points}
  \item{\dots}{other plot arguments}
}


\value{
  No return value, just the plot
}


\seealso{ \code{\link{GenomicControl}}, \code{\link{WGassociation}}}

\examples{
data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")
res<-WGassociation(casco,datSNP,model=c("do","re","log-add"))

# observed vs expected p values for recessive model
qqpval(recessive(res))

}

\keyword{utilities}
