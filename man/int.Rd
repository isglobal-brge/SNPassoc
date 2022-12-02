\name{int}
\alias{int}
\title{
Identify interaction term
}
\description{
This is a special function used for 'haplo.interaction' function.  It
identifies the variable that will interact with the haplotype estimates. 
Using \code{int()} in a formula implies that
the interaction term between this variable and haplotypes is included in 'haplo.glm' function.
}
\usage{
int(x)
}

\arguments{
\item{x}{
A factor variable.
}}

\value{
\code{x}
}

\seealso{
\code{\link{haplo.interaction}}
}

\examples{

library(SNPassoc)
library(haplo.stats)

data(SNPs)
datSNP<-setupSNP(SNPs, 6:40, sep = "")
mod <- haplo.interaction(casco~int(sex)+blood.pre, data = datSNP, 
                         SNPs.sel = c("snp10001","snp10004","snp10005"))

}

\keyword{utilities}