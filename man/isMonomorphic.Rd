\name{is.Monomorphic}
\alias{is.Monomorphic}
\title{Check whether a SNP is Monomorphic}
\description{
 This function verifies when a SNP is Monomorphic
}
\usage{
is.Monomorphic(x) 
}
\arguments{
  \item{x}{any R object}
}
\value{
A logical value TRUE if the SNP is Monomorphic, otherwise a FALSE
}

\examples{

data(SNPs)
is.Monomorphic(SNPs$snp10001)
is.Monomorphic(SNPs$snp100020)
apply(SNPs[,20:30],2,is.Monomorphic)
}

\keyword{internal}
