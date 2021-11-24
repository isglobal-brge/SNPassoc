\name{getSignificantSNPs}
\alias{getSignificantSNPs}

\title{Extract significant SNPs from an object of class 'WGassociation' }

\description{
  Extract significant SNPs from an object of class 'WGassociation' when genomic 
  information is available}

\usage{
getSignificantSNPs(x, chromosome, model, sig = 1e-15)
}

\arguments{
  \item{x}{an object of class 'WGassociation'}
  \item{chromosome}{ chromosome from which SNPs are extracted}
  \item{model}{ genetic model from which SNPs are extracted}
  \item{sig}{ statistical significance level. The default is 1e-15}
}

\value{
  A list with the following components:
  \item{names}{the name of SNPs}
  \item{column}{the columns corresponding to the SNPs in the original data frame}
  ...
}

\seealso{  \code{\link{WGassociation}} }

\examples{
data(resHapMap)
# resHapMap contains the results for a log-additive genetic model

# to get the significant SNPs for chromosome 12
getSignificantSNPs(resHapMap,chromosome=12)
# to get the significant SNPs for chromosome 5
getSignificantSNPs(resHapMap,5)
# to get the significant SNPs for chromosome X at level 1e-8
getSignificantSNPs(resHapMap,5,sig=1e-8)

}


\keyword{ utilities }
