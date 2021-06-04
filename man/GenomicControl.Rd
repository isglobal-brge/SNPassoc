\name{GenomicControl}
\alias{GenomicControl}

\title{Population substructure}

\description{
  This function estimates an inflation (or deflation) factor, lambda, as indicated in the paper
  by Devlin et al. (2001) and corrects the p-values using this factor. 
}

\usage{
GenomicControl(x, snp.sel)
}

\arguments{
  \item{x}{an object of class 'WGassociation'. }
  \item{snp.sel}{SNPs used to compute lambda. Not required.}
}

\details{
 This method is only valid for 2x2 tables. This means that the object of class 'WGassociation' might not
 have fitted the codominant model.

 See reference for further details.
}

\value{
  The same object of class 'WGassociation' where the p-values have been corrected for genomic control.
}

\references{B Devlin, K Roeder, and S.A. Bacanu. Unbiased Methods for Population Based Association Studies. 
            Genetic Epidemiology (2001) 21:273-84
}

\seealso{\code{\link{qqpval}}, \code{\link{WGassociation}}}

\examples{

data(SNPs) 
datSNP<-setupSNP(SNPs,6:40,sep="")
res<-WGassociation(casco,datSNP,model=c("do","re","log-add"))

# Genomic Control 
resCorrected<-GenomicControl(res)
}



\keyword{utilities}
