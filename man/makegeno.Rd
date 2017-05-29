\name{make.geno}
\alias{make.geno}

\title{Create a group of locus objects from some SNPs, assign to 'model.matrix' class.
}

\description{
  This function prepares the CRITICAL element corresponding to matrix of genotypes necessary to be included in 'haplo.glm' function.
}

\usage{
make.geno(data, SNPs.sel)
}

\arguments{
  \item{data}{an object of class 'setupSNP' containing the the SNPs that will be used to estimate the haplotypes.}
  \item{SNPs.sel}{ a vector indicating the names of SNPs that are used to estimate the haplotypes  }
}


\value{
  the same as 'setupGeno' function, from 'haplo.stats' library, returns
}


\seealso{\code{\link{snp}}}

\examples{
## Not run:
data(SNPs)
# first, we create an object of class 'setupSNP'
datSNP<-setupSNP(SNPs,6:40,sep="")
geno<-make.geno(datSNP,c("snp10001","snp10002","snp10003"))
## End(Not run)


}


\keyword{utilities}