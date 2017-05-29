\name{maxstat}
\alias{maxstat}
\alias{maxstat.default}
\alias{maxstat.table}
\alias{maxstat.setupSNP}
\alias{maxstat.matrix}
\alias{print.maxstat}

\title{max-statistic for a 2x3 table}
\description{
  Computes the asymptotic p-value for max-statistic for a 2x3 table
}
\usage{
maxstat(x, \dots)

\method{maxstat}{default}(x, y, \dots)

\method{maxstat}{table}(x, \dots)

\method{maxstat}{setupSNP}(x, y, colSNPs=attr(x,"colSNPs"), ...)

\method{maxstat}{matrix}(x, \dots)

}
\arguments{
  \item{x}{a numeric matrix with 2 rows (cases/controls) and 3 colums (genotypes) or a vector with case/control status
           or an object of class 'setupSNP'.}
  \item{y}{an optional numeric vector containing the information for a given SNP. In this
           case 'x' argument must contain a vector indicarting case/control status. If 'x' argument is an object of 
           class 'setupSNP' this argument migth be the name of the variable containing case/control information.}
  \item{colSNPs}{a vector indicating which columns contain those SNPs to compute max-statistic. By default
                 max-statistic is computed for those SNPs specified when the object of class 'setupSNP' was
                 created.}

  \item{\dots}{further arguments to be passed to or from methods.}
}
\references{
Gonzalez JR, Carrasco JL, Dudbridge F, Armengol L, Estivill X, Moreno V. Maximizing association statistics over 
genetic models (2007). Submitted   

Sladek R, Rocheleau G, Rung J et al. A genome-wide association study identifies novel risk loci for type 2 diabetes (2007). Nature 445, 881-885 

}
\value{
  A matrix with the chi-square statistic for dominant, recessive, log-additive and max-statistic and 
  its asymptotic p-value. 
}
\seealso{
  \code{\link{setupSNP}}
}
\examples{

# example from Sladek et al. (2007) for the SNP rs1111875 
 tt<-matrix(c(77,298,310,122,316,231),nrow=2,ncol=3,byrow=TRUE)
 maxstat(tt)

 data(SNPs)
 maxstat(SNPs$casco,SNPs$snp10001) 
 myDat<-setupSNP(SNPs,6:40,sep="")
 maxstat(myDat,casco)

 

}
\keyword{utilities}