\name{Table.mean.se}
\alias{Table.mean.se}

\title{Descriptive sample size, mean, and standard error}
\description{
  This function computes sample size, mean and standard error of a quantitative trait for
  each genotype (or combination of genotypes)  
}

\usage{
Table.mean.se(var, dep, subset = !is.na(var))
}

\arguments{
  \item{var}{quantitative trait}
  \item{dep}{variable with genotypes or any combination of them}
  \item{subset}{an optional vector specifying a subset of observations to be used in the descriptive analysis}
}

\value{
  \item{tp}{A matrix giving sample size (n), median (me) and standard error (se) for each genotype}
}

\seealso{  \code{\link{Table.N.Per}} }
\examples{
data(SNPs)
# sample size, mean age and standard error for each genotype
Table.mean.se(SNPs$snp10001,SNPs$protein)

# The same table for a subset (males)
Table.mean.se(SNPs$snp10001,SNPs$protein,SNPs$sex=="Male")

# The same table assuming a dominant model
Table.mean.se(dominant(snp(SNPs$snp10001,sep="")),SNPs$protein,SNPs$sex=="Male")


}
\keyword{utilities}

