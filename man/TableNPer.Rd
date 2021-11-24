\name{Table.N.Per}
\alias{Table.N.Per}

\title{ Descriptive sample size and percentage}
\description{
  This function computes sample size and percentage for each category of a categorical trait (e.g.
  case-control status) for each genotype (or combination of genotypes).  
}

\usage{
Table.N.Per(var, dep, subset = !is.na(var))
}

\arguments{
  \item{var}{categorical trait.}
  \item{dep}{variable with genotypes or any combination of them}
  \item{subset}{an optional vector specifying a subset of observations to be used in the descriptive analysis. }
}

\value{
  \item{tp}{A matrix giving sample size (n),and the percentage (\%) for each level of the categorical trait
   for each genotype}
}


\seealso{ \code{\link{Table.mean.se}} }
\examples{
data(SNPs)
#sample size and percentage of cases and controls for each genotype 
Table.N.Per(SNPs$snp10001,SNPs$casco)

# The same table for a subset (males)
Table.N.Per(SNPs$snp10001,SNPs$casco,SNPs$sex=="Male")

# The same table assuming a dominant model
Table.N.Per(dominant(snp(SNPs$snp10001,sep="")),SNPs$casco,SNPs$sex=="Male")


}
\keyword{utilities}

