\name{Bonferroni.sig}
\alias{Bonferroni.sig}

\title{ Bonferroni correction of p values}
\description{
  This function shows the SNPs that are statistically significant after correcting
  for the number of tests performed (Bonferroni correction) for an object of class
  "WGassociation"
}
\usage{
 Bonferroni.sig(x, model = "codominant", alpha = 0.05, 
      include.all.SNPs=FALSE)
}

\arguments{
  \item{x}{an object of class 'WGassociation'.}
  \item{model}{a character string specifying the type of genetic model (mode of inheritance). This 
           indicantes how the genotypes should be collapsed when 'plot.summary' is TRUE. Possible
           values are "codominant", "dominant", "recessive", "overdominant", or "log-additive". 
           The default is "codominant". Only the first words are required, e.g "co", "do", ... . }
  \item{alpha}{nominal level of significance. Default is 0.05}
  \item{include.all.SNPs}{logical value indicating whether all SNPs are considered in the Bonferroni
    correction. That is, the number of performed tests is equal to the number of SNPs or equal to the
    number of SNPs where a p value may be computed. The default value is FALSE indicating that the 
    number of tests is equal to the number of SNPs that are non Monomorphic and the rate of genotyping
    is greater than the percentage indicated in the \code{GeneticModel.pval} function. }     
}
\details{
  After deciding the genetic model, the function shows the SNPs that are statistically significant at 
  \code{alpha} level corrected by the number of performed tests. 
}
\value{
 A data frame with the SNPs and the p values for those SNPs that are statistically significant
 after Bonferroni correction
}


\seealso{ \code{\link{WGassociation}}}
\examples{
data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")
ans<-WGassociation(protein~1,data=datSNP,model="all")
Bonferroni.sig(ans, model="codominant", alpha=0.05, include.all.SNPs=FALSE)

}
\keyword{utilities}
