\name{scanWGassociation}
\alias{scanWGassociation}
\title{Whole genome association analysis}
\description{
  This function also performs a whole genome association analysis as the 
  function \code{\link{WGassociation}} but only computes the p values for
  each SNP corresponding to likelihood ratio test. 
}
\usage{
scanWGassociation(formula, data, model = c("all"), nperm,
   quantitative = is.quantitative(formula, data), genotypingRate = 80)
}
\arguments{
  \item{formula}{either a symbolic description of the model to be fited (a formula object) without the SNP
                or the name of response variable in the case of fitting single models (e.g. unadjusted models). 
                It might have either a continuous variable (quantitative traits) or a 
                factor variable (case-control studies) as the response on the left of the \code{~} 
                operator and terms with additional covariates on the right of the \code{~}. Currently
                only is possible to write \code{~}1 (that is, non-adjusted analysis). See \code{\link{WGassociation}}
                if you are interested in including additional covariates. }
  \item{data}{a required dataframe of class 'setupSNP' containing the variables in the model and the SNPs}
  \item{model}{a character string specifying the type of genetic model (mode of inheritance) for the SNP. 
               This indicates how the genotypes should be collapsed. Possible
               values are "codominant", "dominant", "recessive", "overdominant", "log-additive" or "all". The default
               is "all" that fits the 5 possible genetic models. Only the first words are required, e.g "co", "do", etc.}
  \item{nperm}{number of permutations to simulate the null hypotesis (e.g OR=1), conditioning on the empirical correlation
               structure. Only required to perform a permutation test. Currently this test is only available for binary traits.}
  \item{quantitative}{logical value indicating whether the phenotype (that which is in the
                      left of the operator ~ in 'formula' argument) is quantitative. The function 
                      'is.quantitative' returns FALSE when the phenotype
                      is a variable with two categories (i.e. indicating case-control status). Thus,
                      it is not a required argument but it may be modified by the user.}
  \item{genotypingRate}{ minimum percentage of genotype rate for a given SNP to be included in the analysis.
                  Default is 80\%.}
}
\value{
  An object of class 'WGassociation'. 

  The function 'print' is used to print the results. The p values are saved in the attribute 'pvalues' as a matrix. 
  They may be obtained using attr(,"pvalues") (see examples). The first column indicates whether a problem with 
  genotyping is present.

  The function 'plot' is used to obtain a plot of p values in the -log scale. See plot.WGassociation for
  further details
}

\references{
JR Gonzalez, L Armengol, X Sole, E Guino, JM Mercader, X Estivill, V Moreno.
SNPassoc: an R package to perform whole genome association studies. Bioinformatics,
2007;23(5):654-5.
}

\seealso{ \code{\link{WGassociation}} \code{\link{getSignificantSNPs}} 
    \code{\link{association}} \code{\link{setupSNP}} \code{\link{plot.WGassociation}} 
    \code{\link{permTest}}}

\examples{

# Next steps may be very time consuming. So they are not executed

#data(HapMap)
#myDat<-setupSNP(HapMap, colSNPs=3:9307, sort = TRUE,
#   info=HapMap.SNPs.pos, sep="")
#resHapMap<-scanWGassociation(group~1, data=myDat, model="log")

}
\keyword{utilities}
