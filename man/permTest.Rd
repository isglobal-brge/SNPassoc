\name{permTest}
\alias{permTest}
\alias{print.permTest}
\alias{plot.permTest}

\title{Permutation test analysis}

\description{
  This function extract the p values for permutation approach performed using scanWGassociation function
}
\usage{
   permTest(x, method="minimum", K)
}

\arguments{
  \item{x}{a required object of class 'WGassociation' with the attribute 'permTest'. See details}
  \item{method}{statistic used in the permutation test. The default is 'minimum' but 'rtp' (rank truncated product) is
    also available.}
  \item{K}{number of the K most significant p values from the total number of test performed (e.g number of SNPs) used
            to compute the rank truncated product. This argument is only required when method='rtp'. See references }

}

\details{
 This function extract the p values from an object of class 'WGassociation'. This object migth be obtained using 
 the funcion called 'scanWGassociation' indicating the number of permutations in the argument 'nperm'. 
}


\value{
  An object of class 'permTest'. 

  'print' returns a summary indicating the number of SNPs analyzed, the number of valid SNPs (those non-Monomorphic and 
  that pass the calling rate), the p value after Bonferroni correction, and the p values based on permutation approach. 
  One of them is based on considering the empirical percentil for the minimum p values, and the another one on assuming 
  that the minimum p values follow a beta distribution.
  

  'plot' produces a plot of the empirical distribution for the minimum p values (histogram) and the expected distribution
  assuming a beta distribution. The corrected p value is also showed in the plot. 

  See examples for further illustration about all previous issues. 

}

\references{
Dudbridge F, Gusnanto A and Koeleman BPC. Detecting multiple associations in genome-wide studies.  
Human Genomics, 2006;2:310-317.

Dudbridge F and Koeleman BPC. Efficient computation of significance levels for multiple associations in large studies of 
correlated data, including genomewide association studies. Am J Hum Genet, 2004;75:424-435.

JR Gonzalez, L Armengol, X Sole, E Guino, JM Mercader, X Estivill, V Moreno.
SNPassoc: an R package to perform whole genome association studies. Bioinformatics,
2007;23(5):654-5.
}

\seealso{ \code{\link{scanWGassociation}}}

\examples{

library(SNPassoc)

data(asthma, package = "SNPassoc")
asthma.s <- setupSNP(data=asthma, colSNPs=7:ncol(asthma), sep="")

ans <- WGassociation(casecontrol, data=asthma.s)

}

\keyword{utilities}
