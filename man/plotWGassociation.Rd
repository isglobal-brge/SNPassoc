\name{plot.WGassociation}
\alias{plot.WGassociation}

\title{Function to plot -log p values from an object of class 'WGassociation'}
\description{
  Function to plot -log p values from an object of class 'WGassociation'
}

\usage{
\method{plot}{WGassociation}(x, ...)
}

\arguments{
  \item{x}{ an object of class 'WGassociation'}
  \item{\dots}{other graphical parameters }
}
\details{
  A panel with different plots (one for each mode of inheritance) are plotted. Each of them represents 
  the -log(p value) for each SNP. Two horizontal lines are also plotted. One one them indicates the nominal 
  statistical significance level whereas the other one indicates the statistical 
  significance level after Bonferroni correction.
}

\references{
JR Gonzalez, L Armengol, X Sole, E Guino, JM Mercader, X Estivill, V Moreno.
SNPassoc: an R package to perform whole genome association studies. Bioinformatics,
2007;23(5):654-5.
}


\value{
 No return value, just the plot
}

\examples{

library(SNPassoc)

data(asthma, package = "SNPassoc")
asthma.s <- setupSNP(data=asthma, colSNPs=7:ncol(asthma), sep="")

ans <- WGassociation(casecontrol, data=asthma.s)

plot(ans)

}

\seealso{ \code{\link{association}} \code{\link{setupSNP}} \code{\link{WGassociation}} }



\keyword{utilities}
