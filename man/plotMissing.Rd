\name{plotMissing}
\alias{plotMissing}

\title{ Plot of missing genotypes}

\description{
 Plot a grid showing which genotypes are missing
}
\usage{
plotMissing(x, print.labels.SNPs = TRUE, 
        main = "Genotype missing data", ...) 
}

\arguments{
  \item{x}{an object of class 'setupSNP'}
  \item{print.labels.SNPs}{should labels of SNPs be printed?}
  \item{main}{title to place on plot}
  \item{...}{extra arguments of 'image' function}
}

\details{
    This function uses 'image' function to plot a grid with black pixels where the genotypes
     are missing.  
}

\seealso{ \code{\link{setupSNP}}}

\value{
  No return value, just the plot
}

\examples{
 data(SNPs)
 data(SNPs.info.pos) 
 ans<-setupSNP(SNPs,colSNPs=6:40,sep="")
 plotMissing(ans)
 
 # The same plot with the SNPs sorted by genomic position and 
 # showing the information about chromosomes

 ans<-setupSNP(SNPs,colSNPs=6:40,sort=TRUE,SNPs.info.pos,sep="") 
 plotMissing(ans)
}
\keyword{utilities}