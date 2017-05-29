\name{plot.WGassociation}
\alias{plot.WGassociation}

\title{Function to plot -log p values from an object of class 'WGassociation'}
\description{
  Function to plot -log p values from an object of class 'WGassociation'
}

\usage{
\method{plot}{WGassociation}(x, alpha = 0.05, plot.all.SNPs = FALSE, 
  print.label.SNPs = TRUE, cutPval = c(0, 1e-10, 1), whole, 
  ylim.sup=ifelse(is.null(attr(x,"fast")),1e-40, 1e-30),  
  col.legend = c("red", "gray60"), sort.chromosome=TRUE, 
  centromere, ...)
}

\arguments{
  \item{x}{ an object of class 'WGassociation'}
  \item{alpha}{ statistical significance nominal level. See details }
  \item{plot.all.SNPs}{ are all SNPs plotted? If not, neither monomorphic nor SNPs with genotyping problems
                        are plotted. The default is FALSE.}
  \item{print.label.SNPs}{ are labels of SNPs printed? The default is TRUE}
  \item{cutPval}{when argument 'whole' is TRUE in the 'x' object (e.g. when whole genome analysis is carried
            out) 'cutPval divides the range of p values into intervals and codes these values
            according to which interval they fall (like 'cut' function). The default is c(0, 1e-10, 1). That 
            is, the p values are divided in those less than 1e-10 and those greater than 1e-10.  
}
  \item{whole}{is a whole genome carried out? If TRUE 'dataSNPs.pos' argument in
           'setupSNP' is required.}

  \item{ylim.sup}{superior limit for each panel. This value helps to obtain nicer plots}
  \item{col.legend}{ the color of the bar corresponding to p values plotted in each panel. The default is 
          "red" for those p values less than 1e-10 and "gray60" for those greater than 1e-10. }
  \item{sort.chromosome}{should chromosome be sorted? Set this argument to FALSE when genomic information
               corresponds to different genes.}
  \item{centromere}{ numeric vector specifying the centromere positions. If missing, the default centromere 
          value of human genome are used.}
  \item{\dots}{other graphical parameters }
}
\details{
  A panel with different plots (one for each mode of inheritance) are plotted. Each of them represents 
  the -log(p value) for each SNP. Two horizontal lines are also plotted. One one them indicates the nominal 
  statistical significance level (see 'alpha' argument) whereas the other one indicates the statistical 
  significance level after Bonferroni correction.
  
  A different plot is created when the argument 'whole' the object 'x' is TRUE (see setupSNP function). 
  In that case a plot of p values in the -log sclae is plotted for each SNP and for each chromosome. 

}

\references{
JR Gonzalez, L Armengol, X Sole, E Guino, JM Mercader, X Estivill, V Moreno.
SNPassoc: an R package to perform whole genome association studies. Bioinformatics,
2007;23(5):654-5.
}


\value{
 No return value, just the plot
}

\seealso{ \code{\link{association}} \code{\link{setupSNP}} \code{\link{WGassociation}} }



\keyword{utilities}
