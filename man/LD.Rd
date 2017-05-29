\name{LD}
\alias{LD}
\alias{LD.snp}
\alias{LD.setupSNP}
\alias{LDplot}
\alias{LDtable}

\title{max-statistic for a 2x3 table}
\description{
  Compute pairwise linkage disequilibrium between genetic markers
}
\usage{
LD(g1, \dots)

\method{LD}{snp}(g1, g2, \dots)

\method{LD}{setupSNP}(g1, SNPs, \dots)


LDplot(x, digits = 3, marker, distance, which = c("D", "D'",
    "r", "X^2", "P-value", "n", " "), ...)

LDtable(x, colorcut = c(0, 0.01, 0.025, 0.05, 0.1, 1), 
       colors = heat.colors(length(colorcut)),
       textcol = "black", digits = 3, show.all = FALSE, 
       which = c("D", "D'", "r", "X^2", "P-value", "n"), 
       colorize = "P-value", cex, \dots)


}
\arguments{

  \item{g1}{genotype object or dataframe containing genotype objects}
  \item{g2}{genotype object (ignored if g1 is a dataframe)}  
  \item{SNPs}{columns containing SNPs}  
  \item{x}{LD or LD.data.frame object} 
  \item{digits}{Number of significant digits to display}
   \item{which}{Name(s) of LD information items to be displayed}
\item{colorcut}{P-value cutoffs points for colorizing LDtable}
  \item{colors}{Colors for each P-value cutoff given in 'colorcut' for LDtable}
 \item{textcol}{Color for text labels for LDtable}
  \item{marker}{Marker used as 'comparator' on LDplot.  If omitted separate
          lines for each marker will be displayed}
\item{distance}{Marker location, used for locating of markers on LDplot.}
\item{show.all}{If TRUE, show all rows/columns of matrix. Otherwise omit
          completely blank rows/columns.}
\item{colorize}{LD parameter used for determining table cell colors}
     \item{cex}{Scaling factor for table text. If absent, text will be scaled
          to fit within the table cells.}

     \item{\dots}{Optional arguments ('plot.LD.data.frame' passes these to
          'LDtable' and 'LDplot').}
}

\references{
 genetics R package by Gregory Warnes et al.  (warnes@bst.rochester.edu)
}

\author{
 functions adapted from LD, LDtable and LDplot in package genetics by Gregory Warnes et al.  (warnes@bst.rochester.edu)
}

\value{
  None
}

\seealso{
  \code{\link{setupSNP}}   \code{\link{snp}}
}


\keyword{utilities}
