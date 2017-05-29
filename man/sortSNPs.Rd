\name{sortSNPs}
\alias{sortSNPs}

\title{Sort a vector of SNPs by genomic position}
\description{
  This function sorts a vector with the position of SNPs in a data frame
  using another data frame which contains information about SNPs, their 
  chromosome, and their genomic position
}
\usage{
sortSNPs(data, colSNPs, info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a required data frame with the SNPs}
  \item{colSNPs}{a required vector indicating which columns of 'data' contains genotype data}
  \item{info}{a required data frame with genomic information for the SNPs (chromosome and position).
       The first column must have the SNPs, the second one the chromosome and the third one the genomic 
       position. }
}
\details{
  First of all, the function obtains a vector with the SNPs sorted using the data frame with the genomic 
  positions (see 'dataSNPs.pos' argument).
  Then, the columns which indicate where the information about the 
  genotypes is in our data frame, are sorted using this vector. 
  
  This information is useful when \code{\link{WGassociation}} function is called since it 
  allow the user to summaryze the results with the SNPs sorted by genomic position
}

\value{
  a vector indicating the colums where the SNPs are recorded in our data frame ('data' argument), sorted
  using the genomic positions listed in 'dataSNPs.pos' argument)  
}

\examples{
#
# data(SNPs)
# data(SNPs.info.pos)
# colSNPs.order<-sortSNPs(SNPs,c(6:40),SNPs.info.pos)
#
}
\keyword{utilities}

