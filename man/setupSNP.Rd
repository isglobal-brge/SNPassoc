\name{setupSNP}
\alias{setupSNP}
\alias{summary.setupSNP}
\alias{plot.setupSNP}
\alias{[.setupSNP}
\alias{[[<-.setupSNP}
\alias{[<-.setupSNP}
\alias{$<-.setupSNP} 
\alias{labels.setupSNP}

\title{Convert columns in a dataframe to class 'snp'}

\description{
  \code{setupSNP} Convert columns in a dataframe to class 'snp'

  \code{summary.setupSNP} gives a summary for an object of class 'setupSNP' including 
                     allele names, major allele frequencie, an exact thest of Hardy-Weinberg 
                     equilibrium and percentage of missing genotypes

}

\usage{
setupSNP(data, colSNPs, sort = FALSE, info, sep = "/", ...)
}

\arguments{
  \item{data}{dataframe containing columns with the SNPs to be converted}
  \item{colSNPs}{Vector specifying which columns contain SNPs data}
  \item{sort}{should SNPs be sorted. Default is FALSE}
  \item{info}{if sort is TRUE a dataframe containing information 
       about the SNPs regarding their genomic position and the gene where 
       they are located}
  \item{sep}{character separator used to divide alleles in the genotypes}
  \item{...}{optional arguments}
}
\value{
 a dataframe of class 'setupSNP' containing converted SNP variables. 
 All other variables will be unchanged.
}

\references{
JR Gonzalez, L Armengol, X Sole, E Guino, JM Mercader, X Estivill, V Moreno.
SNPassoc: an R package to perform whole genome association studies. Bioinformatics,
2007;23(5):654-5.
}


\seealso{ \code{\link{snp}}}

\examples{

 data(SNPs)
 myDat<-setupSNP(SNPs,6:40,sep="")


#sorted SNPs and having genomic information
 data(SNPs.info.pos)
 myDat.o<-setupSNP(SNPs,6:40,sep="",sort=TRUE, info=SNPs.info.pos)

# summary
 summary(myDat.o)

# plot one SNP
  plot(myDat,which=2)

}
\keyword{utilities}
