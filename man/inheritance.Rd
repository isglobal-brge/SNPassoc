\name{inheritance}
\alias{inheritance}
\alias{geneticModel}
\alias{codominant}
\alias{dominant}
\alias{recessive}
\alias{overdominant}
\alias{additive}

\title{Collapsing (or recoding) genotypes into different categories (generally two) depending on a given genetic mode of inheritance}
\description{codominant function recodifies a variable having genotypes depending on the allelic frequency
             in descending order. \cr
             
            dominant, recessive, and overdominant functions collapse the three categories of a given SNP 
            into two categories as follows: Let 'AA', 'Aa', and 'aa' be the three genotypes. After determining
            the most frequent allele (let's suppose that 'A' is the major allele) the functions return a vector
            with to categories as follows. dominant: 'AA' and 'Aa-aa'; recessive: 'AA-Aa' and 'aa'; 
            overdominant: 'AA-aa' vs 'Aa'. \cr
            
            additive function creates a numerical variable, 1, 2, 3 corresponding to the three genotypes sorted
            out by descending allelic frequency (this model is referred as log-additive).}
\usage{
codominant(o)
dominant(o)
recessive(o)
overdominant(o)
additive(o)
}

\arguments{
  \item{o}{categorical covariate having genotypes}
}
\value{A snp object collapsing genotypes into different categories depending on a given genetic mode of inheritance}
\examples{
data(SNPs)
dominant(snp(SNPs$snp10001,sep=""))
overdominant(snp(SNPs$snp10001,sep=""))
}
\keyword{utilities}