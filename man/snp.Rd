\name{snp}
\alias{snp}
\alias{is.snp}
\alias{as.snp}
\alias{reorder.snp}
\alias{summary.snp}
\alias{plot.snp}
\alias{dominant.snp}
\alias{codominant.snp}
\alias{recessive.snp}
\alias{additive.snp}
\alias{print.snp}
\alias{[.snp}
\alias{print.summary.snp}

\title{SNP object}
\description{

  \code{snp} creates an snp object

  \code{is} returns \code{TRUE} if \code{x} is of class 'snp'

  \code{as} attempts to coerce its argument into an object of class 'snp'

  \code{reorder} change the reference genotype

  \code{summary} gives a summary for an object of class 'snp' including genotype
                     and allele frequencies and an exact thest of Hardy-Weinberg 
                     equilibrium 

  \code{plot}    gives a summary for an object of class 'snp' including genotype
                     and allele frequencies and an exact thest of Hardy-Weinberg 
                     equilibrium in a plot. Barplot or pie are allowed

  \code{[.snp}   is a copy of [.factor modified to preserve all attributes

}
\usage{

  snp(x, sep = "/", name.genotypes, reorder="common", 
    remove.spaces = TRUE, allow.partial.missing = FALSE) 

  is.snp(x)
 
  as.snp(x, ...)

  \method{additive}{snp}(o)

}

\arguments{
  \item{x}{either an object of class 'snp' or an object to be converted to class 'snp'}
  \item{sep}{character separator used to divide alleles when \code{x} is a vector of strings where 
             each string holds both alleles. The default is "/". See below for details.}         
  \item{name.genotypes}{the codes for the genotypes. This argument may be useful when genotypes 
         are coded using three different codes (e.g., 0,1,2 or hom1, het, hom2)}
  \item{reorder}{how should genotypes within an individual be reordered. Possible values are
   'common' or 'minor'. The default is 
    \code{reorder="common"}. In that case, alleles are sorted within each individual by
    common homozygous.}
  \item{remove.spaces}{logical indicating whether spaces and tabs will be
          removed from the genotypes before processing}
  \item{allow.partial.missing}{logical indicating whether one allele is
          permitted to be missing.  When set to 'FALSE' both alleles
          are set to 'NA' when either is missing.}
  \item{o}{an object of class 'snp' to be coded as a linear covariate: 0,1,2}
  \item{...}{optional arguments}
}
\details{
  SNP objects hold information on which gene or marker alleles were
  observed for different individuals.  For each individual, two alleles
  are recorded.

  The snp class considers the stored alleles to be unordered , i.e., "C/T"
  is equivalent to "T/C". It assumes that the order of the alleles is not important.

  When \code{snp} is called, \code{x} is a character vector, and it is
    assumed that each element encodes both alleles. In this case, if
    \code{sep} is a character string, \code{x} is assumed to be coded
    as "Allele1<sep>Allele2".  If \code{sep} is a numeric value, it is
    assumed that character locations \code{1:sep} contain allele 1 and
    that remaining locations contain allele 2.

  \code{additive.snp} recodes the SNPs for being analyzed as a linear covariate (codes 0,1,2) 
  

}

\references{
JR Gonzalez, L Armengol, X Sole, E Guino, JM Mercader, X Estivill, V Moreno.
SNPassoc: an R package to perform whole genome association studies. Bioinformatics,
2007;23(5):654-5.
}


\value{
 
  The snp class extends "factor" where the levels is a character vector of possible 
   genotype values stored coded by \code{paste( allele1, "", allele2, sep="/")}
}
  


\seealso{
  \code{\link{association}}  
  }

\examples{
# some examples of snp data in different formats

dat1  <- c("21", "21", "11", "22", "21",
                    "22", "22", "11", "11", NA)
ans1  <- snp(dat1,sep="")
ans1

dat2 <- c("A/A","A/G","G/G","A/G","G/G",
                    "A/A","A/A","G/G",NA)
ans2  <- snp(dat2,sep="/")
ans2

dat3 <- c("C-C","C-T","C-C","T-T","C-C",
                    "C-C","C-C","C-C","T-T",NA)
ans3 <- snp(dat3,sep="-")
ans3


dat4 <- c("het","het","het","hom1","hom2",
                    "het","het","hom1","hom1",NA)
ans4 <- snp(dat4,name.genotypes=c("hom1","het","hom2"))
ans4


# summary 
summary(ans3)

# plots

plot(ans3)
plot(ans3,type=pie)
plot(ans3,type=pie,label="SNP 10045")

}
\keyword{utilities}

