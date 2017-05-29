\name{interactionPval}
\alias{interactionPval}
\alias{print.SNPinteraction}
\alias{plot.SNPinteraction}

\title{ Two-dimensional SNP analysis for association studies}
\description{
  Perform a two-dimensional SNP analysis (interaction) for association studies with
  possible allowance for covariate
}
\usage{
interactionPval(formula, data, quantitative = 
          is.quantitative(formula, data), model = "codominant")
}

\arguments{
  \item{formula}{a formula object. It might have either a continuous variable (quantitative traits) or a 
                factor variable (case-control study) as the response on the left of the \code{~} 
                operator and the terms corresponding to the covariates to be adjusted. A crude analysis
                is performed indicating \code{~}1 }
  \item{data}{ a required object of class 'setupSNP'. }
  \item{quantitative}{ logical value indicating whether the phenotype (those which is in the
                      left of the operator ~ in 'formula' argument) is quantitative. The function 
                      'is.quantitative' returns FALSE when the phenotype
                      is a variable with two categories (i.e. indicating case-control status). Thus,
                      it is not a required argument but it may be modified by the user.}
  \item{model}{a character string specifying the type of genetic model (mode of inheritance). This 
           indicates how the genotypes should be collapsed. Possible
           value are "codominant", "dominant", "recessive", "overdominant" or "log-additive". The default
           is "codominant". Only the first words are required, e.g "co", "do", "re", "ov", "log"}
}
\details{
  The 'interactionPval' function calculates, for each pair of SNPs (i,j), the likelihood underling the null model L0, 
  the likelihood under each of the single-SNP, L(i) and L(j), the likelihood under an additive SNP
  model La(i,j), and the likelihood under a full SNP model (including SNP-SNP interaction), Lf(i,j).
  
  The upper triangle in matrix from this function contains the p values for the interaction (epistasis) 
  log-likelihood ratio test, LRT, LRTij = -2 (log Lf(i,j) - log La(i,j))
  
  The diagonal contains the p values from LRT for the crude effect of each SNP, LRTii = -2 (log L(i) - log L0) 
  
  The lower triangle contains the p values from LRT comparing the two-SNP additive likelihood to the best
  of the single-SNP models, LRTji = -2 (log La(i,j) - log max(L(i),L(j))) 
  
  In all cases the models including the SNPs are adjusted by the covariates indicated in the 'formula'
  argument. This method is used either for quantitative traits and dicotomous variables (case-control studies). 
  
}
\value{
  The 'interactionPval' function returns a matrix of class 'SNPinteraction' containing the p values corresponding to the 
  different likelihood ratio tests above describe. 
  
  Methods defined for 'SNPinteraction' objects are provided for print and plot. The plot method uses 'image' to plot a grid 
  of p values.  The upper triangle contains the interaction (epistasis) p values from LRT. The content in the lower triangle 
  is the p values from the LRT comparing the additive model with the best single model. The diagonal contains the main effects 
  pvalues from LRT. The 'plot.SNPinteraction' function also allows the user to plot the SNPs sorted by genomic position and 
  with the information about chromosomes as in the 'plotMissing' function.
  
}

\references{
JR Gonzalez, L Armengol, X Sole, E Guino, JM Mercader, X Estivill, V Moreno.
SNPassoc: an R package to perform whole genome association studies. Bioinformatics,
2007;23(5):654-5.
}

\note{  two-dimensional SNP analysis on a dense grid can take a great deal of computer
     time and memory.
}

\seealso{\code{\link{setupSNP}}}

\examples{

data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")

ansCod<-interactionPval(log(protein)~sex,datSNP)
print(ansCod)
plot(ansCod)

}


\keyword{utilities}