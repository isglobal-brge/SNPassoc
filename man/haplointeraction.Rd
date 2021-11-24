\name{haplo.interaction}
\alias{haplo.interaction}
\alias{print.haploOut}

\title{Haplotype interaction with a covariate}

\description{
  This function computes the ORs (or mean differences if a quantitative trait is analyzed) and their 95\% confidence intervals corresponding to
  an interaction between the haplotypes and a categorical covariate
}

\usage{
haplo.interaction(formula, data, SNPs.sel, quantitative = 
  is.quantitative(formula, data), haplo.freq.min = 0.05, ...)
}

\arguments{
  \item{formula}{ a symbolic description of the model to be fitted (a formula object). 
                It might have either a continuous variable (quantitative traits) or a 
                factor variable (case-control studies) as the response on the left of the \code{~} 
                operator and a term corresponding to the interaction variable on the right indicated using 
                'interaction' function (e.g. \code{~}int(var), where var is a factor variable) and it is 
                required. Terms with additional covariates on the the right of the ~ operator may be 
                added to fit an adjusted model (e.g., \code{~}var1+var2+...+varN+int(var)). }
  \item{data}{ an object of class 'setupSNP' containing the variables in the model and the SNPs that will be used to estimate the
               haplotypes. }
  \item{SNPs.sel}{a vector indicating the names of SNPs that are used to estimate the haplotypes}
  \item{quantitative}{logical value indicating whether the phenotype (which is on the
                      left of the operator \code{~} in 'formula' argument) is quantitative. The function 
                      'is.quantitative' returns FALSE when the phenotype
                      is a variable with two categories (i.e. indicating case-control status). Thus,
                      it is not a required argument but it may be modified by the user. }
  \item{haplo.freq.min}{control parameter for haplo.glm included in 'haplo.glm.control'. This parameter corresponds to the minimum 
                        haplotype frequency for a haplotype to be included in the regression model as its own effect. The
                        haplotype frequency is based on the EM algorithm that estimates haplotype frequencies independently
                        of any trait. }
  \item{...}{additional parameters for 'haplo.glm.control'.}                      
}

\details{
  The function estimates the haplotypes for the SNPs indicated in the 'SNPs.sel' argument. Then, usign 'haplo.glm' function (from 'haplo.stats'
  library) estimates the interaction between these haplotypes and the covariate indicated in the formula by means of 'interaction' function. 
}

\value{
  Three different tables are given. The first one corresponds to the full interaction matrix where the ORs (or mean differences if a quantitative 
  trait is analyzed)  are expressed with respect to the most frequent haplotype and the first category of the covariate. The other two tables 
  show the ORs (or mean differences if a quantitative trait is analyzed) and their 95\% confidence intervals for both marginal models. 
  P values for interaction are also showed in the output. 
  
}

\examples{
# not Run
library(SNPassoc)
library(haplo.stats)

data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")
res <- haplo.interaction(log(protein)~int(sex), data=datSNP,
                    SNPs.sel=c("snp100019","snp10001","snp100029"))
res
}
\keyword{utilities }