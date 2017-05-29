\name{association}
\alias{association}
\alias{print.snpOut}

\title{Association analysis between a single SNP and a given phenotype}

\description{
  This function carries out an association analysis between a single SNP and
  a dependent variable (phenotype) under five different genetic models (inheritance patterns): 
  codominant, dominant, recessive, overdominant and log-additive. The phenotype may be quantitative 
  or categorical. In the second case (e.g. case-control studies) this variable must be of class 'factor' 
  with two levels.
}

\usage{

association(formula, data, model=c("all"), model.interaction=
       c("codominant"), subset, name.snp = NULL, quantitative = 
       is.quantitative(formula,data), genotypingRate= 0, 
       level = 0.95, ...)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fited (a formula object). 
                It might have either a continuous variable (quantitative traits) or a 
                factor variable (case-control studies) as the response on the left of the \code{~} 
                operator and a term corresponding to the SNP on the right. This term must be
                of class \code{snp} (e.g. \code{~}snp(var), where var is a given SNP), and it is 
                required. Terms with additional covariates on the right of the \code{~} operator may be 
                added to fit an adjusted model (e.g., \code{~}var1+var2+...+varN+SNP). The formula 
                allows to incorporate more than one object of class \code{snp}. In that case, 
                the analysis is done for the first SNP which appears in the formula adjusted by the others 
                covariates and other additional SNPs. }
  \item{data}{a required dataframe of class 'setupSNP' containing the variables in the model.}
  \item{model}{a character string specifying the type of genetic model (mode of inheritance) for the SNP. 
               This indicates how the genotypes should be collapsed. Possible
               values are "codominant", "dominant", "recessive", "overdominant", "additive" or "all". The default
               is "all" that fits the 5 possible genetic models. Only the first words are required, e.g "co", "do", etc.}              
  \item{model.interaction}{a character string specifying the type of genetic model (mode of inheritance) assumed for
                           the SNP when it is included in a interaction term. Possible values are "codominant", 
                           "dominant", "recessive", "overdominant". The default is "codominant".}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process}              
  \item{name.snp}{optional label of the SNP variable to be printed.}
  \item{quantitative}{logical value indicating whether the phenotype (that which is in the
                      left of the operator ~ in 'formula' argument) is quantitative. The function 
                      'is.quantitative' returns FALSE when the phenotype
                      is a variable with two categories (i.e. indicating case-control status). Thus,
                      it is not a required argument but it may be modified by the user.}
  \item{genotypingRate}{ minimum percentage of genotype rate for the SNP to be analyzed.
                  Default is 0\% (e.g. all SNPs are analyzed). This parameter should not be
                  changed. It is used in the function 'WGassociation'.}
  \item{level}{signification level for confidence intervals.}
  \item{...}{Other arguments to be passed through glm function}
}

\details{
  This function should be called by the user when we are interested in analyzing an unique SNP.
  It is recommended to use \code{\link{WGassociation}} function when more than one SNP is studied. \cr
    
}

\value{
  For each genetic model (codominant, dominant, recessive, overdominant, and log-additive) the function gives 
  a matrix with sample size and percentages for each genotype, the Odds Ratio and its 95\% confidence interval
  (taking the most frequent homozygous genotype as the reference), the p-value corresponding to the likelihood ratio test obtained 
  from a comparison with the null model, and the Akaike Information Criterion (AIC) of each genetic model. In the case 
  of analyzing a quantitative trait, the function returns a matrix with sample size, mean and standard errors for each genotype, 
  mean difference and its 95\% confidence interval with respect to the most frequent homozygous genotype, 
  the p-value obtained from an overall gene effect and the Akaike Information Criterion (AIC) of each genetic model. 
  
  When an interaction term (a categorical covariate with an SNP) is included in the model, three different tables are given. 
  The first one correponds to the full interaction matrix where the ORs (or mean differences if a quantitative trait is analyzed) 
  are expressed with respect to the non variant genotype and the first category of the covariate. The other two tables show the ORs 
  and their 95\% confidence intervals for both marginal models. P values for interaction and trend are also showed in the output. 
  
}

\references{
JR Gonzalez, L Armengol, X Sole, E Guino, JM Mercader, X Estivill, V Moreno.
SNPassoc: an R package to perform whole genome association studies. Bioinformatics,
2007;23(5):654-5.

Iniesta R, Guino E, Moreno V. Statistical analysis of genetic polymorphisms in epidemiological studies.
Gac Sanit. 2005;19(4):333-41.

Elston RC. Introduction and overview. Statistical methods in genetic epidemiology. Stat Methods Med Res.
2000;9:527-41.
 
}


\seealso{\code{\link{WGassociation}}}

\examples{
data(SNPs)

# first, we create an object of class 'setupSNP'
datSNP<-setupSNP(SNPs,6:40,sep="")

# case-control study, crude analysis
association(casco~snp10001, data=datSNP)

# case-control study, adjusted by sex and arterial blood pressure
association(casco~sex+snp10001+blood.pre, data=datSNP)


# quantitative trait, crude analysis
association(log(protein)~snp10001,data=datSNP)
# quantitative trait, adjusted by sex
association(log(protein)~snp10001+sex,data=datSNP)


#
# Interaction analysis
#

# Interaction SNP and factor
association(log(protein)~snp10001*sex+blood.pre, data=datSNP, 
          model="codominant")

# Interaction SNP and SNP (codominant and codominant)
association(log(protein)~snp10001*factor(snp10002)+blood.pre, 
          data=datSNP, model="codominant")

# Interaction SNP and SNP (dominant and recessive)
association(log(protein)~snp10001*factor(recessive(snp100019))+blood.pre, 
          data=datSNP, model="dominant")


}

\keyword{utilities}

