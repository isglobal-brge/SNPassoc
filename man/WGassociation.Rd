\name{WGassociation}
\alias{WGassociation}
\alias{WGstats}
\alias{print.WGassociation}
\alias{summary.WGassociation}
\alias{labels.WGassociation}
\alias{pvalues}
\alias{pvalues.WGassociation}
\alias{codominant.WGassociation}
\alias{dominant.WGassociation}
\alias{recessive.WGassociation}
\alias{overdominant.WGassociation}
\alias{additive.WGassociation}
\alias{[.WGassociation}

\title{Whole genome association analysis}

\description{
  This function carries out a whole genome association analysis between the SNPs and
  a dependent variable (phenotype) under five different genetic models (inheritance patterns): 
  codominant, dominant, recessive, overdominant and log-additive. The phenotype may be quantitative 
  or categorical. In the second case (e.g. case-control studies) this variable must be of class 'factor' 
  with two levels.

}
\usage{
   WGassociation(formula, data, model = c("all"), 
                 quantitative = is.quantitative(formula, data),
                 genotypingRate = 80, level = 0.95, ...)
}

\arguments{
  \item{formula}{either a symbolic description of the model to be fited (a formula object) without the SNP
                or the name of response variable in the case of fitting single models (e.g. unadjusted models).
                It might have either a continuous variable (quantitative traits) or a 
                factor variable (case-control studies) as the response on the left of the \code{~} 
                operator and terms with additional covariates on the right of the \code{~} operator may be 
                added to fit an adjusted model (e.g., \code{~}var1+var2+...+varN+SNP). See details}
  \item{data}{a required dataframe of class 'setupSNP' containing the variables in the model and the SNPs}
  \item{model}{a character string specifying the type of genetic model (mode of inheritance) for the SNP. 
               This indicates how the genotypes should be collapsed. Possible
               values are "codominant", "dominant", "recessive", "overdominant", "log-additive" or "all". The default
               is "all" that fits the 5 possible genetic models. Only the first words are required, e.g "co", "do", etc.}
  \item{quantitative}{logical value indicating whether the phenotype (that which is in the
                      left of the operator ~ in 'formula' argument) is quantitative. The function 
                      'is.quantitative' returns FALSE when the phenotype
                      is a variable with two categories (i.e. indicating case-control status). Thus,
                      it is not a required argument but it may be modified by the user.}
  \item{genotypingRate}{ minimum percentage of genotype rate for a given SNP to be included in the analysis.
                  Default is 80\%.}
  \item{level}{signification level for confidence intervals. Defaul 95\%.}
  \item{...}{Other arguments to be passed through glm function}
}

\details{
  This function assesses the association between the response variable included in the left side in 
  the 'formula' and the SNPs included in the 'data' argument adjusted by those variables included 
  in the right side of the 'formula'. Different genetic models may be analyzed using 'model' argument.
}


\value{
  An object of class 'WGassociation'. 

  'summary' returns a summary table by groups defined in info (genes/chromosomes).

  'WGstats' returns a detailed output, similar to the produced by \code{\link{association}}.

  'pvalues' and 'print' return a table of p-values for each genetic model for each SNP.
  The first column indicates whether a problem with genotyping is present.

  'plot' produces a plot of p values in the -log scale. See \code{\link{plot.WGassociation}} for
  further details.

  'labels' returns the names of the SNPs analyzed.

  The functions 'codominat', 'dominant', 'recessive', 'overdominant' and 'additive'
  are used to obtain the p values under these genetic models.

  See examples for further illustration about all previous issues. 

}

\references{
JR Gonzalez, L Armengol, X Sole, E Guino, JM Mercader, X Estivill, V Moreno.
SNPassoc: an R package to perform whole genome association studies. Bioinformatics,
2007;23(5):654-5.
}


\seealso{ \code{\link{getSignificantSNPs}} \code{\link{association}} 
          \code{\link{WGstats}} \code{\link{setupSNP}} 
          \code{\link{plot.WGassociation}}}

\examples{
data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")
ansAll<-WGassociation(protein~1,data=datSNP,model="all")

# In that case the formula is not required. You can also write:
# ansAll<-WGassociation(protein,data=datSNP,model="all")


#only codominant and log-additive
ansCoAd<-WGassociation(protein~1,data=datSNP,model=c("co","log-add"))

#for printing p values
print(ansAll)
print(ansCoAd)

#for obtaining a matrix with the p palues
pvalAll<-pvalues(ansAll)
pvalCoAd<-pvalues(ansCoAd)

# when all models are fitted and we are interested in obtaining 
# p values for different genetic models

# codominant model
pvalCod<-codominant(ansAll)

# recessive model
pvalRec<-recessive(ansAll)

# and the same for additive, dominant or overdominant


#summary
summary(ansAll)

#for a detailed report
WGstats(ansAll)

#for plotting the p values
plot(ansAll)

}

\keyword{utilities}
