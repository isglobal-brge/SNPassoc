\name{odds}
\alias{odds}

\title{Extract odds ratios, 95\% CI and pvalues}
\description{
  Extract odds ratios, 95% CI and pvalues from an WGassociation object fitted with WGassociation function. scanWGassociation doesn't fit models so doesn't estimate OR.

}
\usage{
 odds(x, model=c("log-additive", "dominant", "recessive", "overdominant", "codominant"), 
      sorted=c("no","p-value","or"))
}

\arguments{
  \item{x}{an object of class 'WGassociation' output of WGassociation}
  \item{model}{model to be extracted. Only first one is used. The first letter is enough, low or upper case.}
  \item{sorted}{Sort the output by P value or OR.}
}

\value{
  A matrix with OR 95\% CI (lower, upper) and P value for the selected model. For codominant model, the OR and 95\%CI are given for heterozygous and homozigous.
}

\references{
JR Gonzalez, L Armengol, X Sole, E Guino, JM Mercader, X Estivill, V Moreno.
SNPassoc: an R package to perform whole genome association studies. Bioinformatics,
2007;23(5):654-5.
}


\examples{
 data(SNPs)
 datSNP<-setupSNP(SNPs,6:40,sep="")
 ans<-WGassociation(casco~1,data=datSNP,model="all")
 odds(ans)
}


\keyword{utilities}

