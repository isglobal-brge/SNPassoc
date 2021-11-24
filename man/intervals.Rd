\name{intervals}
\alias{intervals}
\alias{intervals.haplo.glm}
\alias{print.intervals}
\alias{summary.haplo.glm}

\title{Print ORs and 95\% confidence intervals for an object of class 'haplo.glm'}

\description{Print ORs and confidence intervals for an object of class 'haplo.glm'}

\usage{
intervals(o, level=.95, ...)
}

\arguments{
  \item{o}{object of class 'haplo.glm'}
  \item{level}{significance level. Default is 95 percent}
  \item{...}{other arguments}

}
\value{ intervals object with ORs and 95\% confidence intervals for an object of class 'haplo.glm'}

\examples{
# Not Run
library(SNPassoc)
library(haplo.stats)

data(asthma, package = "SNPassoc")

asthma.s <- setupSNP(data=asthma, colSNPs=7:ncol(asthma), sep="")
trait <- asthma.s$casecontrol
snpsH <- c("rs714588", "rs1023555",  "rs898070")
genoH <- make.geno(asthma.s, snpsH)

mod <- haplo.stats:: haplo.glm( trait ~ genoH,           
                                family="binomial", 
                                locus.label=snpsH,
                                allele.lev=attributes(genoH)$unique.alleles,
                                control = haplo.glm.control(haplo.freq.min=0.05))   
intervals(mod)
summary(mod)

}
\keyword{utilities}