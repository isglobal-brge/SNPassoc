% --- SNPstat-internal.Rd ---

\name{SNPassoc-internal}
\alias{SNPassoc-internal}
\alias{association.fit}
\alias{extractPval}
\alias{extractPval.i}
\alias{SNPHWE}
\alias{GenotypeRate}
\alias{intervals.dif}
\alias{intervals.or}
\alias{is.quantitative}
\alias{haplo.inter.fit}
\alias{table.corner}
\alias{table.interaction}
\alias{crea.lab}
\alias{togeno}
\alias{orderChromosome}
\alias{expandsetupSNP}
\alias{pvalTest}
\alias{modelTest}
\alias{assoc}
\alias{codominant.default}
\alias{dominant.default}
\alias{recessive.default}
\alias{overdominant.default}
\alias{additive.default}
\alias{trim}
\alias{interleave}


\title{Internal SNPstat functions}
\description{Internal SNPassoc functions}
\value{No return value, internal calls }
\usage{
association.fit(var, dep, adj, quantitative, type, level, 
       nIndiv, genotypingRate = 0, ...) 
extractPval(x)
extractPval.i(i,x,pos,models)
SNPHWE(x)
GenotypeRate(x)
haplo.inter.fit(geno, var2, dep, adj = NULL, fam, 
          haplo.freq.min, ...)
crea.lab(x,pos.ini,cex,dist)
orderChromosome(x)
togeno(f,sep=sep,lab=lab)
expandsetupSNP(o)
pvalTest(dataX,Y,quantitative,type,genotypingRate)
modelTest(X,Y,quantitative,type,genotypingRate)
assoc(y,x,test="lrt",quantitative)
trim(s)
interleave(..., append.source=TRUE, sep=": ", drop=FALSE)
\method{codominant}{default}(o)
\method{dominant}{default}(o)
\method{recessive}{default}(o)
\method{overdominant}{default}(o)
\method{additive}{default}(o)
}

\details{These are not to be called by the user}
\keyword{internal}
