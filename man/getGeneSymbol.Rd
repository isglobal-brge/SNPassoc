\name{getGeneSymbol}
\alias{getGeneSymbol}
\title{Get gene symbol from a list of SNPs}
\usage{
getGeneSymbol(
  x,
  snpCol = 1,
  chrCol = 2,
  posCol = 3,
  db = TxDb.Hsapiens.UCSC.hg19.knownGene
)
}
\arguments{
\item{x}{data.frame containing: SNP name, chromosome and genomic position.}

\item{snpCol}{column of x having the SNP name. Default is 1.}

\item{chrCol}{column of x having the SNP chromosome. Default is 2.}

\item{posCol}{column of x having the SNP position. Default is 3.}

\item{db}{reference genome. Default is 'TxDb.Hsapiens.UCSC.hg19.knownGene'}
}
\value{
a data.frame having initial information and gene symbol
}
\description{
Get gene symbol from a list of SNPs
}
