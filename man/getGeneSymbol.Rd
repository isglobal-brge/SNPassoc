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
\examples{
    \donttest{
        snps = c('rs58108140','rs189107123','rs180734498','rs144762171')
        chr = c('chr1','chr1','chr1','chr1')
        pos = c(10583, 10611, 13302, 13327)
        
        x <- data.frame(snps, chr, pos )
        
        getGeneSymbol(x)
    }
}
