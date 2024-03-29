\name{SNPs}
\alias{SNPs}

\title{SNPs in a case-control study}
\usage{data(SNPs)}
\description{
 SNPs data.frame contains selected SNPs and other clinical covariates for cases and controls 
    in a case-control study

 SNPs.info.pos data.frame contains the names of the SNPs included in the data set 'SNPs' including
    their chromosome and their genomic position
}
\format{
  'SNPs' data.frame contains the following columns:
  \tabular{ll}{
	id \tab identifier of each subject \cr
	casco \tab case or control status: 0-control, 1-case \cr
	sex \tab gender: Male and Female \cr
	blood.pre \tab arterial blood presure \cr
	protein \tab protein levels \cr
	snp10001 \tab SNP 1 \cr
	snp10002 \tab SNP 2 \cr
	  ...        \tab  ...      \cr
	snp100036 \tab SNP 36 \cr
  }
}

\source{
 Data obtained from our department. The reference and details will be supplied after being published.
}

\keyword{datasets}
