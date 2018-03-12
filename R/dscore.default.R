
dscore.character <- function(x, ...){
  snpmart <- biomaRt::useMart("ENSEMBL_MART_SNP", 
                              dataset = "hsapiens_snp", ...)
  snpInfo <- biomaRt::getBM(c("refsnp_id", "chr_name", "chrom_start", 
                     "allele", "minor_allele", "minor_allele_freq"),
                   filters = c("snp_filter"),
                   values = x, mart = snpmart)

  if (nrow(snpInfo)!=length(x)) {
    warning("Genetic score distribution has been computed with these SNPs: \n",
            unique(snpInfo$refsnp_id))
    
  }
    
  mafs<-snpInfo$minor_allele_freq
  ans <- dscore(mafs)
  attr(ans, "MAFs") <- snpInfo
  ans
}


dscore.default <- function(x, ...)
 {
   if(any(x>1 | x<0))
    stop("argument 'x' should be a vector containing probabilities")
   v <- c(0:(2*length(x)))
   ans <- dpoisbinom(v, rep(x, each=2))
   names(ans) <- paste(v, "alleles")
   class(ans) <- "dscore"
   ans 
 }

