dscore.character <- function(x, ...){
 
  snpmart <- biomaRt::useEnsembl(biomart="snp", 
                       dataset="hsapiens_snp")

  
  snpInfo <- biomaRt::getBM(c("refsnp_id", "chr_name", 
                              "chrom_start", 
                              "allele", 
                              "minor_allele", 
                              "minor_allele_freq"),
                            filters = c("snp_filter"),
                            values = x, mart = snpmart)
    
  rownames(snpInfo) <- snpInfo$refsnp_id
  
  if (nrow(snpInfo)!=length(x)) {
    warning("Genetic score distribution has been computed with these SNPs: \n",
            paste(unique(snpInfo$refsnp_id), collapse="; "))
  }
  x.info <- x[x%in%snpInfo$refsnp_id]    
  snpInfo <- snpInfo[x.info,]
  rownames(snpInfo) <- snpInfo$refsnp_id
  mafs<-snpInfo$minor_allele_freq
  ans <- dscore(mafs)
  attr(ans, "MAFs") <- snpInfo
  ans
}
