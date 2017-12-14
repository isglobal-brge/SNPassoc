getLocusInfo <- function(x, annotation, file, candidate, size=1e5, ...){
 if(!inherits(x, "GlmTests"))
   stop("x must be an object of class 'GlmTests")
 if (missing(file))
   stop("Indicate the name and path where output file should be saved")

ps <- p.value(x)
SNP <- names(x)
out <- data.frame(SNP=SNP, pvalue=ps)
out.o <- out[order(out$pvalue),]
if (missing(candidate))
 candidate <- as.character(out.o$SNP[1])

chr <- annotation[candidate, "chromosome"]
pos <- annotation[candidate, "position"] 
mask <- annotation$chromosome == chr &
  annotation$position > pos - size &
  annotation$position < pos + size   
snps.sel <- annotation[mask, "snp.name"]

info <- out[out$SNP%in%snps.sel, 1:2]
names(info) <- c("MarkerName", "P.value")
write.table(info, file=file, sep="\t",
            row.names=FALSE, quote=FALSE)
}
