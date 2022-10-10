`setupSNP` <-
function(data, colSNPs, sort=FALSE, info, sep="/", ...)

{
 if (missing(data)) 
        stop("Required argument data is missing")

 if (is.matrix(data)) 
        data<-as.data.frame(data)

 if (!is.data.frame(data)) 
        stop("Argument data is not a data.frame or matrix")

 if (sort)
   {
     temp <- sortSNPs(data, colSNPs, info)
     pos <- temp$pos 
     info <- temp$dataSorted
     temp <- data[, pos, drop=FALSE]
     #..debug..#   dataSNPs <- mclapply(temp, snp, sep = sep, ...)
     dataSNPs <- lapply(temp, snp, sep = sep, ...)
   }
 else
   { 
       #..debug..#  dataSNPs <- mclapply(data[, colSNPs, drop=FALSE], snp, sep = sep, ...)
       dataSNPs <- lapply(data[, colSNPs, drop=FALSE], snp, sep = sep, ...)
   }

 dataSNPs <- data.frame(dataSNPs) 
 datPhen <- data[ , -colSNPs, drop=FALSE]

 ans<-cbind(datPhen,dataSNPs)

 label.SNPs <- names(dataSNPs)
 class(ans) <- c("setupSNP","data.frame")
 attr(ans,"row.names") <- 1:length(ans[[1]])
 attr(ans,"label.SNPs") <- label.SNPs
 attr(ans,"colSNPs") <- c( (length(ans) - length(label.SNPs)+1) : length(ans))
 if (sort)
   attr(ans,"gen.info") <- info
 ans

}

