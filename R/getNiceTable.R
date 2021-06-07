
getNiceTable <- function(x)
  {
    if (!inherits(x, "WGassociation"))
      stop("object 'x' must be of class 'WGassociation'")
    out <- WGstats(x) 
    temp <- lapply(out, getTableSNP)
    tt <- NULL
    nlines <- NULL
    for (i in 1:length(temp))
     {
      tt.i <- temp[[i]]
      nlines <- c(nlines, nrow(temp[[i]]))
      aux <- rbind(c(names(temp)[i], rep(NA,7)), tt.i)
      tt <- rbind(tt, aux)
     }
    colnames(tt)[c(1,7,8)] <- c("SNP", "CI95%", "p-value")
    colnames(tt) <- gsub("%", "\\\\%", colnames(tt))
    tt2 <- gsub("NA","",  tt)
    ans <- gsub("\\(  -  \\)","",  tt2)
    attr(ans, "nlines") <- nlines
    ans
 }
    


getTableSNP <- function(x)
 {
   ff <- function(x)
    {
     ans <- apply(x,1, function(x) paste("(", paste(x, collapse="-"), ")", sep=""))
     ans[1] <- "NA"
     ans
    }

   part1 <- apply(x[,1:5], 2, format)
   part2 <- ff(format(x[,6:7]))
   part3 <- formatC(x[,8])
   ans <- cbind(rownames(x),part1, part2, part3)
   rownames(ans) <- NULL
   ans
 }
   
   


