dscore<-function(x, ...)
 UseMethod("dscore")


plot.dscore <- function(x, ...) {
  xx <- gsub(" alleles", "", names(x))
  plot(xx, x, type="h", xlab="risk alleles", 
       ylab="exact probability", col="red", lwd=2)
  
}

pscore <- function(q, x, ...){
  gs <- dscore(x, ...)
  ans <- sum(gs[1:(q+1)])
  ans
}