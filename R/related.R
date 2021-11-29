#' Get related samples
#' 
#' @param x An object obtained from SNPrelate package.
#' @return A matrix with related individuals.
#' @examples
#'   library(SNPassoc)
#'   data(SNPs)
#' 
#' @export
related <- function(x) {
  ans <- NULL
  while(nrow(x)>0) {
   xx <- plyr::count(c(x$ID1, x$ID2))
   o <- order(xx$freq, decreasing = TRUE)
   xx <- xx[o,]
   rm.xx <- xx$x[1]
   x <- subset(x, x$ID1 != rm.xx & x$ID2 != rm.xx)
   ans <- c(as.character(rm.xx), ans)
   ans
  }
  ans
}
