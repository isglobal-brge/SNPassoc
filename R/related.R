related <- function(x) {
  ans <- NULL
  while(nrow(x)>0) {
   xx <- plyr::arrange(plyr::count(c(x$ID1, x$ID2)), -freq)
   rm.xx <- xx$x[1]
   x <- subset(x, ID1 != rm.xx & x$ID2 != rm.xx)
   ans <- c(as.character(rm.xx), ans)
   ans
  }
  ans
}
