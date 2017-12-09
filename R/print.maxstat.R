print.maxstat<-function(x,...)
{
 printCoefmat(t(x), dig.tst=3, tst.ind=1:4, P.values = TRUE, na.print="-", ...)
}

head.maxstat <- function(x, n=5, ...) {
  printCoefmat(t(x[,1:n]), dig.tst=3, tst.ind=1:4, P.values = TRUE, na.print="-", ...)
}