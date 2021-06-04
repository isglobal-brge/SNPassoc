dscore.default <- function(x, ...)
 {
   if(all(is.numeric(x)) & any(x>1 | x<0))
    stop("argument 'x' should be a vector containing probabilities (e.g. values in [0,1]")
   v <- c(0:(2*length(x)))
   ans <- dpoisbinom(v, rep(x, each=2))
   names(ans) <- paste(v, "alleles")
   class(ans) <- "dscore"
   ans 
 }

