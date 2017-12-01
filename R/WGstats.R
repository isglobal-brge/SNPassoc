`WGstats` <-
function (object, ...) 
{
    if (!inherits(object, "WGassociation")) 
        stop("object must be an object of class 'WGassociation'")

    if (!is.null(attr(object,"fast")))
       stop("\n summary is implemented only for 'WGassociation' function")

    ans <- attr(object,"tables")
    mostattributes(ans)<-NULL
    ans
}

print.WGstats <- function(x, ...){
  print(x, na.print = "", ...)
}