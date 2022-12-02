`summary.setupSNP` <-
function(object, ...)
{
    if (!inherits(object, "setupSNP")) 
        stop("object must be an object of class 'setupSNP'")
    
    colSNPs<-attr(object,"colSNPs")
    
    if(length(colSNPs)>0) {
        temp <- mclapply(object[,colSNPs, drop=FALSE], expandsetupSNP)
    
        ans <- do.call(rbind,temp)
        out<-as.matrix(ans)
        dimnames(out)[[2]][4] <- "missing (%)"
        print(out, quote=FALSE, na.print="-")
        
    } else {
        class(object)<-"data.frame"
        ans<-summary(object)
        print(ans)
    }
    invisible(ans)
}

