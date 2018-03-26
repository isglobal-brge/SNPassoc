`scanWGassociation` <-
function (formula, data, model = c("all"), nperm, quantitative = is.quantitative(formula, 
    data), genotypingRate = 80) 
{
    if (!inherits(data, "setupSNP")) 
        stop("data must be an object of class 'setupSNP'")
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m0)]
    if (length(grep("~", mf[[2]])) == 0) {
        formula <- as.formula(paste(mf[[2]], "~1", sep = ""))
        formula.1 <- list(formula)
        mode(formula.1) <- "call"
        mf[2] <- formula.1
    }
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    varDep<-mf[,1]
    control.missing<-dimnames(mf)[[1]]
    mt <- attr(mf, "terms")
    temp0 <- as.character(mt)

    ptm<-proc.time()

# Por ahora no
#    adj <- paste(temp0[2], temp0[1], temp0[3])
    if (temp0[3]!=1)
      stop("adjusted analysis is not yet implemented. Try 'WGassociation'")

    cat("Be patient. The program is computing ... \n")

    Terms <- if (missing(data)) 
        terms(formula)
    else terms(formula, data = data)
    ord <- attr(Terms, "order")
    if (any(ord > 1)) 
        stop("interaction term is not implemented")

    colSNPs <- attr(data, "colSNPs")
    if (is.vector(colSNPs) & length(colSNPs) > 0) 
        dataSNPs <- data[control.missing, colSNPs, drop=FALSE]
    else stop("data should have an attribute called 'colSNPs'. Try again 'setupSNP' function")

    type <- charmatch(model, c("codominant", "dominant", "recessive", 
        "overdominant", "log-additive", "all"))
    type <- sort(type)
    if (any(type %in% 6)) 
        type <- 1:5
    if (length(type)==0) 
        stop("model must be 'codominant','dominant','recessive','overdominant', \n                     'log-additive', 'all' or any combination of them")

    SNPs <- attr(data, "label.SNPs")
    lab.model <- c("codominant","dominant","recessive","overdominant","log-additive")

    if (quantitative)
     {
      if(!missing(nperm))
        stop("permutation approach is not yet implemented for quantitative traits")
  
      out <- data.frame(pvalTest(dataSNPs, Y = varDep, quantitative = quantitative,
        type = type, genotypingRate = genotypingRate))
     }

    else
     { 

# VM no calcula bien p si x==0

#      pvalG<-function(x,df){
#        pval <- x
#        pval[pval > 0] <- 1 - pchisq(pval[pval > 0], df)
#        pval
#      }

#      pvalG<-function(x,df) {pchisq(x, df, lower.tail=FALSE)}


# JRG para controlar los "calling rate" y "monomorphics"
       pvalG <- function(x, df) {
           if (x<0 | is.na(x))
             return(x)
           else
             return(pchisq(x, df, lower.tail = FALSE))
        }

      varDep<-as.numeric(as.factor(varDep))-1

# VM
   
#      aux<-cbind(varDep,dataSNPs)
#      o <- lapply(aux,as.numeric)
#      aux2<-matrix(unlist(o),nrow=nrow(aux),ncol=ncol(aux))
#      aux2[is.na(aux2)]<-0

# R 2.15
#      aux2<-matrix(.Internal(unlist(dataSNPs, FALSE, FALSE)),nrow=nrow(dataSNPs),ncol=ncol(dataSNPs))
#	   aux2<-cbind(varDep,aux2) 

      aux2<-matrix(unlist(cbind(varDep,dataSNPs), FALSE, FALSE),nrow=nrow(dataSNPs),ncol=ncol(dataSNPs)+1)
     
      aux2[is.na(aux2)]<-0
      nr<-nrow(aux2)
      nc<-ncol(aux2)

  
      if (type[1]==1)
        ngeno<-3
      else
        ngeno<-2

      out0<-.Fortran("WGassociation",as.integer(nr),as.integer(nc),as.integer(ngeno),
            as.integer(type[1]),as.integer(genotypingRate),as.integer(as.matrix(aux2)),
            stat=as.double(rep(0,nc-1)),PACKAGE = "SNPassoc")

      pval <- sapply(out0$stat, pvalG, df=ngeno-1)
  
      out<-data.frame(rep(NA,nc-1),pval)
      out[out[,2]==-1,1]<-"Monomorphic"
      out[out[,2]==-2,1]<-"Genot Error"
      out[out[,2]==-1 | out[,2]==-2,2]<-NA
      dimnames(out)[[1]]<-attr(data, "label.SNPs")

      if (length(type)>1) {
# VM
      for (i in type[-1]) {
          if (i==1)
           ngeno<-3
          else
           ngeno<-2
         
        out0.i<-.Fortran("WGassociation",as.integer(nr),as.integer(nc),as.integer(ngeno),
            as.integer(i),as.integer(genotypingRate),as.integer(as.matrix(aux2)),
            stat=as.double(rep(0,nc-1)),PACKAGE = "SNPassoc")

        pval <- sapply(out0.i$stat, pvalG, df=ngeno-1)
        pval[pval==-1 | pval==-2]<-NA
        out<-cbind(out,pval)
       }
      }
     }

     if(!missing(nperm)) {
       
       if (nperm<1)
         stop("number of permutation should be greater than 0")

       if (length(type)>1)
         warning("permutation test results only for the first genetic model ")   
       if (type[1]==1)
           ngeno<-3
       else
           ngeno<-2

       perm<-.Fortran("permutation",as.integer(nperm),
            as.integer(nr),as.integer(nc),as.integer(ngeno),as.integer(type[1]),
            as.integer(genotypingRate),as.integer(as.matrix(aux2)),
            stat=as.double(matrix(0,nrow=nc-1,ncol=nperm)),PACKAGE = "SNPassoc")
       ansPerm0<-lapply(perm$stat,pvalG,df=ngeno-1)
       ansPerm0[ansPerm0==-1 | ansPerm0==-2]<-NA    
       pvalPerm<-matrix(unlist(ansPerm0),nrow=nc-1,ncol=nperm)
      }


    if (max(type)==6)
     names(out)<-c("comments",lab.model)
    else
     names(out)<-c("comments",lab.model[type])

    for (i in 2:ncol(out)) out[, i] <- as.numeric(as.character(out[,i]))
 
    cost<-proc.time()-ptm
    cat("The program took", round(cost[3],2), "seconds \n")
     
    attr(out, "label.SNPs") <- attr(data, "label.SNPs")
    attr(out, "models") <- type
    attr(out, "quantitative") <- quantitative
    attr(out, "pvalues") <- out
    attr(out, "gen.info") <- attr(data, "gen.info")
    attr(out, "whole") <- attr(data, "whole")
    attr(out, "colSNPs") <- attr(data, "colSNPs")
    attr(out, "fast")<-TRUE
    if (!missing(nperm))
     attr(out, "pvalPerm") <- pvalPerm
    class(out) <- c("WGassociation","data.frame")
    out
}

