#' Get gene symbol from a list of SNPs
#' 
#' @param x data.frame containing: SNP name, chromosome and genomic position.
#' @param snpCol column of x having the SNP name. Default is 1.
#' @param chrCol column of x having the SNP chromosome. Default is 2.
#' @param posCol column of x having the SNP position. Default is 3.
#' @param db reference genome. Default is 'TxDb.Hsapiens.UCSC.hg19.knownGene'
#' @return a data.frame having initial information and gene symbol


getGeneSymbol <- function(x, snpCol=1, chrCol=2, posCol=3, 
                          db=TxDb.Hsapiens.UCSC.hg19.knownGene) {

    df <- x[,c(snpCol, chrCol, posCol)]
    names(df) <- c("rsid", "chr", "pos")
    
    target <- with(df,
                    GRanges( seqnames = Rle(chr),
                             ranges   = IRanges(pos, end=pos, names=rsid),
                             strand   = Rle(strand("*")) ) )
    
    loc <- locateVariants(target, db, AllVariants())
    names(loc) <- NULL
    info <- as.data.frame(loc)
    info$names <- names(target)[ info$QUERYID ]
    info <- info[ , c("names", "seqnames", "start", "end", "LOCATION", "GENEID", "PRECEDEID", "FOLLOWID")]
    info <- unique(info)
    
    Symbol2id <- as.list( org.Hs.eg.db::org.Hs.egSYMBOL2EG)
    id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
    names(id2Symbol) <- unlist(Symbol2id)
     
    x <- unique( with(info, c(levels(GENEID), levels(PRECEDEID), levels(FOLLOWID))) )
    
    info$GENESYMBOL <- id2Symbol[ as.character(info$GENEID) ]
    info$PRECEDESYMBOL <- id2Symbol[ as.character(info$PRECEDEID) ]
    info$FOLLOWSYMBOL <- id2Symbol[ as.character(info$FOLLOWID) ]
    ans <- info[, c("names", "seqnames", "start", "LOCATION", "GENESYMBOL", "GENEID")]
    names(ans)[3] <- "position"
    temp <- which(names(df)%in%c("chr","pos"))
    out <- merge(ans, df[, -temp, drop=FALSE], by.x="names", by.y="rsid")
    out
}

