#' Get gene symbol from a list of SNPs
#' 
#' @description 
#' The getGeneSymbol function searches for the genes associated with the SNPs
#' for a given chromosome at a given position. To perform the annotation this
#' function needs to have the BioConductor libraries: S4Vectors, GenomicRanges, 
#' IRanges, VariantAnnotation and  org.Hs.eg.db libraries. 
#' As well as the library with the reference genome database, by default needs 
#' TxDb.Hsapiens.ICSC.hg19.KnownGene library
#' @param x data.frame containing: SNP name, chromosome and genomic position.
#' @param snpCol column of x having the SNP name. Default is 1.
#' @param chrCol column of x having the SNP chromosome. Default is 2.
#' @param posCol column of x having the SNP position. Default is 3.
#' @param db reference genome. Default is 'TxDb.Hsapiens.UCSC.hg19.knownGene'
#' @return a data.frame having initial information and gene symbol


getGeneSymbol <- function(x, snpCol=1, chrCol=2, posCol=3, 
                          db="TxDb.Hsapiens.UCSC.hg19.knownGene") {
    
    
    if (!requireNamespace("S4Vectors")) {
        stop( "Package \"S4Vectors\" must be installed to use this function.",
              call. = FALSE )
    }
    if (!requireNamespace("GenomicRanges")) {
        stop( "Package \"db <- eval(parse(text = db))\" must be installed to use this function.",
              call. = FALSE )
    }
    if (!requireNamespace("IRanges")) {
        stop( "Package \"IRanges\" must be installed to use this function.",
              call. = FALSE )
    }
    if (!requireNamespace("VariantAnnotation")) {
        stop( "Package \"VariantAnnotation\" must be installed to use this function.",
              call. = FALSE )
    }
    if (!requireNamespace("org.Hs.eg.db")) {
        stop( "Package \"org.Hs.eg.db\" must be installed to use this function.",
              call. = FALSE )
    }
    
    if (!requireNamespace(db)) {
        stop( "Reference genome must be installed to use this function.",
              call. = FALSE )
    } else {
        library(db, character.only=TRUE)
    }
    
    db <- eval(parse(text = db))    
    
    df <- x[,c(snpCol, chrCol, posCol)]
    names(df) <- c("rsid", "chr", "pos")
    
    target <- with(df,
                   GenomicRanges::GRanges( seqnames = S4Vectors::Rle(chr),
                                           ranges   = IRanges::IRanges(pos, end=pos, names=rsid),
                                           strand   = S4Vectors::Rle(GenomicRanges::strand("*"))))
    
    loc <- VariantAnnotation::locateVariants(target, db,
                                             VariantAnnotation::AllVariants())
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

