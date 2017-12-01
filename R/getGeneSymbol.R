getGeneSymbol <- function(x, snpCol=1, chrCol=2, posCol=3, 
                          db=TxDb.Hsapiens.UCSC.hg19.knownGene, ...) {

#
# x is a data.frame containing three columns:
#    rsid: rs number
#    chr: chromosome (annotated as chr1, chr2, chr3, ..., chrX, chrY)
#    pos: genomic position
# By default they are located in 1st, 2nd and 3rd column - they can be changes by using
# snpCol, chrCol, posCol  

# NOTE: human genome version in argument db should match with the annotation used in 
#  the 'x' argument

# Acknowledgment: Source code modified from Adai's post 

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

Symbol2id <- as.list( get("org.Hs.egSYMBOL2EG", envir=.GlobalEnv ))
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

