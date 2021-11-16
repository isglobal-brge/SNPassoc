dscore.character <- function(x, ...){
 rng <- BSgenome::snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh38::
                               SNPlocs.Hsapiens.dbSNP144.GRCh38, ids=x)
 # genome(rng) <- "GRCh38" # We need this ???
 
 mafs <- GenomicScores::gscores(MafDb.1Kgenomes.phase3.GRCh38::
                                    MafDb.1Kgenomes.phase3.GRCh38, 
                                rng, 
                                pop=c("EUR_AF"))
 
 refAlleles <- as.character(getSeq(BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens, rng))
 altAlleles <- DNA_BASES[(match(refAlleles, DNA_BASES)) %% 4 + 1]

 snpInfo <- data.frame(mafs, refAlleles, altAlleles) 
 mafs <- snpInfo$EUR_AF
 ans <- dscore(mafs)
 attr(ans, "MAFs") <- snpInfo
 ans
}
