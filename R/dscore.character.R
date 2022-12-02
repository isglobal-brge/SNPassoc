dscore.character <- function(x, ...){
 rng <- BSgenome::snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh38::SNPlocs.Hsapiens.dbSNP144.GRCh38, ids=x)
 genome(rng) <- "GRCh38.p2"
 
 mafs <- GenomicScores::gscores(MafDb.1Kgenomes.phase3.GRCh38::MafDb.1Kgenomes.phase3.GRCh38, rng, pop=c("EUR_AF"))
 
 gg <- BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens
 genome(gg) <- "GRCh38.p2"
 
 refAlleles <- as.character(getSeq(gg, rng))
 altAlleles <- DNA_BASES[(match(refAlleles, DNA_BASES)) %% 4 + 1]

 snpInfo <- data.frame(mafs, refAlleles, altAlleles) 
 snpInfo$EUR_AF[is.na(snpInfo$EUR_AF)] <- 0.001
 mafs <- snpInfo$EUR_AF
 ans <- dscore(mafs)
 attr(ans, "MAFs") <- snpInfo
 ans
}