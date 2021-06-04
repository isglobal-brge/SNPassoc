## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(comment="", warning=FALSE, message=FALSE, cache=TRUE)

## ----load_asthma--------------------------------------------------------------
data(asthma, package = "SNPassoc")
str(asthma, list.len=9)
asthma[1:5, 1:8]

## ----prepare_data-------------------------------------------------------------
library(SNPassoc)
asthma.s <- setupSNP(data=asthma, colSNPs=7:ncol(asthma), sep="")

## ----prepare_data2------------------------------------------------------------
idx <- grep("^rs", colnames(asthma))
asthma.s <- setupSNP(data=asthma, colSNPs=idx, sep="")

## ----classSNP-----------------------------------------------------------------
head(asthma.s$rs1422993)
class(asthma.s$rs1422993)

## ----summarySNP---------------------------------------------------------------
summary(asthma.s$rs1422993)

## ----plotsummarySNP, fig.cap="SNP summary. Bar chart showing the basic information of a given SNP"----
plot(asthma.s$rs1422993)

## ----plotsummarySNP2, fig.cap="SNP summary. Pie chart showing the basic information of a given SNP"----
plot(asthma.s$rs1422993, type=pie)

## ----descriptive--------------------------------------------------------------
summary(asthma.s, print=FALSE)

## ----plotMissing, fig.cap="Missing genotypes. Black squares shows missing genotuype information of asthma data example."----
plotMissing(asthma.s, print.labels.SNPs = FALSE)

## ----HWE----------------------------------------------------------------------
hwe <- tableHWE(asthma.s)
head(hwe)

## ----HWE_controls-------------------------------------------------------------
hwe2 <- tableHWE(asthma.s, casecontrol)

#SNPs is HWE in the whole sample but not controls
snpNHWE <- hwe2[,1]>0.05 & hwe2[,2]<0.05
rownames(hwe2)[snpNHWE]
hwe2[snpNHWE,]

## ----remove_SNPs--------------------------------------------------------------
snps.ok <- rownames(hwe2)[hwe2[,2]>=0.001]
pos <- which(colnames(asthma)%in%snps.ok, useNames = FALSE)
asthma.s <- setupSNP(asthma, pos, sep="")

## ----association--------------------------------------------------------------
association(casecontrol ~ rs1422993, data = asthma.s)

## ----maxrs1422993-------------------------------------------------------------
maxstat(asthma.s$casecontrol, asthma.s$rs1422993)

## ----mode_inheritance---------------------------------------------------------
association(casecontrol ~ rs1422993, asthma.s, model="dominant")

## ----adjusted-----------------------------------------------------------------
association(casecontrol ~ rs1422993 + country + smoke, asthma.s)

## ----stratified---------------------------------------------------------------
association(casecontrol ~ rs1422993 + survival::strata(gender), asthma.s)

## ----subset-------------------------------------------------------------------
association(casecontrol ~ rs1422993, asthma.s, 
                subset=country=="Spain")

## ----quantitative-------------------------------------------------------------
association(bmi ~ rs1422993, asthma.s) 

## ----morethan1SNP-------------------------------------------------------------
ans <- WGassociation(casecontrol, data=asthma.s)
head(ans)

## ----morethan1SNPadjusted, eval=FALSE-----------------------------------------
#  ans.adj <- WGassociation(casecontrol ~ country + smoke, asthma.s)
#  head(ans.adj)

## ----installGitHubVersion, eval=FALSE-----------------------------------------
#  devtools::install_github("isglobal-brge/SNPassoc")

## ----plotGWAS, fig.cap="Manhattan-type plots for different genetic models. P-values in -log10 scale to assess the association between case-control status and SNPs in the asthma example.", fig.height=8----
plot(ans)

## ----max-statistic------------------------------------------------------------
ans.max <- maxstat(asthma.s, casecontrol)
ans.max

## ----maxBonferroni------------------------------------------------------------
#minimum P-value across SNPs
min(ans.max["Pr(>z)",])

## ----OR_several_SNPs, results='hide'------------------------------------------
infoTable <- WGstats(ans)

## ----get_info_rs1422993-------------------------------------------------------
infoTable$rs1422993

## ----getNiceTable, eval=FALSE-------------------------------------------------
#  library(xtable)
#  out <- getNiceTable(ans[c("rs1422993", "rs184448")])
#  
#  nlines <- attr(out, "nlines")
#  hlines <- c(-1, -1, 0, cumsum(nlines+1), nrow(out), nrow(out))
#  
#  print(xtable(out, caption='Genetic association using
#                  different genetic models from asthma
#                  data example of rs1422993 and rs184448
#                  SNPs obtained with SNPassoc.',
#               label = 'tab-2SNPs'),
#        tabular.enviroment="longtable", file="tableSNPs",
#        floating=FALSE,  include.rownames = FALSE,
#        hline.after= hlines, sanitize.text.function=identity)

## ----snpxsmoke----------------------------------------------------------------
association(casecontrol ~ dominant(rs1422993)*factor(smoke), 
            data=asthma.s)

## ----snpxsnp------------------------------------------------------------------
association(casecontrol ~ rs1422993*factor(rs184448), 
            data=asthma.s, model.interaction = "dominant" )

## ----gxg_subset---------------------------------------------------------------
ans <- WGassociation(casecontrol, data=asthma.s)
mask <- apply(ans, 1, function(x) min(x, na.rm=TRUE)<0.1)
sig.snps <- names(mask[mask])
sig.snps
idx <- which(colnames(asthma)%in%sig.snps)
asthma.s2 <- setupSNP(asthma, colSNPs = idx, sep="")
ans.int <- interactionPval(casecontrol ~ 1, data=asthma.s2)
ans.int

## ----plotInf, fig.cap="Interaction plot. Interaction plot of SNPs significant al 10\\% significant level (see help of 'interactionPval'  function to see what is represented in the plot). ", fig.height=7, fig.width=7----
plot(ans.int)

## ----haplo_em-----------------------------------------------------------------
library(haplo.stats)
snpsH <- c("rs714588", "rs1023555",  "rs898070")
genoH <- make.geno(asthma.s, snpsH)
em <- haplo.em(genoH, locus.label = snpsH, miss.val = c(0, NA))
em

## ----hap_assoc----------------------------------------------------------------
trait <- asthma.s$casecontrol
mod <- haplo.glm(trait ~ genoH,           
                 family="binomial", 
                 locus.label=snpsH,
                 allele.lev=attributes(genoH)$unique.alleles,
                 control = haplo.glm.control(haplo.freq.min=0.05))   
intervals(mod)

## ----sliding_window_hap-------------------------------------------------------
snpsH2 <- labels(asthma.s)[6:15]
genoH2 <- make.geno(asthma.s, snpsH2)
haplo.score <- list()
for (i in 4:7) {
 trait <- asthma.s$casecontrol
 haplo.score[[i-3]] <- haplo.score.slide(trait, genoH2, 
                          trait.type="binomial",
                          n.slide=i,
                          simulate=TRUE,
                          sim.control=score.sim.control(min.sim=100,
                                       max.sim=200)) 
 }

## ----plotSliding, fig.cap="Sliding window approach. Results obtained of varying haplotype size from 4 up to 7 of  6th to the 15th SNP from asthma data example", fig.height=7, fig.width=7----
par(mfrow=c(2,2))
for (i in 4:7) {
    plot(haplo.score[[i-3]])
    title(paste("Sliding Window=", i, sep=""))
 }

## ----hap_assoc_bestH----------------------------------------------------------
snpsH3 <- snpsH2[4:7]
genoH3 <- make.geno(asthma.s, snpsH3)
mod <- haplo.glm(trait~genoH3,           
                 family="binomial", 
                 locus.label=snpsH3,
                 allele.lev=attributes(genoH3)$unique.alleles,
                 control = haplo.glm.control(haplo.freq.min=0.05))      
intervals(mod)

## ----lrt_hap------------------------------------------------------------------
lrt <- mod$lrt
pchisq(lrt$lrt, lrt$df, lower=FALSE)

## ----lrt_hap_Adj--------------------------------------------------------------
smoke <- asthma.s$smoke
mod.adj.ref <- glm(trait ~ smoke, family="binomial")
mod.adj <- haplo.glm(trait ~ genoH3 + smoke ,           
                 family="binomial", 
                 locus.label=snpsH3,
                 allele.lev=attributes(genoH3)$unique.alleles,
                 control = haplo.glm.control(haplo.freq.min=0.05))

lrt.adj <- mod.adj.ref$deviance - mod.adj$deviance
pchisq(lrt.adj, mod.adj$lrt$df, lower=FALSE)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

