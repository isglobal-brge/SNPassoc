plot.WGassociation <- function(x, ...){
  if (!inherits(x, "WGassociation")) 
    stop("x must be an object of class 'WGassociation'")
  
  xx <- data.frame(SNP=rownames(x), data.frame(x)[,2:ncol(x)])
  names(xx)[2:ncol(x)] <- gsub("log.", "", names(x)[2:ncol(x)] ) 
  # names(xx)[6] <- "additive"
  dat <- tidyr::gather(xx, key="model", value="p.value", -"SNP")
  # dat$model <- factor(dat$model, 
  #                     levels=c("codominant", "dominant",
  #                              "recessive", "overdominant",
  #                              "additive"))
  dat$model <- factor(dat$model, 
                      levels=unique(names(xx)[2:ncol(xx)]))
  
  plt <- ggplot(dat, aes(x=.data$SNP, y=-log10(.data$p.value))) + 
    geom_point() +
    xlab("SNPs") + ylab(expression(-log[10]("p-value"))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(aes(yintercept = -log10(0.05), 
               linetype = "Nominal"), 
               colour="blue") +
    geom_hline(aes(yintercept = -log10(0.05/nrow(xx)),
               linetype = "Bonferroni"), 
               colour="red") + 
   
    scale_linetype_manual(name = "Significance", values = c(2, 2), 
                          guide = guide_legend(override.aes = list(color = c("red", "blue"))))
    
  plt + facet_wrap( ~ model, ncol=1) +
    theme(strip.text.x = element_text(face="bold")) +
    theme(legend.position="top", legend.direction="horizontal") 
}
