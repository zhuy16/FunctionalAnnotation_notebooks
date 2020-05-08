#library(magrittr)
library(clusterProfiler)

data(gcSample)
lapply(gcSample, head)
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))

#' 11.1 Formula interface of compareCluster

mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))

#' 11.2 Visualization of profile comparison

dotplot(ck)
dotplot(formula_res)


dotplot(formula_res)

dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)


