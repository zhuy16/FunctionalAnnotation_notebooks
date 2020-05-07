#' Chapter 6 KEGG analysis
#' 
library(clusterProfiler)
search_kegg_organism('ece', by='kegg_code')


ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
dim(ecoli)
head(ecoli)

#' 6.1 KEGG over-representation test

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

#' 6.2 KEGG Gene Set Enrichment Analysis
#' 
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

#' 6.3 KEGG Module over-representation test

mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa')

#' 6.4 KEGG Module Gene Set Enrichment Analysis
#' 
mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa')

