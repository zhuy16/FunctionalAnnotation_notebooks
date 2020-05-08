#  d <- read.csv("privateNotes/differVL.csv",header = T)
## assume that 1st column is ID
## 2nd column is fold change

#head(d)

## feature 1: numeric vector
#geneList <- d[,2]

## feature 2: named vector
#names(geneList) <- as.character(d[,1])

## feature 3: decreasing order
#geneList <- sort(geneList, decreasing = TRUE)

data(geneList, package="DOSE")
head(geneList)

gene <- names(geneList)[abs(geneList) > 2]
head(gene)

library(magrittr)
library(clusterProfiler)

#' # 3.2 WikiPathways analysis
#'

data(geneList, package="DOSE")

head(geneList)

gene <- names(geneList)[abs(geneList) > 2]

head(gene)

wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")



wp2gene <- read.gmt(wpgmtfile)

head(wp2gene)

#?tidyr::separate

wp2gene[1:3,]

#wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")


head(wp2gene)

wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)

ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)

head(ewp2)

#?GSEA()

#' You may want to convert the gene IDs to gene symbols, which can be done by setReadable function.

library(org.Hs.eg.db)
ewp <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
ewp2 <- setReadable(ewp2, org.Hs.eg.db, keyType = "ENTREZID")
head(ewp)

#' 3.3 Cell Marker
#'


library(vroom)
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))


head(cell_markers)


y <- enricher(gene, TERM2GENE=cell_markers, minGSSize=1)
DT::datatable(as.data.frame(y))

#' 3.4 MSigDb analysis
#'
library(msigdbr)
msigdbr_show_species()

m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame


m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)


em <- enricher(gene, TERM2GENE=m_t2g)

head(em)

library(clusterProfiler)
library(fgsea)

em2 <- GSEA(geneList, TERM2GENE = m_t2g)

head(em2)

#' We can test with other collections, for example, using C3 to test whether the genes are up/down-regulated by sharing specific motif.
#'
#' m_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
m_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em3 <- GSEA(geneList, TERM2GENE = m_t2g)
head(em3)


