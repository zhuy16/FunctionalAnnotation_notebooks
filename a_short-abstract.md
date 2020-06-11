# How to use clusterProfiler to do functional annotation of gene sets
clusterProfiler is a collection of different packages for functional annotation, these packages have a common theme. In the following summary, I try to abstract these functionality with an abstract review and a more detailed summary. All these can be a helpful introduction before going into this wonderful and sophisticated R book illustration https://yulab-smu.github.io/clusterProfiler-book/chapter1.html

Here below the line is a abstract overview of all functions in the clusterProfiler.
---
---
There are common themes of the various functions in the clusterProfiler and related package. This summary attempt to extract the common themes. 

### I. you will need the clusterProfiler functions to perform the functional annotation. 
A these functionns are pvalue dependent, you will use a small list of differentially expressed gennes identified by certain cutoff: 
```
enricher
enrichKEGG
enrichMKEGG 
# KEGG Module is a collection of manually defined function units, sometimes have more straightforward interpretation.
enrichGO 
enrichPathway
enrichMeSH (library(meshes))
enrichDO
enrichNCG
enrichDGN
enrichDGNv
```
as a related function, you could list all related GO terms using 
```
groupGO
```

B these functionns are rank based GSEA test, you will use a whole ranked list of fold-change values, with the names of the list as genes. You don't need p-value to start with and don't need a cut off. 
```
GSEA
gseKEGG
gseMKEGG
gseGO
gsePathway
gseMeSH
gseDO
gseNCG
gseDGN
```

### II, for database files, .gmt files contain information on which gene belong to which gene list. There are several sources of .gmt files. 

download from: http://data.wikipathways.org/20200510/gmt/ 
Or use the msigdbr package:
```
library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)
```
or for cell markers, from this human cell marker database
```
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
   tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
   dplyr::select(cellMarker, geneID) %>%
   dplyr::mutate(geneID = strsplit(geneID, ', '))
cell_markers
```
or for wikipathways, use the [rWikiPathway package](https://bioconductor.org/packages/release/bioc/manuals/rWikiPathways/man/rWikiPathways.pdf)

or for KEGG, the program automatically download from the site. specify the organism = "hsa" or "ko" etc.
[KEGG](http://www.genome.jp/kegg/catalog/org_list.html) or [KEGG Module](http://www.genome.jp/kegg/module.html)

### III, visualization of the results
```
# barplot
barplot(enricherRes, showCategory=20)  #showCategory: the number of go terms you want to show on the plot

#dotplots
edo2 <- gseNCG(geneList, nPerm=10000)
p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
p2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
plot_grid(p1, p2, ncol=2)

# cnetplot
options(repr.plot.width=16, repr.plot.height=24)
cnetplot(edox, foldChange=geneList, showCategory = 10,colorEdge=TRUE,circular=TRUE, node_label="all")
### more plots can be draw in different ways.
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

# Enrichment Map
options(repr.plot.width=8, repr.plot.height=6)
emapplot(xx,pie="count", pie_scale=2, layout="kk",legend_n=4)ridgeplot(edo2) # show the fold change of core genes in the GSEA analysis.

upsetplot(kk2) # this showes the overlaping of different genesets. (for GSEA it showes 
ridgeplot(edo2)
gseaplot2(edo2, geneSetID = 2, title = edo2$Description[2]) # geneSetID are the row number in the edo2 object

```
### IV, multiple gene set visualization (themes)
```
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

# visualization
dotplot(ck,showCategory = 20) #showCategory: the number of go terms you want to show on the plot

# visualization with different panels:
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)
```
---
---

**The following will be a more detailed summary** of each section.
---
## I. you will need the clusterProfiler functions to perform the functional annotation. 
1. Pvalue dependent
```
ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa')


ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)
               
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)               

```

2. Rank dependent 
```
ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
head(ewp2)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'mcc',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
```               

```
library(magrittr)
library(clusterProfiler)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
```

In addition, similar analysis using other databases, 
ReactomePA(Yu and He 2016) uses Reactome as a source of pathway data. The function call of **enrichPathway** and **gsePathway** in **ReactomePA**.

of gene list or whole expression profile using **MeSH annotation**. Data source from gendoo, gene2pubmed and RBBH are all supported. User can selecte interesting category to test. All 16 categories are supported. The analysis supports >70 species listed in MeSHDb BiocView.
```
library(meshes)
data(geneList, package="DOSE")
de <- names(geneList)[1:100]
x <- enrichMeSH(de, MeSHDb = **"MeSH.Hsa.eg.db"**, database='gendoo', category = 'C')
head(x)
```

Functional analysis using NGS data (eg, RNA-Seq and ChIP-Seq) can be performed by linking coding and non-coding regions to coding genes via **ChIPseeker**(Yu, Wang, and He 2015) package, which can annotates genomic regions to their nearest genes, host genes, and flanking genes respectivly. In addtion, it provides a function, seq2gene, that simultaneously considering host genes, promoter region and flanking gene from intergenic region that may under control via cis-regulation. This function maps genomic regions to genes in a many-to-many manner and facilitate functional analysis. For more details, please refer to ChIPseeker(Yu, Wang, and He 2015).
For details of ChIPseeker analysis, visit https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

Chapter 4 Disease analysis
DOSE(Yu et al. 2015) supports Disease Ontology (DO) Semantic and Enrichment analysis. The **enrichDO** function is very useful for identifying disease association of interesting genes, and function **gseDO** function is designed for gene set enrichment analysis of DO.

In addition, DOSE also supports enrichment analysis of **Network of Cancer Gene (NCG)**(A. et al. 2016) and **Disease Gene Network**(Janet et al. 2015), please refer to the DOSE vignettes.
```
x <- enrichDO(gene          = gene,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              **universe***      = names(geneList),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
```
```
enrichNCG
```
**DisGeNET**(Janet et al. 2015) is an integrative and comprehensive resources of gene-disease associations from several public data sources and the literature. It contains gene-disease associations and snp-gene-disease associations.

The enrichment analysis of disease-gene associations is supported by the `enrichDGN` function and analysis of snp-gene-disease associations is supported by the `enrichDGNv` function.

```
snp <- c("rs1401296", "rs9315050", "rs5498", "rs1524668", "rs147377392",
         "rs841", "rs909253", "rs7193343", "rs3918232", "rs3760396",
         "rs2231137", "rs10947803", "rs17222919", "rs386602276", "rs11053646",
         "rs1805192", "rs139564723", "rs2230806", "rs20417", "rs966221")
dgnv <- enrichDGNv(snp)
```

For Gsea test:
```
y <- gseDO(geneList,
           nPerm         = 100,
           minGSSize     = 120,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)

           
ncg <- gseNCG(geneList,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = FALSE)
ncg <- gseNCG(geneList,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = FALSE)
```

## II, for database files, .gmt files contain information on which gene belong to which gene list. There are several sources of .gmt files. 
One can 

download from: http://data.wikipathways.org/20200510/gmt/ 

```
wp2gene <- read.gmt("./wikipathways-20200510-gmt-Homo_sapiens.gmt")

wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")

wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
```
This one donesn't have gene symbols, therefore could be problematic 

You need to convert the ENTREZID to symbol. 

```
#library(clusterProfiler)

#library(org.Hs.eg.db)

geneLUT <- bitr(wpid2gene$gene,  fromType = "ENTREZID",
                        toType = c("ENSEMBL", "SYMBOL"),
                        OrgDb = org.Hs.eg.db)

### this will generate a lookuptable for the ENTREZID, the value as SYMBOL.

x=geneLUT$SYMBOL

names(x)=geneLUT$ENTREZID

y=x[wpid2gene$gene] # this will generate 

wpid2gene$gene=y

wpid2gene=na.omit(wpid2gene)
```
Or use the msigdbr package:
```
library(msigdbr)

msigdbr_show_species()

m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)

head(m_t2g)

wpid2gene=m_t2g

colnames(wpid2gene)=c("term","gene")

head(wpid2gene)
```
# III VISUALIZATION
```
# barplot: Bar plot is the most widely used method to visualize enriched terms. It depicts the enrichment scores (e.g. p values) and gene count or ratio as bar height and color.

barplot(enricherRes, showCategory=20)  #showCategory: the number of go terms you want to show on the plot
# dotplot
edo2 <- gseNCG(geneList, nPerm=10000)
p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
p2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
plot_grid(p1, p2, ncol=2)


# cnetplot

# Both the barplot and dotplot only displayed most significant enriched terms, while users may want to know **which genes are involved in these significant terms**. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, **we developed cnetplot function to extract the complex association**. The cnetplot depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network. GSEA result is also supported with **only core enriched genes displayed.**

options(repr.plot.width=16, repr.plot.height=24)
cnetplot(edox, foldChange=geneList, showCategory = 10,colorEdge=TRUE,circular=TRUE, node_label="all")
### more plots can be draw in different ways.
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

##### Enrichment Map
options(repr.plot.width=8, repr.plot.height=6)
emapplot(xx,pie="count", pie_scale=2, layout="kk",legend_n=4)ridgeplot(edo2) # show the fold change of core genes in the GSEA analysis.

upsetplot(kk2) # this showes the overlaping of different genesets. (for GSEA it showes 
ridgeplot(edo2)
gseaplot2(edo2, geneSetID = 2, title = edo2$Description[2]) # geneSetID are the row number in the edo2 object
```

# IV chapter 11, Biological theme comparison

clusterProfiler was developed for biological theme comparison(Yu et al. 2012), and it provides a function, **compareCluster**, to automatically calculate enriched functional categories of each gene clusters.
```
data(gcSample)

# gcSample is a list of gene vectors of ENTREZID named "X1","X2","X3"...
lapply(gcSample, head)

The input for geneCluster parameter should be a named list of gene IDs. To speed up the compilation of this document, we set use_internal_data = TRUE.

ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck)) # 
```

### or groups of genes can be identified from interaction of factors.
**compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")**

```
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))
```
### visualization
```
dotplot(
  object,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 10,
  size = NULL,
  split = NULL,
  font.size = 12,
  title = "",
  ...
)

dotplot(ck)

dotplot(ck,showCategory = 20) #showCategory: the number of go terms you want to show on the plot

visualization with different panels:
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)
```
---
