# Define a function to do automatic annotation. 
## fc is a vector of fold change values of all genes, named each value with it's gene ID (not as a symbl, but as 'ENTREZID')
## then you can run this function below for GSEA tests:

FuncAnno=function(fc){
    cutoff=1
# 1. KEGG
    #' 
    ## -----------------------------------------------------------------------------
    gse_KEGG=gseKEGG(geneList=fc,organism='hsa',nPerm=1000,minGSSize=120,pvalueCutoff = cutoff,verbose=FALSE)
    gse_KEGG <- setReadable(gse_KEGG, 'org.Hs.eg.db', 'ENTREZID')
    # head(gse_KEGG)
# 234. GO
    gse_GO_MF=setReadable(gseGO(geneList=fc,OrgDb=org.Hs.eg.db,ont="MF",nPerm=1000,minGSSize=100,maxGSSize=500,,pvalueCutoff = cutoff,verbose=FALSE), 'org.Hs.eg.db', keyType="ENTREZID")
    gse_GO_CC=setReadable(gseGO(geneList=fc,OrgDb=org.Hs.eg.db,ont="CC",nPerm=1000,minGSSize=100,maxGSSize=500,,pvalueCutoff = cutoff,verbose=FALSE), 'org.Hs.eg.db', keyType="ENTREZID")
    gse_GO_BP=setReadable(gseGO(geneList=fc,OrgDb=org.Hs.eg.db,ont="BP",nPerm=1000,minGSSize=100,maxGSSize=500,pvalueCutoff = cutoff,verbose=FALSE), 'org.Hs.eg.db', keyType="ENTREZID")
# 5. hallmark
    ## -----------------------------------------------------------------------------
    h=msigdbr(species = "Homo sapiens", category = "H") 
    wpid2gene=h %>% dplyr::select(gs_name, entrez_gene) 
    gse_H <- GSEA(fc, TERM2GENE = wpid2gene,pvalueCutoff = cutoff)
    gse_H <- setReadable(gse_H, 'org.Hs.eg.db', 'ENTREZID')
# 6. Trans
    # wpid2gene <- read.gmt("de/c3.tft.v7.1.entrez.gmt")
    c3=msigdbr(species = "Homo sapiens", category = "C3") 
    wpid2gene=c3 %>% dplyr::select(gs_name, entrez_gene) 
    gse_Trans=GSEA(fc, TERM2GENE = wpid2gene,pvalueCutoff = cutoff, verbose=FALSE)
    gse_Trans=setReadable(gse_Trans, 'org.Hs.eg.db', keyType="ENTREZID")
    
# 7. Disease ontology'
    gse_DO=gseDO(fc,nPerm=100,minGSSize=120,pvalueCutoff = cutoff,verbose=FALSE)
    gse_DO=setReadable(gse_DO, 'org.Hs.eg.db', keyType="ENTREZID")
    
# 8. Disease Gene Networks.
    gse_DGN <- gseDGN(fc,nPerm = 100, minGSSize = 120,pvalueCutoff = cutoff, verbose = FALSE)
    gse_DGN <-setReadable(gse_DGN, 'org.Hs.eg.db', keyType="ENTREZID")    
# 12. immune signature
    m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
      dplyr::select(gs_name, entrez_gene)
    gse_Imm <- GSEA(fc, TERM2GENE = m_t2g,pvalueCutoff = cutoff)
       
    enrichList=list(
        gse_KEGG,
        gse_GO_BP,
        gse_GO_CC,
        gse_GO_MF,
        gse_H,
        gse_Trans,
        gse_DO,
        gse_DGN,
        gse_Imm
        )

    names(enrichList)=c(
        "gse_KEGG",
        "gse_GO_BP",
        "gse_GO_CC",
        "gse_GO_MF",
        "gse_H",
        "gse_Trans",
        "gse_DO",
        "gse_DGN",
        "gse_Imm"
    )
    return(enrichList)
}

# for a complete task including GO enrichments, here is another function.
## in this case the 'de_file' is a list of 3 items. the 2nd and 3rd item are used here. 
## the 2nd element is the fc vector as on the top. 
## the 3rd item is data frame containing all the differentially expressed genes (symbls) as the rowname. 
FuncAnn=function(de_file){ 
    de2gse=readRDS(de_file)
    #de_file
    library(clusterProfiler)
    library(DOSE)
    library(msigdbr)
    library(enrichplot)
    library(org.Hs.eg.db)

    #' 
    ## -----------------------------------------------------------------------------
    
    fc=de2gse[[2]]
    de=de2gse[[3]]
    deid=bitr(rownames(de),fromType="SYMBOL",toType=c("ENTREZID"),OrgDb=org.Hs.eg.db)
    deg=deid[,2]

#     head(fc)
#     head(deg)
    #' 
    ## -----------------------------------------------------------------------------
    xid=bitr(names(fc),fromType="SYMBOL",toType=c("ENTREZID"),OrgDb=org.Hs.eg.db)
    names(fc)=xid[,2]
#     fc
    #' # 1. KEGG
    #' 
    ## -----------------------------------------------------------------------------
    gse_KEGG=gseKEGG(geneList=fc,organism='hsa',nPerm=1000,minGSSize=120,pvalueCutoff=0.01,verbose=FALSE)
    gse_KEGG <- setReadable(gse_KEGG, 'org.Hs.eg.db', 'ENTREZID')
    # head(gse_KEGG)
    enrich_KEGG=enrichKEGG(gene=deg,organism='hsa',pvalueCutoff=0.05)
    enrich_KEGG=setReadable(enrich_KEGG, 'org.Hs.eg.db', keyType="ENTREZID")
    # head(enrich_Kegg)


    #' # 2. Wikipath, 
    library(msigdbr)
    wp2gene <- read.gmt("de/wikipathways-20200610-gmt-Homo_sapiens.gmt")
    wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
    wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
    wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

    gse_Wiki <- GSEA(fc, TERM2GENE = wpid2gene, pvalueCutoff=0.01, TERM2NAME = wpid2name)

    gse_Wiki <- setReadable(gse_Wiki, 'org.Hs.eg.db', keyType="ENTREZID")

    enrich_Wiki <- enricher(deg, TERM2GENE = wpid2gene,  pvalueCutoff = 0.05, TERM2NAME = wpid2name)
    enrich_Wiki <- setReadable(enrich_Wiki, 'org.Hs.eg.db', keyType="ENTREZID")
#     head(gse_Wiki)
#     head(enrich_Wiki)
    #' # 3. GO
    gse_GO_MF=setReadable(gseGO(geneList=fc,OrgDb=org.Hs.eg.db,ont="MF",nPerm=1000,minGSSize=100,maxGSSize=500,pvalueCutoff =0.01,verbose=FALSE), 'org.Hs.eg.db', keyType="ENTREZID")
    gse_GO_CC=setReadable(gseGO(geneList=fc,OrgDb=org.Hs.eg.db,ont="CC",nPerm=1000,minGSSize=100,maxGSSize=500,pvalueCutoff =0.01,verbose=FALSE), 'org.Hs.eg.db', keyType="ENTREZID")
    gse_GO_BP=setReadable(gseGO(geneList=fc,OrgDb=org.Hs.eg.db,ont="BP",nPerm=1000,minGSSize=100,maxGSSize=500,pvalueCutoff =0.001,verbose=FALSE), 'org.Hs.eg.db', keyType="ENTREZID")
#?gseGO
    enrich_GO_MF=setReadable(enrichGO(gene=deg,universe=names(fc),OrgDb=org.Hs.eg.db,ont="MF",pAdjustMethod="BH",pvalueCutoff=0.05,readable=TRUE), 'org.Hs.eg.db', keyType="ENTREZID")

    enrich_GO_CC=setReadable(enrichGO(gene=deg,universe=names(fc),OrgDb=org.Hs.eg.db,ont="CC",pAdjustMethod="BH",pvalueCutoff=0.05,readable=TRUE), 'org.Hs.eg.db', keyType="ENTREZID")

    enrich_GO_BP=setReadable(enrichGO(gene=deg,universe=names(fc),OrgDb=org.Hs.eg.db,ont="BP",pvalueCutoff=0.1,readable=TRUE), 'org.Hs.eg.db', keyType="ENTREZID")
#     head(gse_GO_BP)
#     head(enrich_GO_BP)
    # head(gse_GO_CC)
    # head(enrich_GO_CC)
    # head(gse_GO_MF)
    # head(enrich_GO_MF)
    #' 
    #' # 4. Hallmark
    #' 
    ## -----------------------------------------------------------------------------
    h=msigdbr(species = "Homo sapiens", category = "H") 
    wpid2gene=h %>% dplyr::select(gs_name, entrez_gene) 

    gse_H <- GSEA(fc, TERM2GENE = wpid2gene,pvalueCutoff =0.01)
    gse_H <- setReadable(gse_H, 'org.Hs.eg.db', 'ENTREZID')
    enrich_H <- enricher(deg, TERM2GENE = wpid2gene,pvalueCutoff = 0.05)
    enrich_H <- setReadable(enrich_H, 'org.Hs.eg.db', keyType="ENTREZID")
#     head(gse_H)
#     head(enrich_H)
    #' # 5. Canonical Pathway'
    wpid2gene <- read.gmt("de/c2.cp.v7.1.entrez.gmt")
    gse_CP <- GSEA(fc, TERM2GENE = wpid2gene,pvalueCutoff=0.01)
    gse_CP <- setReadable(gse_CP, 'org.Hs.eg.db', keyType="ENTREZID")

    enrich_CP <- enricher(deg, TERM2GENE = wpid2gene,pvalueCutoff = 0.05)
    enrich_CP <- setReadable(enrich_CP, 'org.Hs.eg.db', keyType="ENTREZID")

    #' # 6. MeSH
    #' 
    ## -----------------------------------------------------------------------------
    library(meshes)
    gse_MESH=gseMeSH(fc, MeSHDb = "MeSH.Hsa.eg.db", database='gendoo', pvalueCutoff=0.05, category = 'C')
    gse_MESH=setReadable(gse_MESH, 'org.Hs.eg.db', keyType="ENTREZID")
    enrich_MESH=enrichMeSH(deg, MeSHDb = "MeSH.Hsa.eg.db", pvalueCutoff = 0.05, database='gendoo', category = 'C')
    enrich_MESH=setReadable(enrich_MESH, 'org.Hs.eg.db', keyType="ENTREZID")
   
    wpid2gene <- read.gmt("de/c3.tft.v7.1.entrez.gmt")
    gse_Trans=GSEA(fc, TERM2GENE = wpid2gene, pvalueCutoff=0.0000001, verbose=FALSE)
    gse_Trans=setReadable(gse_Trans, 'org.Hs.eg.db', keyType="ENTREZID")

#     head(gse_Trans)
    enrich_Trans <- enricher(deg, pvalueCutoff=0.05,TERM2GENE = wpid2gene,)
    enrich_Trans <- setReadable(enrich_Trans, 'org.Hs.eg.db', keyType="ENTREZID")

    #' # 7. Disease ontology'
    gse_DO=gseDO(fc,nPerm=100,minGSSize=120,pvalueCutoff =0.05,verbose=FALSE)
    gse_DO=setReadable(gse_DO, 'org.Hs.eg.db', keyType="ENTREZID")

    enrich_DO=enrichDO(gene=deg,ont="DO",pvalueCutoff =0.05,universe=names(fc),minGSSize=5,maxGSSize=500,readable=FALSE)
    enrich_DO=setReadable(enrich_DO, 'org.Hs.eg.db', keyType="ENTREZID")
    
 
    #' # 8. Disease Gene Networks.

    gse_DGN <- gseDGN(fc,nPerm = 100, minGSSize = 120, pvalueCutoff = 0.05, verbose = FALSE)
    gse_DGN <-setReadable(gse_DGN, 'org.Hs.eg.db', keyType="ENTREZID")
    

    enrich_DGN=enrichDGN(deg,pvalueCutoff = 0.05)
    head(setReadable(enrich_DGN, 'org.Hs.eg.db', keyType="ENTREZID"))
     
    #' 
    enrichList=list(
        gse_KEGG,
        gse_CP,
        gse_Wiki,
        gse_GO_BP,
        gse_GO_CC,
        gse_GO_MF,
        gse_H,
        gse_MESH,
        gse_Trans,
        gse_DO,
        gse_DGN,

        enrich_KEGG,
        enrich_CP,
        enrich_Wiki,
        enrich_GO_BP,
        enrich_GO_CC,
        enrich_GO_MF,
        enrich_H,
        enrich_MESH,
        enrich_Trans,
        enrich_DO,
        enrich_DGN
    )
    #names(enrichList[1])
    names(enrichList)=c(
        "gse_KEGG",
        "gse_CP",
        "gse_Wiki",
        "gse_GO_BP",
        "gse_GO_CC",
        "gse_GO_MF",
        "gse_H",
        "gse_MESH",
        "gse_Trans",
        "gse_DO",
        "gse_DGN",

        "enrich_KEGG",
        "enrich_CP",
        "enrich_Wiki",
        "enrich_GO_BP",
        "enrich_GO_CC",
        "enrich_GO_MF",
        "enrich_H",
        "enrich_MESH",
        "enrich_Trans",
        "enrich_DO",
        "enrich_DGN"
    )
    #outFilename="HvVL.negative"

    #saveRDS(enrichList,outFilename)
    return(enrichList)
}
# nrow(enrichList[[20]])
