# Define a function to do automatic annotation. 
## fc is a vector of fold change values of all genes, named each value with it's gene ID (not as a symbl, but as 'ENTREZID')
## then you can run this function below:

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

