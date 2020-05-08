# Chapter 14 Useful utilities
# 14.1 bitr: Biological Id TranslatoR

x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)


library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db")
head(ids)

#' 14.1.1 bitr_kegg: converting biological IDs using KEGG API

data(gcSample)
hg <- gcSample[[1]]
head(hg)

eg2np <- bitr_kegg(hg, fromType='kegg', toType='ncbi-proteinid', organism='hsa')
head(eg2np)

bitr_kegg("Z5100", fromType="kegg", toType='ncbi-geneid', organism='ece')

bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')

bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')

#' 14.2 setReadable: translating gene IDs to human readable symbols
#' 
library(org.Hs.eg.db)
library(clusterProfiler)

data(geneList, package="DOSE")
de <- names(geneList)[1:100]
x <- enrichKEGG(de)
## The geneID column is ENTREZID
head(x, 3)

y <- setReadable(x, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
## The geneID column is translated to symbol
head(y, 3)

