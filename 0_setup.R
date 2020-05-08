install.packages("devtools")
devtools::install_github("rstudio/bookdown")
bookdown::publish_book(render = 'local')
getwd()
list.files()
setwd("bookdown-demo/")
bookdown::publish_book(render = 'local')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# not good, BiocManager::install("DOSE")

library(DOSE)

list.files()

BiocManager::install("GSEABase")
BiocManager::install("BiocGenerics")
packageVersion("BiocGenerics")
BiocManager::install("AnnotationDbi")


packageVersion("BioGenerics")   
#' BiocManager::install("org.Hs.eg.db") didn't install the newest version. 

install.packages("vroom")
install.packages("msigdbr")

BiocManager::install("pathview")

#' install from BiocManager always used old versions, so I download and installed it from source.
#' 
install.packages("yunhua/DOSE_3.14.0.tar.gz", repos = NULL, type="source")
packageVersion("DOSE")   
devtools::install_github('https://github.com/YuLab-SMU/clusterProfiler',dependencies = T)
devtools::install_github("GuangchuangYu/enrichplot")


# packageVersion("AnnotationDbi")
# [1] ‘1.44.0’
# > packageVersion("GSEABase")
# [1] ‘1.44.0’
# > packageVersion("org.Hs.eg.db")
# [1] ‘3.7.0’
# > packageVersion("vroom")
# [1] ‘1.2.0’
# > packageVersion("msigdbr")
# [1] ‘7.0.1’
# > packageVersion("clusterProfiler")
# [1] ‘3.17.0’
# > packageVersion("enrichplot")
# [1] ‘1.9.1’

# BiocManager::install("clusterProfiler.dplyr")
# install.packages("clusterProfiler.dplyr")
devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
install.packages("ggstance")

wget https://www.bioconductor.org/packages/release/bioc/src/contrib/DOSE_3.14.0.tar.gz
wget https://bioconductor.org/packages/release/bioc/src/contrib/fgsea_1.14.0.tar.gz
install.packages("DOSE_3.14.0.tar.gz", repos = NULL, type="source")
install.packages("fgsea_1.14.0.tar.gz", repos = NULL, type="source")