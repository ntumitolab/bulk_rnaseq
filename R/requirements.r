install.packages("ggplot2")
install.packages("ggnewscale")
install.packages("glue")
# install.packages("argparse")
install.packages("readxl")
install.packages("stringr")
install.packages("enrichR")
install.packages("msigdbr")
install.packages("ggridges")
install.packages("docopt")
install.packages("rjson")
install.packages("here")

BiocManager::install("impute")
install.packages("WGCNA")
install.packages(c("impute","dynamicTreeCut","qvalue","flashClust","Hmisc"))



BiocManager::install("vctrs")
# optional
install.packages('ggfortify')


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")
BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")
BiocManager::install("enrichplot")

# optional (batch effect removal)
BiocManager::install("sva")
BiocManager::install("RUVSeq")

install.packages("metaRNASeq")

