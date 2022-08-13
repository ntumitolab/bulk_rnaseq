install.packages("remotes")
install.packages("VennDiagram")

remotes::install_github('alexvpickering/crossmeta')


install.packages("survival","samr","combinat")
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("limma")
biocLite("edgeR")
biocLite("DESeq2")
biocLite("Biobase")