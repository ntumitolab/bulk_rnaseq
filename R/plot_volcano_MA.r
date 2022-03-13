library(ggplot2)
library(scales)
library(glue)
library(dplyr)
library(grid) # grob
library(ggrepel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./plotting.r")

INDIR <- "../../tables/S1_DEGs"
OUTDIR <- "../../plots/DEG_IGF1"
COMPR <- c("84Q_15Q", "84Q+IGF1_15Q", "84Q+IGF1_84Q")

# col names
X <- "baseMean"
Y <- "log2FoldChange"
P <- "padj"
XLABELS <- "mean of normalized counts"
YLABELS <- "log2 fold change"
PLABELS <- "p-value"

n <- sapply(COMPR, plotMAs, inDir=INDIR, outDir=OUTDIR,
            x=X, y=Y, p=P, 
            xlab=XLABELS, ylab=YLABELS, yThres=log2(1.5), pThres=0.1)

n <- sapply(COMPR, plotVolcano, 
            inDir=INDIR, outDir=OUTDIR,
            anotFn="../../data/gene_symbols.tsv", 
            y=Y, p=P, xlab=YLABELS, ylab=PLABELS, yThres=log2(1.5), pThres=0.1)
