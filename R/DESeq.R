#!/usr/bin/env Rscript
library(DESeq2)
library(glue)
library(dplyr)
library(ggplot2)
library(stringr)

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./plotting.r")
ADJ_METHOD <- "BH"
LFC_THRES <- 1
PVAL_THRES <- 0.05
TABLE_RESULT_DIR <- "../results/2022March/DEG/tables"
PLOT_RESULT_DIR <- "../results/2022March/DEG/plots"
comparisons <- list()
comparisons[[1]] <- c("2", "1")
comparisons[[2]] <- c("4", "1")
comparisons[[3]] <- c("5", "1")
comparisons[[4]] <- c("4", "2")
comparisons[[5]] <- c("5", "2")
comparisons[[6]] <- c("5", "4")
outliers <- c("")  # "R031", "R040", "R046", "R050", "R054"

sampleInfo <- read.table("../data/all_info/labels_5g.csv", sep='\t', header = TRUE)
# sampleInfo <- sampleInfo[-match(outliers, sampleInfo$sample),]
countTable <- read.table("../data/mg_counts_2203.tsv", sep='\t', header = TRUE)
countsData <- as.matrix(countTable[,2:dim(countTable)[2]])

colnames(countsData) <- str_match(colnames(countsData), "(.*)A")[,2]
rownames(countsData) <- countTable$gene

# countsData <- countsData[, -match(outliers, colnames(countsData))]

sampleInfo$group <- factor(sampleInfo$group)
sampleInfo$batch <- factor(sampleInfo$batch)
dds <- DESeqDataSetFromMatrix(countsData[, sampleInfo$sample], 
                              colData = sampleInfo, 
                              design = ~batch + group)
vsd <- vst(dds)
vsdResult <- plotPCA(vsd, intgroup=c("group"), returnData=T)
jpeg(file.path(PLOT_RESULT_DIR, "PCA_counts_vst_rvout.jpeg"), 
     width = 9, height = 9, units = "in",
     res = 600, quality = 100)
plotPCAc(vsdResult, condName = "group")
dev.off()

vsd <- vst(dds)
vsdResult <- plotPCA(vsd, intgroup=c("batch"), returnData=T)
jpeg(file.path(PLOT_RESULT_DIR, "PCA_counts_vst_batch_rvout.jpeg"), 
     width = 9, height = 9, units = "in",
     res = 600, quality = 100)
plotPCAc(vsdResult, condName = "batch")
dev.off()

dds <- DESeq(dds)
for (i in 1:length(comparisons))
{
  print(glue("{comparisons[[i]][1]} vs {comparisons[[i]][2]}"))
  res <- results(dds, contrast = c("group", comparisons[[i]]), pAdjustMethod = ADJ_METHOD)
  notna <- res[!is.na(res$padj), ]
  notna <- notna[(abs(notna$log2FoldChange) > LFC_THRES) & (notna$padj < PVAL_THRES), ]
  jpeg(file.path(PLOT_RESULT_DIR, glue("MA_{comparisons[[i]][1]}_{comparisons[[i]][2]}.jpeg")),
       width = 9, height = 9, units = "in",
       res = 300, quality = 100)
  DESeq2::plotMA(res, ylim=c(-3,3))
  dev.off()
  write.table(res, file.path(TABLE_RESULT_DIR, glue("DE_{comparisons[[i]][1]}_{comparisons[[i]][2]}.tsv")), sep = "\t", row.names = TRUE)
  write.table(data.frame(DGE =row.names(notna)), file.path(TABLE_RESULT_DIR, glue("DEG_list_{comparisons[[i]][1]}_{comparisons[[i]][2]}.tsv")), sep = "\t", row.names = TRUE)
}

X <- "baseMean"
Y <- "log2FoldChange"
P <- "padj"
XLABELS <- "mean of normalized counts"
YLABELS <- "log2 fold change"
PLABELS <- "adjusted p-value"

INDIR <- "../results/2022March/DEG/tables"
OUTDIR <- "../results/2022March/DEG/plots"
annotFilePath <- "../data/gene_symbol_2.tsv"
COMPR <- c("2_1", "4_1", "5_1", "4_2", "5_2", "5_4")

n <- sapply(COMPR, plotMAs, inDir=INDIR, outDir=OUTDIR,
            x=X, y=Y, p=P, 
            xlab=XLABELS, ylab=YLABELS, yThres=LFC_THRES, pThres=PVAL_THRES)

n <- sapply(COMPR, plotVolcano, 
            inDir=INDIR, outDir=OUTDIR,
            anotFn=annotFilePath, 
            y=Y, p=P, xlab=YLABELS, ylab=PLABELS, yThres=LFC_THRES, pThres=PVAL_THRES)
