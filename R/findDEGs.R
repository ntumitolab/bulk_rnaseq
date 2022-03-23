#!/usr/bin/env Rscript
library(DESeq2)
library(glue)
library(dplyr)
library(ggplot2)
library(stringr)
library(argparse)
library(tools)
library(readxl)


handleCountData <- function(file.dir, verbose){
  files <- list.files(path = file.dir)
  if (length(files) != 1){
    if (verbose) {
      print("Multiple files are detected in the input directory, trying to merge them into a data.frame..")
    }
    first.df <- T
    for (f in files){
      new.df <- read.csv(file.path(file.dir, f), header = T)
      colnames(new.df) <- c('gene', tools::file_path_sans_ext(f))
      if (first.df){
        merged.df <- new.df
        first.df <- F
      } else {
        merged.df <- merge(merged.df, new.df, by='gene', all=TRUE)
      }
    }
    row.names(merged.df) <- merged.df[,1]
    merged.df <- merged.df[,-c(1)]
    
  } else {
    merged.df <- read.csv(file.path(file.dir, files[[1]]), row.names = 1, header = T)
  }
  return(merged.df)
}

get.comparison <- function(groups, ct=NA){
  
}

# args for building a command line tool
parser <- argparse::ArgumentParser()

parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-i", "--input", help="Inpur folder storing raw count data")
parser$add_argument("-o", "--output", help="Output folder")
parser$add_argument("-l", "--label", help="A file containing samples' labels")
parser$add_argument("-p", "--pval", help="The p-value cutoff used to determine the DEGs", 
                    default=0.05,)
parser$add_argument("-f", "--foldchange", 
                    help="The absolute log2 fold change cutoff used to determine the DEGs", 
                    default=1,)
parser$add_argument("-m", "--method", default="D")
parser$add_argument("-a", "--adjMethod", 
                    help="The name of p-value adjusting method", 
                    default="BH")
parser$add_argument("-s", "--skip", 
                    help="A file listing the names of the skipped raw count data", 
                    default="BH")

args <- parser$parse_args()

# check and create the output folders
ifelse(!dir.exists(file.path(args$input, "tables")), 
       dir.create(file.path(args$input, "tables"), 
                  recursive = F), FALSE)
ifelse(!dir.exists(file.path(args$output, "plots")), 
       dir.create(file.path(args$output, "plots"), 
                  recursive = F), FALSE)

TABLE_RESULT_DIR <- file.path(args$output, "tables")
PLOT_RESULT_DIR <- file.path(args$output, "plots")

# deal with the label file
if (file_ext(args$label) == "csv") {
  sample.info <- read.csv(args$label)
} else if (file_ext(args$label) == "tsv") {
  sample.info <- read.delim(args$label)
} else if (file_ext(args$label) == "xlsx") {
  sample.info <- read.ex(args$label)
} else {
  sample.info <- read.table(args$label, header = T, sep = ",")
}

if (~("group" %in% colnames(sample.info))){
  stop("The label file should contain samples' groups")
}

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
