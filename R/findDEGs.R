#!/usr/bin/env Rscript

'Do DE analysis using raw count data
Usage:
    findDEGs.R [--pval=<pval> --logfc=<logfc> --adjMethod=<adjMethod>] <input> <output> <metadata>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --pval=<pval>  Width of the output [default: 0.05]
    --logfc=<logfc>  The absolute log2 fold change cutoff used to determine the DEGs [default: 1]
    --adjMethod=<adjMethod> The name of p-value adjusting method [default: BH]

Arguments:
    input  inpur folder storing raw count data
    output  output folder
    metadata  A file containing samples information
' -> doc


library(docopt)
suppressMessages(library(DESeq2))
suppressMessages(library(glue))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(tools))
suppressMessages(library(readxl))


read_data <- function(file_name){
  loaded_data <- switch(file_ext(file_name),
                        csv = read.csv(file_name, row.names = 1),
                        tsv = read.delim(file_name, row.names = 1),
                        xlsx = readxl::read_xlsx(file_name),
                        xls = readxl::read_xls(file_name),
                        txt = read.csv(file_name, row.names = 1, header = T, sep = "\t"),
                        read.table(file_name, header = T, sep = ",")
  )
  return(loaded_data)
}


handle_count_data <- function(filedir, verbose){
  files <- list.files(path = filedir)
  if (length(files) != 1){
    if (verbose) {
      print("Multiple files are detected in the input directory, trying to merge them into a data.frame..")
    }
    first.df <- T
    for (f in files){
      new.df <- read.csv(file.path(filedir, f), header = T)
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
    if (verbose) {
      print("Only one file is detected in the input folder")
    }
    
    merged.df <- read_data(file.path(filedir, files[[1]]))
  }
  return(merged.df)
}

get_comparison <- function(groups, ct=NA){
  n_g = length(groups)
  comparisons <- list()
  
  if (is.na(ct)){
    comb_idx = combn(n_g, 2)
    n_comb = dim(comb_idx)[2]
    
    for (comb_i in 1: n_comb) {
      idx = comb_idx[, comb_i]
      comparisons[[length(comparisons) + 1]] <- c(groups[idx[2]], groups[idx[1]])
    }
  } else {
    if (!(ct %in% groups)) {
      stop("Assigned control label is not in the groups")
    }
    for (grp in groups) {
      if (grp == ct){
        next
      }
      comparisons[[length(comparisons) + 1]] <- c(grp, ct)
    }
  }
  comparisons
}

GROUP_COLORS <- c("#3367c7", "#f08300", "#17b055", "#5e17b0", "#9c1f1f")

plotPCAc <- function(pcaTable, condName){
  percentVar <- attr(pcaTable, "percentVar")
  xLabel <- glue("PC1 {round(percentVar[1] * 100, digits=2)} % variance")
  yLabel <- glue("PC2 {round(percentVar[2] * 100, digits=2)} % variance")
  p <- ggplot(pcaTable) + 
    theme_bw() +
    geom_point(aes_string(x="PC1", 
                          y="PC2", 
                          color=condName), 
               size=4) +
    geom_vline(xintercept=0, size=1, alpha=0.8) +
    geom_hline(yintercept=0, size=1, alpha=0.8) +
    labs(x = xLabel, y = yLabel) +
    scale_color_manual(name = condName, values=GROUP_COLORS[1:nlevels(pcaTable$group)])
  
  return(p)
}

# args for building a command line tool
args <- docopt(doc)

# check and create the output folders
ifelse(!dir.exists(file.path(args$output, "tables")), 
       dir.create(file.path(args$output, "tables"), 
                  recursive = T), FALSE)
ifelse(!dir.exists(file.path(args$output, "plots")), 
       dir.create(file.path(args$output, "plots"), 
                  recursive = T), FALSE)

TABLE_RESULT_DIR <- file.path(args$output, "tables")
PLOT_RESULT_DIR <- file.path(args$output, "plots")

# deal with the label file

print("loading count data")
counts_data <- as.matrix(handle_count_data(args$input, T))

print("loading metadata..")
print(args$metadata)
metadata <- read_data(args$metadata)
print("metadata loaded")
# 
# 
# if (~("group" %in% colnames(metadata))){
#   stop("The label file should contain samples' groups")
# }

metadata$group <- factor(metadata$group)
# metadata$batch <- factor(metadata$batch)

comparisons <- get_comparison(levels(metadata$group))
# outliers <- args$skip  # not implemented yet
# metadata <- metadata[-match(outliers, metadata$sample),]


# colnames(countsData) <- str_match(colnames(countsData), "(.*)A")[,2]
# countsData <- countsData[, -match(outliers, colnames(countsData))]

print("creating DESeqDataSet..")
dds <- DESeqDataSetFromMatrix(counts_data[, metadata$sample], 
                              colData = metadata, 
                              design = ~group)
vsd <- vst(dds)
vsdResult <- plotPCA(vsd, intgroup=c("group"), returnData=T)
jpeg(file.path(PLOT_RESULT_DIR, "PCA_counts_vst_rvout.jpeg"), 
     width = 9, height = 9, units = "in",
     res = 600, quality = 100)
plotPCAc(vsdResult, condName = "group")
dev.off()
# 
# vsd <- vst(dds)
# vsdResult <- plotPCA(vsd, intgroup=c("batch"), returnData=T)
# jpeg(file.path(PLOT_RESULT_DIR, "PCA_counts_vst_batch_rvout.jpeg"), 
#      width = 9, height = 9, units = "in",
#      res = 600, quality = 100)
# plotPCAc(vsdResult, condName = "batch")
# dev.off()

dds <- DESeq(dds)
for (i in 1:length(comparisons))
{
  print(glue("{comparisons[[i]][1]} vs {comparisons[[i]][2]}"))
  res <- results(dds, contrast = c("group", comparisons[[i]]), pAdjustMethod = args$adjMethod)
  notna <- res[!is.na(res$padj), ]
  notna <- notna[(abs(notna$log2FoldChange) > args$logfc) & (notna$padj < args$pval), ]
  write.table(res, file.path(TABLE_RESULT_DIR, glue("DE_{comparisons[[i]][1]}_{comparisons[[i]][2]}.tsv")), sep = "\t", row.names = TRUE)
  write.table(data.frame(DGE =row.names(notna)), file.path(TABLE_RESULT_DIR, glue("DEG_list_{comparisons[[i]][1]}_{comparisons[[i]][2]}.tsv")), sep = "\t", row.names = TRUE)
}
# 
# X <- "baseMean"
# Y <- "log2FoldChange"
# P <- "padj"
# XLABELS <- "mean of normalized counts"
# YLABELS <- "log2 fold change"
# PLABELS <- "adjusted p-value"
# 
# INDIR <- "../results/2022March/DEG/tables"
# OUTDIR <- "../results/2022March/DEG/plots"
# annotFilePath <- "../data/gene_symbol_2.tsv"
# COMPR <- c("2_1", "4_1", "5_1", "4_2", "5_2", "5_4")
# 
# n <- sapply(COMPR, plotMAs, inDir=INDIR, outDir=OUTDIR,
#             x=X, y=Y, p=P, 
#             xlab=XLABELS, ylab=YLABELS, yThres=LFC_THRES, pThres=PVAL_THRES)
# 
# n <- sapply(COMPR, plotVolcano, 
#             inDir=INDIR, outDir=OUTDIR,
#             anotFn=annotFilePath, 
#             y=Y, p=P, xlab=YLABELS, ylab=PLABELS, yThres=LFC_THRES, pThres=PVAL_THRES)
