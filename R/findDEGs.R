#!/usr/bin/env Rscript

'Do DE analysis using raw count data
Usage:
    findDEGs.R [--batch=<batch> --pval=<pval> --logfc=<logfc> --adjMethod=<adjMethod>] <input> <output> <metadata>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --batch=<batch>  To consider the batch information in the metadata if this gets the column name of batches [default: None]
    --pval=<pval>  Width of the output [default: 0.05]
    --logfc=<logfc>  The absolute log2 fold change cutoff used to determine the DEGs [default: 1]
    --adjMethod=<adjMethod> The name of p-value adjusting method [default: BH]

Arguments:
    input  input folder storing raw count data
    output  output folder
    metadata  A file containing samples information
' -> doc


suppressMessages(library(docopt))
suppressMessages(library(DESeq2))
suppressMessages(library(glue))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(grid))
suppressMessages(library(stringr))
suppressMessages(library(tools))
suppressMessages(library(readxl))
suppressMessages(library(scales))


GROUP_COLORS <- c("#3367c7", "#f08300", "#17b055", "#5e17b0", "#9c1f1f")

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


plotPCAc <- function(pcaTable, condName, fileName){
  jpeg(fileName, 
       width = 9, height = 9, units = "in",
       res = 600, quality = 100)
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
  print(p)
  dev.off()
}

plotMAs <- function(comp, inDir, outDir, x, y, p, xlab, ylab, yThres=log2(1.5), pThres=0.1, logScale=T){
  t <- read.table(file.path(inDir, glue("DE_{comp[1]}_{comp[2]}.tsv")), sep='\t')
  t <- t[complete.cases(t),]
  t$de <- 0
  t$de[(t[, y] > yThres)&(t[, p] < pThres)] <- 1
  t$de[(t[, y] < (-yThres))&(t[, p] < pThres)] <- (-1)
  t$de <- as.factor(t$de)
  deLevels <- levels(t$de)
  color_values <- c()
  color_labels <- c()
  if (-1 %in% deLevels){
    color_values <- c(color_values, "#0466c8")
    color_labels <- c(color_labels, "Down regulated")
  }
  color_values <- c(color_values, "#999999")
  color_labels <- c(color_labels, "Not DE genes")
  if ( 1 %in% deLevels){
    color_values <- c(color_values, "#e63946")
    color_labels <- c(color_labels, "Up regulated")
  }
  
  degname <- strsplit(inDir, "/")
  degname <- degname[[1]][length(degname[[1]])]
  
  p <- ggplot(t, aes_string(x=x, y=y, color="de")) + 
    geom_point(size=1) +
    scale_color_manual(name = "DE genes", values=color_values, labels=color_labels) + 
    geom_hline(yintercept=0, size=1.5, alpha=0.4) +
    labs(x = xlab, y = ylab, title = paste(strsplit(comp, "_")[[1]], collapse = " vs "))
  if (logScale){
    p + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))
  }
  ggsave(file.path(outDir, glue("{degname}_{comp[1]}_{comp[2]}.jpeg")), width = 9, height = 9, units = "in", device="jpeg")
}

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

plotVolcano <- function(comp, inDir, outDir, anotFn, y, p, xlab, ylab, yThres=log2(1.5), pThres=0.1, pSig=1e-15){
  print(file.path(inDir, glue("DE_{comp[1]}_{comp[2]}.tsv")))
  t <- read.table(file.path(inDir, glue("DE_{comp[1]}_{comp[2]}.tsv")), sep='\t', header=1)
  if (dim(str_match(row.names(t), "(.*)\\..*")) >= 2){
    t["Gene"] <- str_match(row.names(t), "(.*)\\..*")[,2]  # strip the version id
  }
  annot <- read.table(anotFn, sep='\t', header=1)
  t <- t[complete.cases(t),]
  t$de <- 0
  t$de[(t[, y] > yThres)&(t[, p] < pThres)] <- 1
  t$de[(t[, y] < (-yThres))&(t[, p] < pThres)] <- (-1)
  t$de <- as.factor(t$de)
  
  dat <- data.frame(DE = t$de)
  dat <- dat %>% 
    group_by(DE) %>%
    summarise(no_rows = length(DE))
  n_upreg <- dat[dat$DE==1,]$no_rows
  n_downreg <- dat[dat$DE==-1,]$no_rows
  
  t <- merge(t, annot, by="Gene")
  deLevels <- levels(t$de)
  color_values <- c()
  color_labels <- c()
  if (-1 %in% deLevels){
    color_values <- c(color_values, "#0466c8")
    color_labels <- c(color_labels, "Down regulated")
  }
  color_values <- c(color_values, "#999999")
  color_labels <- c(color_labels, "Not DE genes")
  if ( 1 %in% deLevels){
    color_values <- c(color_values, "#e63946")
    color_labels <- c(color_labels, "Up regulated")
  }
  
  degname <- strsplit(inDir, "/")
  degname <- degname[[1]][length(degname[[1]])]
  grobUp <- grobTree(textGrob(glue("{n_upreg}"), x=0.9, y=0.9, gp=gpar(col="#e63946", fontsize=18)))
  grobDown <- grobTree(textGrob(glue("{n_downreg}"), x=0.1, y=0.9, gp=gpar(col="#0466c8", fontsize=18)))
  p <- ggplot(t, aes_string(x=y, y=p, color="de")) + 
    geom_point(size=1) +
    scale_color_manual(name = "DE genes", values=color_values, labels=color_labels) + 
    # geom_vline(xintercept=0, size=1.5, alpha=0.4) +
    geom_vline(xintercept=yThres, size=1, alpha=0.3) +
    geom_vline(xintercept=-yThres, size=1, alpha=0.3) +
    geom_hline(yintercept=pThres, size=1, alpha=0.3) +
    labs(x = xlab, y = ylab, title = paste(strsplit(comp, "_")[[1]], collapse = " vs ")) +
    scale_y_continuous(trans=reverselog_trans(10)) + 
    annotation_custom(grobUp) +
    annotation_custom(grobDown) + 
    geom_text_repel(
      data = t[(t[, p] < pSig)&(complete.cases(t))&(!duplicated(t)), ],
      aes(label = Symbol),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      show.legend = F
    )
  ggsave(file.path(outDir, glue("volcano_{degname}_{comp[1]}_{comp[2]}.jpeg")), width = 9, height = 9, units = "in", device="jpeg")
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

pca_factors <- c("group")
design_formula <- "~"

metadata$group <- factor(metadata$group)

if (args$batch != "None"){
  if (args$batch %in% colnames(metadata)){
    metadata[, c(args$batch)] <- (metadata %>% select(args$batch) %>% mutate_all(factor))
    pca_factors <- append(pca_factors, args$batch)
    design_formula <- paste(design_formula, args$batch, sep = "")
  } else {
    stop("The batch columns is not contained in the metadata")
  }
}

design_formula <- paste(design_formula, ifelse(nchar(design_formula) <= 1, "group", "+group"), sep = "")

comparisons <- get_comparison(levels(metadata$group))
# outliers <- args$skip  # not implemented yet
# metadata <- metadata[-match(outliers, metadata$sample),]
# colnames(countsData) <- str_match(colnames(countsData), "(.*)A")[,2]
# countsData <- countsData[, -match(outliers, colnames(countsData))]

print("creating DESeqDataSet..")
dds <- DESeqDataSetFromMatrix(counts_data[, metadata$sample], 
                              colData = metadata, 
                              design = as.formula(design_formula))

vsd <- vst(dds)

for (f in pca_factors){
  print(glue("plotting PCA plot for factor: {f}"))
  vsdResult <- plotPCA(vsd, intgroup=c(f), returnData=T)
  plotPCAc(vsdResult, condName = f, file.path(PLOT_RESULT_DIR, glue("PCA_counts_vst_{f}.jpeg")))
}

print("Using thresholds: ")
print(glue("p-val < {as.numeric(args$pval)} and abs(logFC) > {as.numeric(args$logfc)}"))

dds <- DESeq(dds)
for (i in 1:length(comparisons))
{
  print(glue("{comparisons[[i]][1]} vs {comparisons[[i]][2]}"))
  res <- results(dds, contrast = c("group", comparisons[[i]]), pAdjustMethod = args$adjMethod)
  notna <- res[!is.na(res$padj), ]
  notna <- notna[(abs(notna$log2FoldChange) > as.numeric(args$logfc)) & (notna$padj < as.numeric(args$pval)), ]
  write.table(res, file.path(TABLE_RESULT_DIR, glue("DE_{comparisons[[i]][1]}_{comparisons[[i]][2]}.tsv")), sep = "\t", row.names = TRUE)
  write.table(data.frame(DGE =row.names(notna)), file.path(TABLE_RESULT_DIR, glue("DEG_list_{comparisons[[i]][1]}_{comparisons[[i]][2]}.tsv")), sep = "\t", row.names = TRUE)
}
# 
X <- "baseMean"
Y <- "log2FoldChange"
P <- "padj"
XLABELS <- "mean of normalized counts"
YLABELS <- "log2 fold change"
PLABELS <- "adjusted p-value"
annotFilePath <- "./gene_data/gene_symbols.tsv"

n <- sapply(comparisons, plotMAs, inDir=TABLE_RESULT_DIR, outDir=PLOT_RESULT_DIR,
            x=X, y=Y, p=P,
            xlab=XLABELS, ylab=YLABELS, yThres=as.numeric(args$logfc), pThres=as.numeric(args$pval))

n <- sapply(comparisons, plotVolcano,
            inDir=TABLE_RESULT_DIR, outDir=PLOT_RESULT_DIR,
            anotFn=annotFilePath,
            y=Y, p=P, xlab=YLABELS, ylab=PLABELS, yThres=as.numeric(args$logfc), pThres=as.numeric(args$pval))
