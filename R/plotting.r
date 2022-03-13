library(ggplot2)
library(ggrepel)
library(scales)
library(glue)
library(dplyr)
library(grid) # grob

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

plotMAs <- function(comp, inDir, outDir, x, y, p, xlab, ylab, yThres=log2(1.5), pThres=0.1, logScale=T){
  t <- read.table(file.path(inDir, glue("DE_{comp}.tsv")), sep='\t')
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
  ggsave(file.path(outDir, glue("{degname}_{comp}.jpeg")), width = 9, height = 9, units = "in", device="jpeg")
}

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

plotVolcano <- function(comp, inDir, outDir, anotFn, y, p, xlab, ylab, yThres=log2(1.5), pThres=0.1, pSig=1e-15){
  print(file.path(inDir, glue("DE_{comp}.tsv")))
  t <- read.table(file.path(inDir, glue("DE_{comp}.tsv")), sep='\t', header=1)
  t["Gene"] <- row.names(t)
  
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
  ggsave(file.path(outDir, glue("volcano_{degname}_{comp}.jpeg")), width = 9, height = 9, units = "in", device="jpeg")
}
