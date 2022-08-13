rm(list=ls())

library(WGCNA)
library(flashClust)
library(DESeq2)
library(readxl)
library(tidyverse)
library(magrittr)
library(here)
library(pheatmap)

print_res <- 450

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}


plot_jpeg <- function(p, file_name, width=12, height=8){
  px_ft <- 120
  jpeg(file_name,height=height*px_ft,width=width*px_ft, quality = 100, res = print_res)
  p
  dev.off()
}

plot_power_est <- function(RpowerTable, powers1, cex1=0.7){
  cex1 = 0.7
  op <- par(mfrow = c(1,2), mar=c(3,3,3,3))
  plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], 
       xlab = "soft threshold (power)", 
       ylab = "scale free topology model fit, signes R^2", type = "n")
  text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], 
       labels = powers1, 
       cex = cex1, col = "red")
  abline(h = 0.95, col = "red")
  plot(RpowerTable[,1], 
       RpowerTable[,5], xlab = "soft threshold (power)", ylab = "mean connectivity", type = "n")
  text(RpowerTable[,1], RpowerTable[,5], labels = powers1, cex = cex1, col = "red")
}

cal_best_power <- function(data, 
                           network_type="signed", 
                           draw_power=F,
                           file_name=NaN){
  
  if (dim(data)[1] > dim(data)[2]){
    warning("Data should be a sample by gene df")
  }
  
  powers1 <- c(seq(1, 10, by=1), seq(12, 20, by=2))
  sft <- pickSoftThreshold(data, powerVector = powers1, networkType = network_type)
  RpowerTable <- sft[[2]]
  if (draw_power) {
    if (!is.nan(file_name)) {
      plot_jpeg(plot_power_est(RpowerTable, powers1), 
                file_name, width=24, height = 12)
    } else {
      plot_power_est(RpowerTable, powers1)
    }
  }
  gc()
  sft
  
}

draw_data_dist <- function(data){
  mdata <- as_tibble(t(data))
  mdata <- mdata %>%
    tidyr::pivot_longer(
      .,
      col = all_of(names(mdata))
    )
  
  p <- mdata %>%
    ggplot(., aes(x = name, y = value)) +             # x = treatment, y = RNA Seq count
    geom_violin() +                                   # violin plot, show distribution
    geom_point(alpha = 0.2) +                         # scatter plot
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)          # Rotate treatment text
    ) +
    labs(x = "Treatment Groups", y = "RNA Seq Counts")
  p
}

get_good_data <- function(data){
  gsg = goodSamplesGenes(data, verbose = 3)
  
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", 
                       paste(names(data)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", 
                       paste(rownames(data)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    data = data[gsg$goodSamples, gsg$goodGenes]
  }
  
  data
}

select_median <- function(data, num=5000) {
  
  keepGenesExpr <- rank(-colMedians(data))<=num
  data_sg <- data[,keepGenesExpr]
  data_sg
}

select_median_2 <- function(data) {
  gene_median <- colMedians(data)
  data[, gene_median > mean(gene_median)]
}


calc_cont_fac_cor <- function(ME, sample_orders, factor_df) {
  moduleTraitCor <- cor(ME[sample_orders,], factor_df, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, dim(ME)[1])
  colnames(moduleTraitCor) <- paste(colnames(moduleTraitCor), "cor", sep = "_")
  colnames(moduleTraitPvalue) <- paste(colnames(moduleTraitPvalue), "pval", sep = "_")
  merge(moduleTraitCor, moduleTraitPvalue, by=0)
}


get_TOM <- function(data, power, type="signed") {
  adjacency_mat <- adjacency(data, power=power, type=type)
  diag(adjacency_mat) <- 0
  TOMsimilarity(adjacency_mat, TOMType=type)
}


get_module_split_choices <- function(gene_cluster, dissTOM, file_name="Module_choices.pdf", cut_method="hybrid"){
  mColorh <- NULL
  for (ds in 0:4){
    if (cut_method == "hybrid"){
    tree <- cutreeHybrid(dendro = gene_cluster, pamStage=FALSE,
                        minClusterSize = (30-3*ds), cutHeight = 0.99,
                        deepSplit = ds, distM = dissTOM)
    } else if (cut_method == "dynamic") {
    tree <- cutreeDynamic(dendro = gene_cluster, distM = dissTOM, deepSplit = ds, 
                          pamRespectsDendro = FALSE, minClusterSize = (30-3*ds))
    }
    mColorh<- cbind(mColorh,labels2colors(tree$labels));
  }
  
  plot_jpeg(plotDendroAndColors(gene_cluster, mColorh, paste("dpSplt =",0:4), main = "",dendroLabels=FALSE),
            file_name, height=10, width=25);
  
  mColorh
}

get_ME <- function(data, colors){
  PC1 <- moduleEigengenes(data, colors=colors)
  MEs <- PC1$eigengenes
  distPC1 <- 1-abs(cor(MEs,use="p"))
  distPC1 <- ifelse(is.na(distPC1), 0, distPC1)
  pcTree <- hclust(as.dist(distPC1),method="a")
  MDS <- cmdscale(as.dist(distPC1),2)
  PC_colors <- names(table(colors))
  list(PC1=PC1, MEs=MEs, distPC1=distPC1, pcTree=pcTree, MDS=MDS, PC_colors=PC_colors)
}

calc_MM <- function(data, ME, colors) {
  geneModuleMembership <- signedKME(data, ME)
  colnames(geneModuleMembership) <- paste("PC", colors,".cor",sep="")
  MMPvalue <- corPvalueStudent(as.matrix(geneModuleMembership),dim(data)[[1]])
  colnames(MMPvalue) <- paste("PC",colors,".pval",sep="")
  genes <- colnames(data)
  kMEtable = cbind(genes, colors)
  for (i in 1:length(colors))
    kMEtable = cbind(kMEtable, geneModuleMembership[,i], MMPvalue[,i])
  colnames(kMEtable)=c("gene","Module",sort(c(colnames(geneModuleMembership),
                                              colnames(MMPvalue)))) 
  kMEtable
}

calc_GS <- function(data, y){
  GS <- as.numeric(WGCNA::cor(y, data, use="p"))
  p.Standard <- corPvalueFisher(GS, nSamples =length(y) )
  p.Standard2 <- p.Standard
  p.Standard2[is.na(p.Standard)] <- 1
  q.Standard <- qvalue(p.Standard2)$qvalues
  data.frame(gene=colnames(data),PearsonCorrelation=GS, p.Standard, q.Standard)
}



save_network <- function(MEs, data, TOM, file_dir, threshold=0.5){
  for (i in 1:length(MEs)){
    colors <- c(substring(names(MEs)[i], 3));
    genes <- colnames(data)
    inModule <- is.finite(match(colors, colors))
    modGenes <- genes[inModule]
    modTOM <- TOM[inModule,inModule]
    dimnames(modTOM) <- list(modGenes,modGenes)
    cyt <- exportNetworkToCytoscape(modTOM,
                                    edgeFile = here::here(file_dir, paste("edges_", 
                                                                          paste(colors, collapse="_"), 
                                                                          ".txt", sep="")),
                                    nodeFile = here::here(file_dir, paste("nodes_", 
                                                                          paste(colors, collapse="_"), 
                                                                          ".txt", sep="")),
                                    weighted = TRUE, 
                                    threshold = threshold, 
                                    nodeNames = modGenes, 
                                    nodeAttr = colors[inModule])
  }
}


options(stringsAsFactors = F)
enableWGCNAThreads()

# for CC
output_path <- here::here("../../DOX/results/3_wgcna/CC/")
sampleInfo <- read.csv(here("../../DOX/data/merged/wgcna_1/meta.csv"))
countsData <- read.table(here("../../DOX/data/merged/wgcna_1/count_br.csv"), 
                         sep=',', header = TRUE, row.names = 1)
dds <- DESeqDataSetFromMatrix(countsData[, sampleInfo$sample], 
                              colData = sampleInfo, 
                              design = ~dosage)


# for CO
output_path <- here::here("../../DOX/results/3_wgcna/CO/")
sampleInfo <- readxl::read_excel("../../DOX/data/meta/S-HECA10.xlsx")
sampleInfo <- as.data.frame(sampleInfo)
countsData <- read.table(here("../../DOX/data/studies/S-HECA10/counts.csv"), 
                         sep=',', header = TRUE, row.names = 1)

dds <- DESeqDataSetFromMatrix(countsData[, sampleInfo$sample], 
                              colData = sampleInfo, 
                              design = ~dosage + time + time:dosage)

# for MC
output_path <- here::here("../../DOX/results/3_wgcna/MC/")
sampleInfo <- readxl::read_excel("../../DOX/data/merged/wgcna_3/meta.xlsx")
sampleInfo <- as.data.frame(sampleInfo)
countsData <- read.table(here("../../DOX/data/merged/wgcna_3/count_br.csv"), 
                         sep=',', header = TRUE, row.names = 1)

dds <- DESeqDataSetFromMatrix(countsData[, sampleInfo$sample], 
                              colData = sampleInfo, 
                              design = ~dosage + time + time:dosage)


# for HT
output_path <- here::here("../../DOX/results/3_wgcna/HT/")
sampleInfo <- read.csv("../../DOX/data/merged/wgcna_4/meta.csv")

countsData <- read.table(here("../../DOX/data/merged/wgcna_4/count_br.csv"), 
                         sep=',', header = TRUE, row.names = 1)

dds <- DESeqDataSetFromMatrix(countsData[, sampleInfo$sample], 
                              colData = sampleInfo, 
                              design = ~dosage)



sampleInfo$dosage <- factor(sampleInfo$dosage)
sampleInfo$time <- factor(sampleInfo$time)

annot <- sampleInfo[, c(3, 6, 8)]
row.names(annot) <- sampleInfo[, 2]
dose <- ifelse(annot$dosage == "dox0.6uM", 0.6, ifelse(annot$dosage == "dox0.2uM", 0.2, 0))
dose <- ifelse(annot$dosage == "dox", 1, 0)
dose <- annot$dosage

time <- as.integer(unlist(strsplit(as.vector(annot$time), split = "hrs")))
colnames(annot) <- c("dosage", "time (hrs)", "cell type")

# data processing
vsd <- vst(dds)
vsdData <- log2(vsd@assays@data@listData[[1]] + 1)
vsdData <- t(quantile_normalisation(vsdData))

# filtering
vsdData <- get_good_data(vsdData)


# plot sample dendrogram
sampleTree = hclust(dist(vsdData), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")


# filter genes by their exp median
vsdData_sg <- select_median(vsdData, 10000)

# get soft threshold
sft <- cal_best_power(vsdData_sg, draw_power = T, file_name = here(output_path, "power_est.jpeg"))

write.csv(sft$fitIndices, here(output_path, "sft_search.csv"))
power <- sft$powerEstimate
# power <- 12 # if no good power is found

# get TOM and gene cluster
TOM <- get_TOM(vsdData_sg, power=power, type="signed")
dissTOM <- 1-TOM
geneTree <- flashClust(as.dist(dissTOM), method="average") 
plot_jpeg(plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity",
          labels=FALSE,hang=0.04), 
          width = 48, height = 36,
          here::here(output_path, "TOM_diss_cluster.jpeg"))

# choose the best module here
m_choices <- get_module_split_choices(geneTree, dissTOM, 
                                      file_name = here::here(output_path, "module_choices.jpeg"))
unmergedColors <- m_choices[, 2]


# after choosing the module, merge those too close to each other
ME_list <- get_ME(vsdData_sg, unmergedColors)


plot(ME_list$pcTree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col = "red")

merge <- mergeCloseModules(vsdData_sg, unmergedColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge$colors

ME_list_mg <- get_ME(vsdData_sg, mergedColors)

write.csv(ME_list_mg$MEs, here::here(output_path, "MEs.csv"))

plot_jpeg(plotDendroAndColors(geneTree, cbind(unmergedColors, mergedColors),
                              c("Unmerged", "Merged"),
                              dendroLabels = FALSE, hang = 0.03,
                              addGuide = TRUE, guideHang = 0.05),
          here(output_path, "merged_color.jpeg"), width = 36, height = 24)

# save the network (threshold = 90% quantile)
th <- quantile(TOM, 0.99)[["99%"]]
gc()
save_network(ME_list_mg$MEs, vsdData_sg, TOM, 
             here::here(output_path, "modules"), 
             threshold = th)

# calculate MM


kMEtable <- calc_MM(vsdData_sg, ME_list_mg$MEs, ME_list_mg$PC_colors)

write.csv(kMEtable, here::here(output_path, "kMEtable.csv"), row.names=FALSE)

# calculate GS


jpeg(here::here(output_path, "heatmap.jpeg"), width=960*2, height=960*2, quality=100, res=print_res)
ph <- pheatmap(ME_list_mg$MEs,
               cluster_col=T,
               cluster_row=T,
               show_rownames=F,
               show_colnames=T,
               fontsize=4,
               annotation_row =annot[c(1, 3)])
grid::grid.newpage()
grid::grid.draw(ph$gtable)
dev.off()

gs_dose <- calc_GS(vsdData_sg, dose)
gs_time <- calc_GS(vsdData_sg, time)
write.csv(gs_dose, here::here(output_path, "GS_dose.csv"), row.names=FALSE)
write.csv(gs_time, here::here(output_path, "GS_time.csv"), row.names=FALSE)

# if the factor is continuous
fac_cor <- calc_cont_fac_cor(ME_list_mg$MEs, 
                             row.names(annot), 
                             as.data.frame(list("dose"=dose, "time"=time)))
write.csv(fac_cor, here::here(output_path, "factor_cors.csv"))

con_trait <- standardScreeningNumericTrait(ME_list_mg$MEs, dose)
write.csv(con_trait, here::here(output_path, "GS_dox.csv"))

# if the factor is binary
bin_trait <- standardScreeningBinaryTrait(ME_list_mg$MEs, sampleInfo$dosage == "dox")
write.csv(bin_trait, here::here(output_path, "GS_dox.csv"))


# (optional) get final module plot
jpeg("Final_modules.png",height=8,width=12, quality = 100)
plotDendroAndColors(geneTree, mergedColors, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors")
dev.off() 


# (optional) some other plots
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)

plot(ME_list_mg$MDS, col= ME_list_mg$PC_colors, main="MDS plot", cex=2, pch=19)
ordergenes = geneTree$order
plotMat(scale(log(vsdData_sg)[,ordergenes]) ,clabels = mergedColors[ordergenes], rlabels=
           row.names(vsdData_sg), ccols=mergedColors[ordergenes])
for (which.module in names(table(mergedColors))){
  ME = ME_1A[, paste("ME",which.module, sep="")]
  barplot(ME, col=which.module, main="", cex.main=2,
          ylab="eigengene expression",xlab="array sample")
}
