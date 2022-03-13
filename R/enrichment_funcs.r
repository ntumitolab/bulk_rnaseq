library(glue)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(ggplot2)
library(DOSE)

GOsubontologies <- c("ALL", "BP", "CC", "MF")
expressionDirection <- c("down", "up",  "all")
usedDB <- "org.Hs.eg.db"
org_idx <- "hsa"
org_name <- "human"
if (org_name == "human"){
  orgDb = org.Hs.eg.db
}

split_DE <- function (comparison, inDir, log2thres = log2(2), pvalthres = 0.1, method = "DESeq", isadjust = TRUE){
  df <- read.table(file.path(inDir, glue("DE_{comparison}.tsv")), sep='\t')
  df <- df[complete.cases(df), ]
  upRegGenes <- df[(df$log2FoldChange >= log2thres)&(df$padj < pvalthres),]
  upRegGenes <- upRegGenes[complete.cases(upRegGenes), ]
  downRegGenes <- df[(df$log2FoldChange <= -log2thres)&(df$padj < pvalthres),]
  downRegGenes <- downRegGenes[complete.cases(downRegGenes), ]
  allRegGenes <- df[(abs(df$log2FoldChange) >= log2thres)&(df$padj < pvalthres),]
  allRegGenes <- allRegGenes[complete.cases(allRegGenes), ]
  write.table(upRegGenes, file.path(inDir, glue("DEG_list_{comparison}_up.tsv")), sep='\t', row.names = T)
  write.table(downRegGenes, file.path(inDir, glue("DEG_list_{comparison}_down.tsv")), sep='\t', row.names = T)
  write.table(allRegGenes, file.path(inDir, glue("DEG_list_{comparison}_all.tsv")), sep='\t', row.names = T)
}


findGenes <- function(genes){
  A <- tryCatch({
    return(bitr(genes, fromType="ENSEMBL", toType=c("ENSEMBL", "ENTREZID"), OrgDb=usedDB))
  }, error = function(e) {
    print(e)
    return(NULL)
  }, finally={
  })
  return(A)
}

makeDf <- function(fileName, printInfo = T){
  DE_data <- read.table(fileName, sep='\t')
  if (org_idx == "mmu"){
    genes <- substr(row.names(DE_data), 1, 18)
  } else if (org_idx == "hsa") {
    genes <- substr(row.names(DE_data), 1, 15)
  }
  if(printInfo){
    print(glue("gene length: {length(genes)}"))
  }
  DEG.df <- findGenes(genes)
  if (is.null(DEG.df)) {
    if(printInfo){
      print("No deg.df")
    }
    return(NULL)
  } else {
    if(printInfo){
      print(DEG.df$ENTREZID)
      print(glue("ENTREZID length: {length(DEG.df$ENTREZID)}"))
    }
    return(DEG.df)
  }
}

drawEnrichMaps <- function(enrichment, geneList, fileName, outDirFigs, orgDb, ont=NULL){
  enrichx <- setReadable(enrichment, orgDb, keyType = "ENTREZID")
  jpeg(file.path(outDirFigs, ifelse(is.null(ont), glue("cnet_{fileName}.jpeg"), glue("cnet_{fileName}_{ont}.jpeg"))),
       width = 9, height = 9, units = "in",
       res = 300, quality = 100)
  plt <- enrichplot::cnetplot(enrichx, foldChange=geneList)
  print(plt)
  dev.off()
  
  jpeg(file.path(outDirFigs, ifelse(is.null(ont), glue("heatplot_{fileName}.jpeg"), glue("heatplot_{fileName}_{ont}.jpeg"))),
       width = 20, height = 8, units = "in",
       res = 300, quality = 100)
  plt <- enrichplot::heatplot(enrichx, showCategory=50, foldChange=geneList)
  print(plt)
  dev.off()
  
  sim <- pairwise_termsim(enrichment)
  if (dim(sim)[1] > 1){
    print(sim)
    print(dim(sim))
    jpeg(file.path(outDirFigs, ifelse(is.null(ont), glue("emap_{fileName}.jpeg"), glue("emap_{fileName}_{ont}.jpeg"))),
         width = 8, height = 8, units = "in",
         res = 300, quality = 100)
    plt <- enrichplot::emapplot(sim, layout="kk")
    print(plt)
    dev.off()
  }
  
}

saveGOtable <- function(subontology, DEG.df, usedDB, pCutoff, qCutoff, pAdjustMethod, inDir, outDir, outDirFigs, comparison, tableFnPrefix){
  ego <- enrichGO(gene = DEG.df$ENTREZID, 
                  OrgDb = usedDB,
                  pAdjustMethod = pAdjustMethod, 
                  ont = subontology,
                  pvalueCutoff = pCutoff,
                  qvalueCutoff = qCutoff,
                  readable = TRUE)
  if (!is.null(ego) && dim(ego)[1] != 0)
  {
    write.table(ego, file=file.path(outDir, glue("table_{tableFnPrefix}_{subontology}.tsv")), sep="\t")
    if (is.null(outDirFigs) == F){
      drawEnrichMaps(ego, makeGSEAData(inDir, glue("DE_{comparison}.tsv"), useEntrez = T), 
                     orgDb=orgDb, 
                     ont = subontology,
                     fileName=tableFnPrefix, 
                     outDirFigs=outDirFigs)
    }
  }
}

getGOEnrichTables <- function(comparison, inDir, outDir, outDirFigs=NULL, pCutoff=0.05, qCutoff=0.2, org=org_name){
  for (ed in expressionDirection){
    print(glue("converting {comparison}_{ed} ..."))
    DEG.df <- makeDf(file.path(inDir, glue("DEG_list_{comparison}_{ed}.tsv")))
    if (is.null(DEG.df)) {
      next
    }
    sapply(GOsubontologies, 
           saveGOtable, 
           DEG.df=DEG.df, 
           usedDB=usedDB, 
           pCutoff=pCutoff, 
           pAdjustMethod = "BH",
           qCutoff=qCutoff, 
           inDir=inDir,
           outDir=outDir, 
           outDirFigs=outDirFigs,
           comparison=comparison,
           tableFnPrefix=glue("{comparison}_{ed}"))
  }
}

getKEGGEnrichTables <- function(comparison, inDir, outDir, outDirFigs=NULL, pCutoff=0.05, qCutoff=0.2, org=org_name){
  for (ed in expressionDirection){
    print(glue("converting {comparison}_{ed} ..."))
    DEG.df <- makeDf(file.path(inDir, glue("DEG_list_{comparison}_{ed}.tsv")))
    if (is.null(DEG.df)) {
      next
    }
    ekegg <- enrichKEGG(DEG.df$ENTREZID, 
                        organism = org_idx, 
                        keyType = "kegg", 
                        pvalueCutoff = pCutoff,
                        pAdjustMethod = "BH",
                        minGSSize = 10, 
                        maxGSSize = 500,
                        qvalueCutoff = qCutoff, 
                        use_internal_data = F)
    if (!is.null(ekegg) && dim(ekegg)[1] != 0){
      write.table(ekegg, file.path(outDir, glue("/table_{comparison}_{ed}.tsv")), sep="\t")
      if (is.null(outDirFigs) == F){
        drawEnrichMaps(ekegg, makeGSEAData(inDir, glue("DE_{comparison}.tsv"), useEntrez = T), orgDb=orgDb, 
                       fileName=glue("{comparison}_{ed}"), outDirFigs=outDirFigs)
      }
    }
    else {
      print("ekegg is null")
    }
  }
}

getReactomeEnrichTables <- function(comparison, inDir, outDir, outDirFigs=NULL, pCutoff=0.05, qCutoff=0.2, org=org_name){
  for (ed in expressionDirection){
    print(glue("converting {comparison}_{ed} ..."))
    DEG.df <- makeDf(file.path(inDir, glue("DEG_list_{comparison}_{ed}.tsv")))
    if (is.null(DEG.df)) {
      next
    }
    erea <- enrichPathway(DEG.df$ENTREZID, 
                          organism = org, 
                          pvalueCutoff = pCutoff, 
                          pAdjustMethod = "BH",
                          minGSSize = 10, 
                          maxGSSize = 500, 
                          qvalueCutoff = qCutoff,
                          readable = TRUE)
    if (!is.null(erea) && dim(erea)[1] != 0){
      write.table(erea, file.path(outDir, glue("/table_{comparison}_{ed}.tsv")), sep="\t")
      if (is.null(outDirFigs) == F){
        drawEnrichMaps(erea, makeGSEAData(inDir, glue("DE_{comparison}.tsv"), useEntrez = T), orgDb=orgDb, 
                       fileName=glue("{comparison}_{ed}"), outDirFigs=outDirFigs)
      }
    }
    else {
      print("erea is null")
    }
  }
}

makeGSEAData <- function(inDir, fileName, useEntrez=F){
  df <- read.table(file.path(inDir, fileName), sep='\t')
  if (org_idx == "mmu"){
    genes <- substr(row.names(df), 1, 18)
  } else if (org_idx == "hsa") {
    genes <- substr(row.names(df), 1, 15)
  }
  gseInp <- df$log2FoldChange
  names(gseInp) <- genes
  if (useEntrez){
    mapped <- clusterProfiler::bitr(genes, "ENSEMBL", "ENTREZID", org.Hs.eg.db)
    gseInp <- gseInp[mapped$ENSEMBL]
    names(gseInp) <- mapped$ENTREZID
  }
  gseInp <- sort(gseInp, decreasing = T)
  return(gseInp)
}

plotSingleGSEA <- function(id, gseResult, outDirFigs, fileName, ont=NULL){
  jpeg(file.path(outDirFigs, ifelse(is.null(ont), glue("GSEA_{fileName}_{id}.jpeg"), glue("GSEA_{fileName}_{ont}_{id}.jpeg"))),
       width = 12, height = 9, units = "in",
       res = 300, quality = 100)
  plt <- enrichplot::gseaplot2(gseResult, geneSetID = id, title = gseResult$Description[id])
  print(plt)
  dev.off()
  return(NULL)
}

savegseGOresult <- function(ont, gseData, usedDB, pCutoff, outDirTable, 
                            outDirFigs, fileName, nTop, saveplot=T){
  gseGO.res <- gseGO(gseData, 
                     usedDB, 
                     keyType = "ENSEMBL", 
                     ont=ont, 
                     pvalueCutoff=pCutoff,
                     pAdjustMethod = "BH")
  if (!is.null(gseGO.res) && dim(gseGO.res)[1] != 0){
    write.table(gseGO.res, file.path(outDirTable, glue("table_{fileName}_{ont}.tsv")), sep="\t")
    if (saveplot){
      jpeg(file.path(outDirFigs, glue("ridgeplot_{fileName}_{ont}.jpeg")),
           width = 9, height = 9, units = "in",
           res = 300, quality = 100)
      p <- enrichplot::ridgeplot(gseGO.res) + ggplot2::labs(x = "enrichment distribution")
      print(p)
      dev.off()
      ns <- sapply(seq(min(dim(gseGO.res)[1], nTop)), 
                   plotSingleGSEA, gseResult=gseGO.res, outDirFigs=outDirFigs, fileName=fileName, ont=ont)
    }
    
  }
  else {
    print("gseGO.res is null")
  }
  return(gseGO.res)
}

savegseKEGGresult <- function(gseData, organism, pCutoff, outDirTable, outDirFigs, fileName, nTop, saveplot=T){
  gseKEGG.res <- gseKEGG(gseData, 
                         organism=organism, 
                         keyType = "ncbi-geneid", 
                         pvalueCutoff=pCutoff,
                         pAdjustMethod = "BH")
  if (!is.null(gseKEGG.res) && dim(gseKEGG.res)[1] != 0){
    write.table(gseKEGG.res, file.path(outDirTable, glue("table_{fileName}.tsv")), sep="\t")
    if (saveplot){
      print(file.path(outDirFigs, glue("ridge_{fileName}.jpeg")))
      jpeg(file.path(outDirFigs, glue("ridge_{fileName}.jpeg")),
           width = 9, height = 9, units = "in",
           res = 300, quality = 100)
      p <- enrichplot::ridgeplot(gseKEGG.res) + ggplot2::labs(x = "enrichment distribution")
      print(p)
      dev.off()
      ns <- sapply(seq(min(dim(gseKEGG.res)[1], nTop)), 
                   plotSingleGSEA, gseResult=gseKEGG.res, outDirFigs=outDirFigs, fileName=fileName)
    }
    
  }
  else {
    print("gseKEGG.res is null")
  }
  return(gseKEGG.res)
}

savegseReactomeresult <- function(gseData, organism, pCutoff, outDirTable, outDirFigs, fileName, nTop, saveplot=T){
  gseReac.res <- gsePathway(gseData, 
                           organism=organism, 
                           pvalueCutoff=pCutoff,
                           pAdjustMethod = "BH")
  if (!is.null(gseReac.res) && dim(gseReac.res)[1] != 0){
    write.table(gseReac.res, file.path(outDirTable, glue("table_{fileName}.tsv")), sep="\t")
    if (saveplot){
      jpeg(file.path(outDirFigs, glue("ridge_{fileName}.jpeg")),
           width = 9, height = 9, units = "in",
           res = 300, quality = 100)
      p <- enrichplot::ridgeplot(gseReac.res) + ggplot2::labs(x = "enrichment distribution")
      print(p)
      dev.off()
      ns <- sapply(seq(min(dim(gseReac.res)[1], nTop)), 
                   plotSingleGSEA, gseResult=gseReac.res, outDirFigs=outDirFigs, fileName=fileName)
    }
    
  }
  else {
    print("gseReac.res is null")
  }
  return(gseReac.res)
}

plotgseGO <- function(comparison, inDir, outDirTable, outDirFigs, usedDB=orgDb, 
                      pCutoff=0.05, nTop=30, saveplot=T){
  gseData <- makeGSEAData(inDir, glue("DE_{comparison}.tsv"))
  sapply(GOsubontologies[2:4], 
         savegseGOresult, 
         gseData=gseData, 
         usedDB=usedDB, 
         pCutoff=pCutoff,
         outDirTable=outDirTable, 
         outDirFigs=outDirFigs, 
         fileName=comparison,
         nTop=nTop,
         saveplot=saveplot)
}

plotgseKEGG <- function(comparison, inDir, outDirTable, outDirFigs, usedDB=orgDb, 
                        pCutoff=0.05, nTop=30, saveplot=T){
  gseData <- makeGSEAData(inDir, glue("DE_{comparison}.tsv"), useEntrez=T)
  savegseKEGGresult(
    gseData=gseData, 
    organism = org_idx,
    pCutoff=pCutoff,
    outDirTable=outDirTable, 
    outDirFigs=outDirFigs, 
    fileName=comparison,
    nTop=nTop,
    saveplot=saveplot
  )
  return(gseData)
}

plotgseReactome <- function(comparison, inDir, outDirTable, outDirFigs, usedDB=orgDb, 
                            pCutoff=0.05, nTop=30, saveplot=T){
  gseData <- makeGSEAData(inDir, glue("DE_{comparison}.tsv"), useEntrez=T)
  savegseReactomeresult(
          gseData=gseData, 
          organism = org_name,
          pCutoff=pCutoff,
          outDirTable=outDirTable, 
          outDirFigs=outDirFigs, 
          fileName=comparison,
          nTop=nTop,
          saveplot=saveplot
  )
}