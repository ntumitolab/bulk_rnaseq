library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)


GOsubontologies <- c("ALL", "BP", "CC", "MF")
msigdbr_gs <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
exp_dirs <- c("down", "up",  "all")



split_DE <- function (comparison, inDir, log2thres = log2(2), pvalthres = 0.05, method = "DESeq", isadjust = TRUE){
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


findGenes <- function(genes, usedDB){
  A <- tryCatch({
    return(bitr(genes, fromType="ENSEMBL", toType=c("ENSEMBL", "ENTREZID"), OrgDb=usedDB))
  }, error = function(e) {
    print(e)
    return(NULL)
  }, finally={
  })
  return(A)
}

makeDf <- function(fileName, usedDB, org_short="hsa", printInfo = T){
  DE_data <- read.table(fileName, sep='\t')
  if (org_short == "mmu"){
    genes <- substr(row.names(DE_data), 1, 18)
  } else if (org_short == "hsa") {
    genes <- substr(row.names(DE_data), 1, 15)
  }
  if(printInfo){
    print(glue("gene length: {length(genes)}"))
  }
  DEG.df <- findGenes(genes, usedDB)
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
  enrichx <- DOSE::setReadable(enrichment, orgDb, keyType = "ENTREZID")
  jpeg(file.path(outDirFigs, ifelse(is.null(ont), glue("cnet_{fileName}.jpeg"), glue("cnet_{fileName}_{ont}.jpeg"))),
       width = 1920, height = 1920, units = "px",
       res = 150, quality = 100)
  plt <- enrichplot::cnetplot(enrichx, showCategory=50, foldChange=geneList)
  print(plt)
  dev.off()
  
  jpeg(file.path(outDirFigs, ifelse(is.null(ont), glue("heatplot_{fileName}.jpeg"), glue("heatplot_{fileName}_{ont}.jpeg"))),
       width = 1920, height = 960, units = "px",
       res = 150, quality = 100)
  plt <- enrichplot::heatplot(enrichx, showCategory=50, foldChange=geneList)
  print(plt)
  dev.off()
  
  sim <- pairwise_termsim(enrichment)
  if (dim(sim)[1] > 1){
    print(sim)
    print(dim(sim))
    jpeg(file.path(outDirFigs, ifelse(is.null(ont), 
                                      glue("emap_{fileName}.jpeg"), 
                                      glue("emap_{fileName}_{ont}.jpeg"))),
         width = 1920, height = 1920, units = "px",
         res = 150, quality = 100)
    plt <- enrichplot::emapplot(sim, showCategory=50, layout="kk")
    print(plt)
    dev.off()
  }
  
}


### GSEA ###
makeGSEAData <- function(inDir, fileName, org_short="hsa", useEntrez=F){
  df <- read.table(file.path(inDir, fileName), sep='\t')
  if (org_short == "mmu"){
    genes <- substr(row.names(df), 1, 18)
  } else if (org_short == "hsa") {
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
  jpeg(file.path(outDirFigs, ifelse(is.null(ont), 
                                    glue("GSEA_{fileName}_{id}.jpeg"), 
                                    glue("GSEA_{fileName}_{ont}_{id}.jpeg"))),
       width = 12, height = 9, units = "in",
       res = 300, quality = 100)
  plt <- enrichplot::gseaplot2(gseResult, geneSetID = id, title = gseResult$Description[id])
  print(plt)
  dev.off()
  return(NULL)
}



enrichAnalyzer <- setRefClass("enrichAnalyzer", fields = list(name = "character",
                                                              organism = "character",
                                                              input_dir = "character",
                                                              output_dir = "character",
                                                              output_dir.figs = "character",
                                                              q_cutoff = "numeric",
                                                              p_cutoff = "numeric",
                                                              adjust_method = "character",
                                                              usedDB = "character",
                                                              org_short = "character",
                                                              species = "character", 
                                                              orgDB = "ANY"),
                              
                              methods = list(
                                get_enrich_obj = function(...) {
                                  NULL
                                },
                                iter = function(comparison, ed, DEG.df, ...) {
                                  .self$get_table(comparison=comparison,
                                                  DEG.df=DEG.df,
                                                  tbl_filename=glue("{comparison}_{ed}"),
                                                  subclass=NULL)
                                },
                                convert_org_name = function(){
                                  if (organism == "human") {
                                    usedDB <<- "org.Hs.eg.db"
                                    org_short <<- "hsa"
                                    species <<- "Homo sapiens"
                                    orgDB <<- org.Hs.eg.db
                                  } else if (organism == "mouse") {
                                    usedDB <<- "org.Mm.eg.db"
                                    org_short <<- "mmu"
                                    species <<- "Mus musculus"
                                  } else {
                                    stop("The organism is not valid for now")
                                  }
                                },
                                get_table = function(subclass, comparison, tbl_filename, ...) {
                                  
                                  if (is.null(subclass)){
                                    e_obj <- get_enrich_obj(...)
                                  } else {
                                    e_obj <- get_enrich_obj(subclass=subclass,...)
                                  }
                                  
                                  if (!is.null(e_obj) && dim(e_obj)[1] != 0) {
                                    write.table(e_obj, file=file.path(output_dir, 
                                                                      ifelse(!is.null(subclass),
                                                                             glue("table_{tbl_filename}_{subclass}.tsv"),
                                                                             glue("table_{tbl_filename}.tsv"))), 
                                                sep="\t")
                                    if (!is.null(output_dir.figs)){
                                      drawEnrichMaps(e_obj, 
                                                     makeGSEAData(input_dir, 
                                                                  org_short=org_short,
                                                                  glue("DE_{comparison}.tsv"), 
                                                                  useEntrez = T), 
                                                     orgDb=orgDB, 
                                                     ont = subclass,
                                                     fileName=tbl_filename, 
                                                     outDirFigs=output_dir.figs)
                                    }
                                  }
                                }, 
                                
                                get_all_tables = function(comparison) {
                                  convert_org_name()
                                  for (ed in exp_dirs) {
                                    print(glue("processing {comparison}_{ed} ..."))
                                    DEG.df <- makeDf(file.path(input_dir, 
                                                                glue("DEG_list_{comparison}_{ed}.tsv")),
                                                     usedDB = usedDB,
                                                     org_short = org_short)
                                    if (is.null(DEG.df)) {
                                      next
                                    }
                                    iter(comparison=comparison, ed=ed, DEG.df=DEG.df)
                                  }
                                }))

GSEAAnalyzer <- setRefClass("GSEAAnalyzer", fields = list(name = "character",
                                                          organism = "character",
                                                          input_dir = "character",
                                                          output_dir = "character",
                                                          output_dir.figs = "character",
                                                          q_cutoff = "numeric",
                                                          p_cutoff = "numeric",
                                                          nTop = "numeric",
                                                          adjust_method = "character",
                                                          usedDB = "character",
                                                          org_short = "character",
                                                          species = "character", 
                                                          orgDB = "ANY"),
                            methods = list(
                              get_gsea_obj = function(...) {
                                NULL
                              },
                              iter = function(comparison, gsea_data) {
                                .self$get_all_tables(subclass=NULL,
                                                     comparison=comparison,
                                                     gsea_data=gsea_data)
                                
                              },
                              convert_org_name = function(){
                                if (organism == "human") {
                                  usedDB <<- "org.Hs.eg.db"
                                  org_short <<- "hsa"
                                  species <<- "Homo sapiens"
                                  orgDB <<- org.Hs.eg.db
                                } else if (organism == "mouse") {
                                  usedDB <<- "org.Mm.eg.db"
                                  org_short <<- "mmu"
                                  species <<- "Mus musculus"
                                } else {
                                  stop("The organism is not valid for now")
                                }
                              },
                              get_all_tables = function(subclass, comparison, ...) {
                                if (is.null(subclass)){
                                  gsea <- get_gsea_obj(...)
                                } else {
                                  gsea <- get_gsea_obj(subclass=subclass, ...)
                                }
                                
                                
                                if (!is.null(gsea) && dim(gsea)[1] != 0){
                                  write.table(gsea, file.path(output_dir, 
                                                              ifelse(!is.null(subclass),
                                                                     glue("table_{comparison}_{subclass}.tsv"),
                                                                     glue("table_{comparison}.tsv"))), 
                                              sep="\t")
                                  if (!is.null(output_dir.figs)){
                                    jpeg(file.path(output_dir.figs, 
                                                   ifelse(!is.null(subclass),
                                                          glue("ridgeplot_{comparison}_{subclass}.jpeg"),
                                                          glue("ridgeplot_{comparison}.jpeg"))),
                                         width = 9, height = 9, units = "in",
                                         res = 300, quality = 100)
                                    p <- enrichplot::ridgeplot(gsea) + ggplot2::labs(x = "enrichment distribution")
                                    print(p)
                                    dev.off()
                                    ns <- sapply(seq(min(dim(gsea)[1], nTop)), 
                                                 plotSingleGSEA, 
                                                 gseResult=gsea, 
                                                 outDirFigs=output_dir.figs, 
                                                 fileName=comparison, 
                                                 ont=subclass)
                                  }
                                  
                                }
                                else {
                                  print("GSEA result is null")
                                }
                                return(gsea)
                              },
                              do_gsea = function(comparison) {
                                convert_org_name()
                                useEntrez <- name %in% c("msigdb", "wp", "kegg", "reactome")
                                
                                gsea_data <- makeGSEAData(input_dir, 
                                                          glue("DE_{comparison}.tsv"),
                                                          useEntrez = useEntrez)
                                
                                iter(comparison, gsea_data)
                              }
                            ))
