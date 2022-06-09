library(glue)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(DOSE)
library(msigdbr)


source("./utils.R")

goGSEA <- setRefClass("goGSEA", contains = "GSEAAnalyzer",
                          methods = list(
                            get_gsea_obj = function(subclass, gsea_data) {
                              clusterProfiler::gseGO(gsea_data, 
                                                     OrgDb = orgDB, 
                                                     keyType = "ENSEMBL", 
                                                     ont=subclass, 
                                                     pvalueCutoff=p_cutoff,
                                                     pAdjustMethod = adjust_method)
                            },
                            iter = function(comparison, gsea_data) {
                              sapply(GOsubontologies[2:4], 
                                     comparison=comparison,
                                     .self$get_all_tables, 
                                     gsea_data=gsea_data)
                              
                            }
                          ))

keggGSEA <- setRefClass("keggGSEA", contains = "GSEAAnalyzer",
                      methods = list(
                        get_gsea_obj = function(subclass, gsea_data) {
                          clusterProfiler::gseKEGG(gsea_data, 
                                                   organism = org_short, 
                                                   keyType = "ncbi-geneid", 
                                                   pvalueCutoff=p_cutoff,
                                                   pAdjustMethod = adjust_method)
                        }
                      ))


reactomeGSEA <- setRefClass("reactomeGSEA", contains = "GSEAAnalyzer",
                        methods = list(
                          get_gsea_obj = function(subclass, gsea_data) {
                            ReactomePA::gsePathway(gsea_data, 
                                                   organism=organism, 
                                                   pvalueCutoff=p_cutoff,
                                                   pAdjustMethod = adjust_method)
                          }
                        ))

wpGSEA <- setRefClass("wpGSEA", contains = "GSEAAnalyzer",
                        methods = list(
                          get_gsea_obj = function(subclass, gsea_data) {
                            clusterProfiler::gseWP(gsea_data, 
                                                    organism = species, 
                                                    pvalueCutoff=p_cutoff,
                                                    pAdjustMethod = adjust_method)
                          }
                        ))

msigdbGSEA <- setRefClass("msigdbGSEA", contains = "GSEAAnalyzer",
                      methods = list(
                        get_gsea_obj = function(subclass, gsea_data) {
                          m_t2g <- msigdbr(species = species, category = subclass) %>% 
                            dplyr::select(gs_name, entrez_gene)
                          
                          clusterProfiler::GSEA(gsea_data, 
                                                TERM2GENE=m_t2g,
                                                pvalueCutoff=p_cutoff,
                                                pAdjustMethod = adjust_method)
                        },
                        iter = function(comparison, gsea_data) {
                          sapply(msigdbr_gs, 
                                 comparison=comparison,
                                 .self$get_all_tables, 
                                 gsea_data=gsea_data)
                          
                        }
                      ))


gsea_analyzer <- list("go" = goGSEA, 
                      "kegg" = keggGSEA, 
                      "reactome" = reactomeGSEA,
                      "wp" = wpGSEA,
                      "msigdb" = msigdbGSEA
)


##' A general interface for doing GSEA
do_GSEA <- function(comparisons, 
                    input_dir, 
                    output_dir, 
                    type="go",
                    organism="human", 
                    nTop=30,
                    adjust_method="BH",
                    p_cutoff=0.05) {
  if (is.null(gsea_analyzer[[type]])){
    stop("The type is not valid")
  }
  
  output_dir.tbl = file.path(output_dir, type, "tables")
  output_dir.fig = file.path(output_dir, type, "plots")
  
  ifelse(!dir.exists(output_dir.tbl), 
         dir.create(output_dir.tbl, recursive = T), FALSE)
  ifelse(!dir.exists(output_dir.fig), 
         dir.create(output_dir.fig, recursive = T), FALSE)
  
  
  analyzer <- gsea_analyzer[[type]](name = type,
                                    organism = organism,
                                    input_dir = input_dir,
                                    output_dir = output_dir.tbl,
                                    output_dir.figs = output_dir.fig,
                                    p_cutoff = p_cutoff,
                                    nTop=nTop,
                                    adjust_method = adjust_method)
  sapply(comparisons, analyzer$do_gsea)
  
}

