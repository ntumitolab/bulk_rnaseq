library(glue)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(DOSE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./utils.R")


goAnalyzer <- setRefClass("goAnalyzer", contains = "enrichAnalyzer",
                          methods = list(
                            get_enrich_obj = function(subclass, DEG.df) {
                              clusterProfiler::enrichGO(gene = DEG.df$ENTREZID, 
                                                        OrgDb = usedDB,
                                                        pAdjustMethod = adjust_method, 
                                                        ont = subclass,
                                                        pvalueCutoff = p_cutoff,
                                                        qvalueCutoff = q_cutoff,
                                                        readable = TRUE)
                            },
                            iter = function(comparison, ed, DEG.df) {
                              sapply(GOsubontologies, 
                                     .self$get_table, 
                                     DEG.df=DEG.df,
                                     comparison=comparison,
                                     tbl_filename=glue("{comparison}_{ed}"))
                              
                            }
                          ))


keggAnalyzer <- setRefClass("keggAnalyzer", contains = "enrichAnalyzer",
                          methods = list(
                            get_enrich_obj = function(subontology, DEG.df) {
                              clusterProfiler::enrichKEGG(gene = DEG.df$ENTREZID, 
                                                          organism = org_short, 
                                                          pAdjustMethod = adjust_method, 
                                                          keyType = "kegg", 
                                                          minGSSize = 10, 
                                                          maxGSSize = 500,
                                                          pvalueCutoff = p_cutoff,
                                                          qvalueCutoff = q_cutoff,
                                                          use_internal_data = F)
                            
                            },
                            iter = function(comparison, ed, DEG.df) {
                              .self$get_table(comparison=comparison,
                                              DEG.df=DEG.df,
                                              tbl_filename=glue("{comparison}_{ed}"),
                                              subclass=NULL)
                            }
                          ))

reactomeAnalyzer <- setRefClass("reactomeAnalyzer", contains = "enrichAnalyzer",
                            methods = list(
                              get_enrich_obj = function(subontology, DEG.df) {
                                ReactomePA::enrichPathway(DEG.df$ENTREZID, 
                                                          organism = organism, 
                                                          pvalueCutoff = p_cutoff, 
                                                          pAdjustMethod = adjust_method,
                                                          minGSSize = 10, 
                                                          maxGSSize = 500, 
                                                          qvalueCutoff = q_cutoff,
                                                          readable = TRUE)
                                
                              }
                            ))

msigdbAnalyzer <- setRefClass("msigdbAnalyzer", contains = "enrichAnalyzer",
                              methods = list(
                                get_enrich_obj = function(subclass, DEG.df) {
                                  m_t2g <- msigdbr(species = species, category = subclass) %>% 
                                    dplyr::select(gs_name, entrez_gene)
                                  
                                  clusterProfiler::enricher(gene = DEG.df$ENTREZID, 
                                                            pAdjustMethod = adjust_method, 
                                                            pvalueCutoff = p_cutoff,
                                                            qvalueCutoff = q_cutoff,
                                                            TERM2GENE=m_t2g)
                                },
                                iter = function(comparison, ed, DEG.df) {
                                  sapply(msigdbr_gs, 
                                         .self$get_table, 
                                         DEG.df=DEG.df,
                                         comparison=comparison,
                                         tbl_filename=glue("{comparison}_{ed}"))
                                  
                                }
                              ))

wpAnalyzer <- setRefClass("wpAnalyzer", contains = "enrichAnalyzer",
                          methods = list(
                            get_enrich_obj = function(subclass, DEG.df) {
                              clusterProfiler::enrichWP(DEG.df$ENTREZID, 
                                                        organism = species, 
                                                        pvalueCutoff = p_cutoff, 
                                                        pAdjustMethod = adjust_method,
                                                        qvalueCutoff = q_cutoff,
                                                        minGSSize = 10, 
                                                        maxGSSize = 500)
                            }
                          ))


davidAnalyzer <- setRefClass("davidAnalyzer", contains = "enrichAnalyzer",
                             methods = list(
                               get_enrich_obj = function(subclass, DEG.df) {
                                 clusterProfiler::enrichDAVID(DEG.df$ENTREZID, 
                                                           organism = species, 
                                                           pvalueCutoff = p_cutoff, 
                                                           pAdjustMethod = adjust_method,
                                                           qvalueCutoff = q_cutoff,
                                                           minGSSize = 10, 
                                                           maxGSSize = 500)
                               }
                             ))


enrich_analyzer <- list("go" = goAnalyzer, 
                        "kegg" = keggAnalyzer, 
                        "reactome" = reactomeAnalyzer,
                        "msigdb" = msigdbAnalyzer,
                        "wp" = wpAnalyzer
                        )

##' A general interface for doing enrichment analysis
do_enrichment <- function(comparisons, 
                          input_dir, 
                          output_dir, 
                          type="go",
                          organism="human", 
                          adjust_method="BH",
                          p_cutoff=0.05, 
                          q_cutoff=0.2) {
  if (is.null(enrich_analyzer[[type]])){
    stop("The type is not valid")
  }
  
  output_dir.tbl = file.path(output_dir, type, "tables")
  output_dir.fig = file.path(output_dir, type, "plots")
  
  ifelse(!dir.exists(output_dir.tbl), 
         dir.create(output_dir.tbl, recursive = T), FALSE)
  ifelse(!dir.exists(output_dir.fig), 
         dir.create(output_dir.fig, recursive = T), FALSE)
  
  analyzer <- enrich_analyzer[[type]](name = type,
                                      organism = organism,
                                      input_dir = input_dir,
                                      output_dir = output_dir.tbl,
                                      output_dir.figs = output_dir.fig,
                                      q_cutoff = q_cutoff,
                                      p_cutoff = p_cutoff,
                                      adjust_method = adjust_method)
  sapply(comparisons, analyzer$get_all_tables)
  
}

