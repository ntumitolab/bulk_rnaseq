
# should use valid gene id here
library(stringr)

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./enrichment_funcs.r")
comparisons <- c("2_1", "4_1", "5_1", "4_2", "5_2", "5_4")
DEG_IN_DIR <- "../results/DEG/tables"
ERC_OUT_DIR <- "../results/enrichments"
GSEA_OUT_DIR <- "../results/GSEA"

# data splitting
sapply(comparisons, split_DE, inDir=DEG_IN_DIR)

# enrichment analysis
sapply(comparisons, getGOEnrichTables, 
       inDir=DEG_IN_DIR, 
       outDir=file.path(ERC_OUT_DIR, "go/tables"),
       outDirFigs=file.path(ERC_OUT_DIR, "go/plots"))

sapply(comparisons, getKEGGEnrichTables, 
       inDir=DEG_IN_DIR, 
       outDir=file.path(ERC_OUT_DIR, "kegg/tables"), 
       outDirFigs=file.path(ERC_OUT_DIR, "kegg/plots"))

sapply(comparisons, getReactomeEnrichTables, 
       inDir=DEG_IN_DIR, 
       outDir=file.path(ERC_OUT_DIR, "reactome/tables"),
       outDirFigs=file.path(ERC_OUT_DIR, "reactome/plots"))

sapply(comparisons, getDOEnrichTables, 
       inDir=DEG_IN_DIR, 
       outDir=file.path(ERC_OUT_DIR, "do/tables"),
       outDirFigs=file.path(ERC_OUT_DIR, "do/plots"))

# GSEA analysis
gseaGOs <- sapply(comparisons, plotgseGO, 
                  inDir=DEG_IN_DIR, 
                  outDirTable=file.path(GSEA_OUT_DIR, "go/tables"), 
                  outDirFigs=file.path(GSEA_OUT_DIR, "go/plots"))

gseaKEGGs <- sapply(comparisons, plotgseKEGG, 
                    inDir=DEG_IN_DIR, 
                    outDirTable=file.path(GSEA_OUT_DIR, "kegg/tables"), 
                    outDirFigs=file.path(GSEA_OUT_DIR, "kegg/plots"))

gseaReactomes <- sapply(comparisons, plotgseReactome, 
                        inDir=DEG_IN_DIR, 
                        outDirTable=file.path(GSEA_OUT_DIR, "reactome/tables"), 
                        outDirFigs=file.path(GSEA_OUT_DIR, "reactome/plots"))

gseaDOs <- sapply(comparisons, plotgseDO, 
                  inDir=DEG_IN_DIR, 
                  outDirTable=file.path(GSEA_OUT_DIR, "do/tables"), 
                  outDirFigs=file.path(GSEA_OUT_DIR, "do/plots"))