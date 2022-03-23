
# should use valid gene id here
library(stringr)

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./enrichment_funcs.r")
source("./gsea_funcs.r")
comparisons <- c("2_1", "4_1", "5_1", "4_2", "5_2", "5_4")
DEG_IN_DIR <- "../example/tables"
ERC_OUT_DIR <- "../results/enrichments"
GSEA_OUT_DIR <- "../results/GSEA"

# data splitting
sapply(comparisons, split_DE, inDir=DEG_IN_DIR)

# enrichment analysis
for (db in c("go", "kegg", "reactome", "wp", "msigdb")) {
  do_enrichment(comparisons, 
                input_dir=DEG_IN_DIR, 
                output_dir=ERC_OUT_DIR,
                type=db)
}

do_enrichment(comparisons, 
              input_dir=DEG_IN_DIR, 
              output_dir=ERC_OUT_DIR,
              type="msigdb")

# GSEA analysis
for (db in c("go", "kegg", "reactome", "wp", "msigdb")) {
  do_GSEA(comparisons, 
          input_dir=DEG_IN_DIR, 
          output_dir=GSEA_OUT_DIR,
          type=db)
}

do_GSEA(comparisons, 
              input_dir=DEG_IN_DIR, 
              output_dir=ERC_OUT_DIR,
              type="msigdb")