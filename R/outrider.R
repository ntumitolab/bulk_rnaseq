library(OUTRIDER)
library(biomaRt)
library(matrixStats)
library(here)

ctsFile <- here("../Research/code/NTUH/data/mg_count_july_br_combatseq.csv")
ctsTable <- read.csv(ctsFile)

if(interactive()){
  
  mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  g = getBM(c("hgnc_symbol", "ensembl_gene_id"), 
            "ensembl_gene_id", 
            row.names(ctsTable), mart)
  names(g) = c("hgnc_symbol", "gene")
  g <- merge(x = ctsTable, y = g, by="gene")
  write.table(g, glue("f:/ntuh/hg38/Annot_merged_count.tsv"), sep='\t')
}
g <- g[complete.cases(g[, c("hgnc_symbol")]), ]
g <- g[(!duplicated(g[, c("hgnc_symbol")])), ]
row.names(g) <- g[, c("hgnc_symbol")]
g <- g[, -c(1, 1068)]

write.csv(g, here("../Research/code/NTUH/data/mg_count_july_for_outrider.csv"))
g <- read.csv("../Research/code/NTUH/data/mg_count_july_for_outrider.csv", row.names = 1)
g <- g[rank(apply(g, 1, var)) > dim(g)[1] - 3000,]

ods <- OutriderDataSet(countData=g)
ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
ods <- OUTRIDER(ods)
res <- results(ods, padjCutoff=0.1, zScoreCutoff=2)
write.csv(res, "../../user/Research/code/NTUH/results/2022July/outrider_top3000_genes.csv")
write.csv(aberrant(ods, by="sample"), "../../user/Research/code/NTUH/results/2022July/outrider_aberrant_sammples.csv")
write.csv(aberrant(ods, by="gene"), "../../user/Research/code/NTUH/results/2022July/outrider_aberrant_genes.csv")

tail(sort(aberrant(ods, by="sample")))
plotAberrantPerSample(ods, padjCutoff=0.05)