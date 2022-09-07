#!/usr/bin/env Rscript

'Do enrichment analysis using raw count data
Usage:
    remove_batch_effect.R [--method=<method> --group] <input> <metadata> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --method=<method>  The name of p-value adjusting method [default: combatseq]
    --group  Consider group column in the metadata

Arguments:
    input  input csv file name
    metadata  metadata file name
    output  output folder
' -> doc


library(biomaRt)
library(glue)
library(sva)
library(RUVSeq)
library(limma)
library(readxl)
library(here)
suppressMessages(library(docopt))

cleanY <- function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

args <- docopt(doc)
df <- read.csv(file.path(args$input), row.names = 1)
df <- df[rowSums(df) != 0,]

meta <- read.csv(args$metadata, row.names = 1)
batch <- meta$batch
# batch <- paste(batch, meta$condition)
if (args$group){
  group <- meta$group
  mod <- model.matrix(~group)
  mod0 <- cbind(mod[,1])
}


if (tolower(args$method) == "combatseq") {
  if (args$group){
    adjusted <- sva::ComBat_seq(as.matrix(df), batch=batch, group=group)
  } else {
    adjusted <- sva::ComBat_seq(as.matrix(df), batch=batch)
  }
} else if (tolower(args$method) == "limma") {
  limma_fit <- limma::removeBatchEffect(x=log(as.matrix(df) + 1), batch = factor(batch))
  adjusted <- exp(limma_fit) - 1
} else if (method == "ruv") { 
  # not finished
  controls <- (group == 1)
  ruvfit <- RUVg(df, cIdx= controls, k=1)
  batch_ruv_cp <- ruvfit$W
} else {
  stop("Invalid method")
}

write.csv(adjusted, args$output)

# sva
# controls <- (group == 1)  # for supervised methods
# svafit <- sva::svaseq(as.matrix(df), mod = mod, mod0=mod0, n.sv=2)
# batch_sva_us <- svafit$sv

# svafit <- sva::svaseq(as.matrix(df), mod = mod, mod0=mod0, controls = controls, n.sv=2)
# batch_sva_s <- svafit$sv
# 
# plot(batch_sva_us, col=factor(batch), cex=1,
#      xlab="SV1", ylab="SV2")

# combat
# sva::ComBat()