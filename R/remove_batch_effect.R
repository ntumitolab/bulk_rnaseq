install.packages('ggfortify')
BiocManager::install("RUVSeq")

library(biomaRt)
library(glue)
library(sva)
library(RUVSeq)
library(limma)
library(readxl)


library(ggfortify)

cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

root_dir <- file.path("../../NTUH/data")
df <- read.csv(file.path(root_dir, "mg_count_may.csv"), row.names = 1)
df <- df[rowSums(df) != 0,]
meta <- as.data.frame(readxl::read_excel(file.path(root_dir, "all_info/labels_may_mg.xlsx")))
batch <- meta$batch
group <- meta$group
mod <- model.matrix(~group)
mod0 <- cbind(mod[,1])
controls <- (group == 1)  # for supervised methods

# PCA
pca_res <- prcomp(t(df), scale. = TRUE)
autoplot(pca_res, data=meta, colour = 'group')

# limma
limma_fit <- limma::removeBatchEffect(x=log(as.matrix(df) + 1), batch = factor(batch))
pca_res <- prcomp(t(limma_fit), scale. = TRUE)
autoplot(pca_res, data=meta, colour = 'group')
write.csv(exp(limma_fit) - 1, file.path(root_dir, "mg_count_may_br_limma.csv"))

# ruv
ruvfit <- RUVg(df, cIdx= controls, k=1)
batch_ruv_cp <- ruvfit$W

# sva
svafit <- sva::svaseq(as.matrix(df), mod = mod, mod0=mod0, n.sv=2)
batch_sva_us <- svafit$sv

# svafit <- sva::svaseq(as.matrix(df), mod = mod, mod0=mod0, controls = controls, n.sv=2)
# batch_sva_s <- svafit$sv

plot(batch_sva_us, col=factor(batch), cex=1,
     xlab="SV1", ylab="SV2")

# combat
# sva::ComBat()

adjusted <- sva::ComBat_seq(as.matrix(df), batch=batch)
pca_res <- prcomp(t(adjusted), scale. = TRUE)
autoplot(pca_res, data=meta, colour = 'group')

write.csv(adjusted, file.path(root_dir, "mg_count_may_br_combatseq.csv"))
