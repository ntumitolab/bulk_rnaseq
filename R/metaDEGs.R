library(here)
library(metaRNASeq)
library(glue)


save_result <- function(p_list, prefix, updir, downdir, gene_names){
  
  result_df <- as.data.frame(list("raw_pval"=p_list$rawpval, 
                                  "adj_pval"=p_list$adjpval), 
                             row.names = gene_names)
  up_deg <- subset(result_df, updir & result_df$adj_pval < 0.05)
  down_deg <- subset(result_df, downdir & result_df$adj_pval < 0.05)
  write.csv(result_df, glue("../../DOX/results/1_degs/metadeg/{prefix}_pval.csv"))
  write.table(row.names(up_deg), glue("../../DOX/results/1_degs/metadeg/{prefix}_up_deg.csv"), 
              row.names = F, col.names = F)
  write.table(row.names(down_deg), glue("../../DOX/results/1_degs/metadeg/{prefix}_down_deg.csv"), 
              row.names = F, col.names = F)
}

fcval <- read.csv(here("../../DOX/results/1_degs/metadeg/MC_fc.csv"), row.names = 1)
updir <- rowSums(fcval) == dim(fcval)[2]
downdir <- rowSums(fcval) == -dim(fcval)[2]

rawpval <- read.csv(here("../../DOX/results/1_degs/metadeg/MC_pval.csv"), row.names = 1)
rawpval[is.na(rawpval)] <- 1


n_samples <- c(36, 36, 36, 10, 16, 8, 8, 8, 16, 9, 9, 6, 6) # HT c(9, 30) CC c(12, 8, 8, 12, 12, 8, 12, 6)
fisher_p <- fishercomb(rawpval, BHth = 0.05)
invnorm_p <- metaRNASeq::invnorm(rawpval, n_samples, BHth = 0.05)

save_result(fisher_p, "MC_fisher", updir, downdir, gene_names=row.names(rawpval))
save_result(invnorm_p, "MC_invnorm", updir, downdir, gene_names=row.names(rawpval))








hist(rawpval[,2], breaks=100, col="grey", main="Fisher method",xlab = "Raw p-values (meta-analysis)")

hist(fishcomb$rawpval, breaks=100, col="grey", main="Fisher method",xlab = "Raw p-values (meta-analysis)")