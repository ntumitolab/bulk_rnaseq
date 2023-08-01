library(here)
library(metaRNASeq)
library(rmeta)
library(glue)

IO_DIR <- "../../DOX/results/1_degs/metadeg_compare/"
prefix = "IPSC"
save_result <- function(p_list, prefix, updir, downdir, gene_names){
  
  result_df <- as.data.frame(list("raw_pval"=p_list$rawpval, 
                                  "adj_pval"=p_list$adjpval), 
                             row.names = gene_names)
  up_deg <- subset(result_df, updir & result_df$adj_pval < 0.05)
  down_deg <- subset(result_df, downdir & result_df$adj_pval < 0.05)
  write.csv(result_df, file.path(IO_DIR, glue("{prefix}_pval.csv")))
  write.table(row.names(up_deg), 
              file.path(IO_DIR, glue("{prefix}_metadeg"), glue("{prefix}_up_deg.csv")), 
              row.names = F, col.names = F)
  write.table(row.names(down_deg), 
              file.path(IO_DIR, glue("{prefix}_down_deg.csv")), 
              row.names = F, col.names = F)
}


fcval <- read.csv(here(file.path(IO_DIR, glue("{prefix}_fc.csv"))), row.names = 1)
updir <- rowSums(fcval) == dim(fcval)[2]
downdir <- rowSums(fcval) == -dim(fcval)[2]

rawpval <- read.csv(here(file.path(IO_DIR, glue("{prefix}_pval.csv"))), row.names = 1)
rawpval[is.na(rawpval)] <- 1

fisher_p <- fishercomb(rawpval, BHth = 0.05)
save_result(fisher_p, glue("{prefix}_fisher"), updir, downdir, gene_names=row.names(rawpval))

n_samples <- c(36, 36, 36, 10, 16, 16, 9, 9, 6, 6) # HT c(9, 30, 61) CC c(12, 8, 8, 12, 12, 8, 12, 6)
invnorm_p <- metaRNASeq::invnorm(rawpval, n_samples, BHth = 0.05)
save_result(invnorm_p, glue("{prefix}_invnorm"), updir, downdir, gene_names=row.names(rawpval))

# plotting
hist(rawpval[,2], breaks=100, col="grey", main="Fisher method",xlab = "Raw p-values (meta-analysis)")
hist(fishcomb$rawpval, breaks=100, col="grey", main="Fisher method",xlab = "Raw p-values (meta-analysis)")

# calculate meta logFC using rmeta
fc_se <- read.csv(here(file.path(IO_DIR, glue("{prefix}_fc_se.csv"))), row.names = 1)

get_new_fc_df <- function(fc_se){
  fc_part <- fc_se[,seq(1, dim(fc_se)[2] / 2)]
  se_part <- fc_se[,seq(dim(fc_se)[2] / 2 + 1, dim(fc_se)[2])]
  combined_fc <- c()
  for (i in seq(1, dim(fc_part)[1])) {
    fc <- as.numeric(fc_part[i,])
    fc <- fc[!is.na(fc)]
    se <- as.numeric(se_part[i,])
    se <- se[!is.na(se)]
    if (length(se) == 0){
      combined_fc <- c(combined_fc, NA)
    } else {
      metafc <- meta.summaries(fc, se, "fixed", logscale = T)$summary
      combined_fc <- c(combined_fc, metafc)
    }
  }
  new_fc <- data.frame(list(logFC=combined_fc))
  row.names(new_fc) <- row.names(fc_se)
  new_fc
}
new_fc <- get_new_fc_df(fc_se)
write.csv(new_fc, here(file.path(IO_DIR, glue("{prefix}_new_fc.csv"))))

