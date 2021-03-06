#!/usr/bin/env Rscript

'Do enrichment analysis using raw count data
Usage:
    enrichment_main.R (--enrich | --gsea) [--dePval=<dePval> --pval=<pval> --qval=<qval> --adjMethod=<adjMethod> --list] <input> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --enrich  Perform enrichment analysis
    --gsea  Perform GSEA
    --list  If specified, the input data is treated as a list of genes, otherwise it will be treated as a DEG table
    --dePval=<dePval>  The p-value used to select DEGs [default: 0.05]
    --pval=<pval>  The p-value cutoff used to run the analysis [default: 0.05]
    --qval=<qval>  The q-value cutoff used to run the analysis [default: 0.2]
    --adjMethod=<adjMethod> The name of p-value adjusting method [default: BH]

Arguments:
    input  inpur folder storing DEG results
    output  output folder
' -> doc


# should use valid gene id here
suppressMessages(library(docopt))
suppressMessages(library(stringr))


source("./enrichment_funcs.r")
source("./gsea_funcs.r")

args <- docopt(doc)

main <- function(args){
  if (args$list){
    print("Using gene list files")
    comparisons <- list.files(args$input, full.names	= F)
    print(comparisons)
  } else {
    print("Using DEG table files")
    comparisons <- gsub("^DE_|.tsv", "", grep("^DE_", list.files(path = args$input), value = TRUE))
    sapply(comparisons, split_DE, inDir=args$input, pvalthres=args$dePval)
  }
  if (args$enrich) {
    for (db in c("go", "kegg", "reactome", "wp", "msigdb")) {
      do_enrichment(comparisons, 
                    input_dir=args$input, 
                    output_dir=args$output,
                    is_list=args$list,
                    adjust_method=args$adjMethod,
                    p_cutoff=as.numeric(args$pval),
                    q_cutoff=as.numeric(args$qval),
                    type=db)
    }
  } else if (args$gsea) {
    for (db in c("go", "kegg", "reactome", "wp", "msigdb")) {
      do_GSEA(comparisons, 
              input_dir=args$input, 
              output_dir=args$output,
              adjust_method=args$adjMethod,
              p_cutoff=as.numeric(args$pval),
              type=db)
    }
  }
}

main(args)
