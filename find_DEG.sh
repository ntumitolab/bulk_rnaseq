Rscript findDEGs.R ../data/AC16_1 ../results/AC16_1/DEG ../data/metadata/AC16_1.xlsx
Rscript enrichment_main.r --enrich ../results/AC16_1/DEG/tables ../results/AC16_1/enrichments
Rscript enrichment_main.r --gsea ../results/GSE106297/DEG/tables ../results/GSE106297/GSEA

Rscript enrichment_main.r --gsea ../results/GSE157282/DEG/tables ../results/GSE157282/GSEA

Rscript enrichment_main.r --gsea ../results/AC16_1/DEG/tables ../results/AC16_1/GSEA