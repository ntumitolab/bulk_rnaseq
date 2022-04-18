# Bulk RNAseq

## Usage

### Without Docker
0. Make sure you have R and Python installed
1. Navigate to this directory
2. Install the requirements.r in the R folder
```bash
Rscript .R/requirements.r
```
troubleshooting: if you haven't add the path to Rscript to the environment path, you will need to do that first.
3. If you want to generate DEG analysis output from raw count data, use this command:
```bash
Rscript ./R/findDEGs.R [--batch=<batch> --pval=<pval> --logfc=<logfc> --adjMethod=<adjMethod>] <input> <output> <metadata>
```
use Rscript ./R/findDEGs.R -h to get more information

4. If you want to generate enrichment analysis / gsea output from DEG analysis, use this command:
```bash
Rscript ./R/enrichment_main.R (--enrich | --gsea) [--dePval=<dePval> --pval=<pval> --qval=<qval> --adjMethod=<adjMethod>] <input> <output>
```

use Rscript ./R/enrichment_main.R -h to get more information