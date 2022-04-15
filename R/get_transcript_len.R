#!/usr/bin/env Rscript

'Generate an exon table from gff3 or gtf file
Usage:
    get_transcript_len.R <input> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.

Arguments:
    input  input gtf or gff3 file name
    output  output file name
' -> doc
suppressMessages(library(docopt))
suppressMessages(library(GenomicFeatures))

args <- docopt(doc)

txdb <- makeTxDbFromGFF(file.path(args$input))
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))

write.csv(exonic.gene.sizes, file.path(args$output))