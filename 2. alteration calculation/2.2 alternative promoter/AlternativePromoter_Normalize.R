 #! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(doParallel))
parser <- ArgumentParser(description='Normalize alternative promoter usage for each trasncript, TPM of Transcripts shared same promoter / TPM of all transcripts of one gene')
parser$add_argument('-i', '--input_file', type='character', required=TRUE,
    help='input TPM-by-promoter_${dataset}.txt')
parser$add_argument('-o', '--output_file', type='character', required=TRUE,
    help='output file for normalized promoter usage, AlternativePromoter_${dataset}_normalized.txt')
args <- parser$parse_args()

input <- args$input_file
output <- args$output_file

Alt.promoter <- read.csv(input,header = TRUE, sep = "\t", row.names = 1)

Alt.promoter$gene <- as.character(lapply(strsplit(rownames(Alt.promoter),"|",fixed = TRUE),function(x) x[1]))

Alt.promoter.base <- aggregate(. ~ gene, data = Alt.promoter, sum)

Alt.promoter.tmp <- as.data.frame(Alt.promoter[,"gene"])
colnames(Alt.promoter.tmp) <- "gene"

Alt.promoter.base.all <- left_join(Alt.promoter.tmp,Alt.promoter.base, by = c("gene"="gene"))

Alt.promoter.normalized <- Alt.promoter[,-ncol(Alt.promoter)]/Alt.promoter.base.all[,-1]

Alt.promoter.normalized[is.na(Alt.promoter.normalized)] <- 0
colnames(Alt.promoter.normalized) <- gsub(".","-",fixed=TRUE,colnames(Alt.promoter.normalized))
write.table(Alt.promoter.normalized,output,sep = "\t",quote=FALSE)
