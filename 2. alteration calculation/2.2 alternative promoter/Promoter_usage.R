#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(doParallel))

parser <- ArgumentParser(description='Alternative promoter usage: major, minor and inactivate(TPM < 1).')
parser$add_argument('-n', '--input_normalized', type='character', required=TRUE,
    help='input promoter normalized usage')
parser$add_argument('-t', '--input_TPM', type='character', required=TRUE,
    help='input TPM-by-promoter_${dataset}.txt')
parser$add_argument('-s', '--sample', type='character', required=TRUE,
    help='input 1 sample')
parser$add_argument('-o', '--output', type='character', required=TRUE,
    help='output dir(do not end with /)')
args <- parser$parse_args()

      input_normalized <- args$input_normalized
      input_TPM <- args$input_TPM
      Altpromoter_normalized <- read.table(input_normalized,sep = "\t", header = TRUE,check.names=FALSE,row.names = 1)
      Altpromoter_activity <- read.table(input_TPM,sep = "\t", header = TRUE,check.names=FALSE,row.names = 1)
     
      output <- args$output
      dir.create(output)

      sample <- args$sample 
      genes <- unique(as.character(lapply(strsplit(rownames(Altpromoter_normalized),".",fixed = TRUE), function(x) x[1])))

        single_sample_promoter_usage <- data.frame(matrix(nrow=3,ncol = length(genes)))
        colnames(single_sample_promoter_usage) <- genes
        rownames(single_sample_promoter_usage) <- paste0(c("inactivate","major","minor"),"_",sample)
        for(gene in genes){
	message(gene)
          single_ratio <- Altpromoter_normalized[grep(gene,rownames(Altpromoter_normalized)),sample,drop=FALSE]
          single_promoter_activity <- Altpromoter_activity[grep(gene,rownames(Altpromoter_activity)),sample,drop=FALSE]
          if(length(single_ratio[single_promoter_activity>=1])==0){
	    inactivate <- 1
            major <- 0
            minor <- 0
          } else {
	  inactivate <- sum(single_ratio[single_promoter_activity<1])
          major <- max(single_ratio[single_promoter_activity>=1])
          minor <- sum(single_ratio[single_promoter_activity>=1 & single_promoter_activity<max(single_promoter_activity[single_promoter_activity>=1])])
          }
          single_sample_promoter_usage[paste0("inactivate_",sample),gene] <- inactivate
          single_sample_promoter_usage[paste0("major_",sample),gene] <- major
          single_sample_promoter_usage[paste0("minor_",sample),gene] <- minor
        }
        single_sample_promoter_usage$sample <- sample
        single_sample_promoter_usage$promoter <- c("inactivate","major","minor")
        write.table(single_sample_promoter_usage,paste0(output,"/",sample,".txt"),sep = "\t",quote = FALSE)
