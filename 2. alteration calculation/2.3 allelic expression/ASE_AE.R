 #! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(doParallel))
parser <- ArgumentParser(description='Summarize ASE allele frequency from ASEcount. AE = |0.5-ref/total|')
parser$add_argument('-i', '--input_dir', type='character', required=TRUE,
    help='input directory containing all vcf wait to be summarized(do not end with /)')
parser$add_argument('-o', '--output_file', type='character', required=TRUE,
    help='output file for AE count(do not end with /)')
args <- parser$parse_args() 

  input_dir <- args$input_dir
  output_file <- args$output_file
  files <- dir(input_dir)
  
  cl <- makeCluster(16) #not to overload your computer
  registerDoParallel(cl) 

  i=1
  final_SNP_AE <- as.data.frame(matrix(factor(0),ncol = 1,nrow = 0))
  colnames(final_SNP_AE) <- c("SNP")
  
  message('Number of files to summarize:', length(files))
  while(i<=length(files)) {
  SNP <- read.table(paste0(input_dir,"/",files[i]),header = TRUE)
  sample <- as.character(lapply(strsplit(files[i],".",fixed = 1),function(x) x[1]))
  SNP_passed <- SNP[SNP$totalCount > 15,]
  
  SNP_rowname <- paste(SNP_passed$contig,SNP_passed$position,SNP_passed$variantID,SNP_passed$refAllele,SNP_passed$altAllele, sep = "_")
  rownames(SNP_passed) <- SNP_rowname
  
  SNP_AE <- data.frame("SNP"=SNP_rowname,"AE"=abs(0.5-SNP_passed$refCount/SNP_passed$totalCount))
  colnames(SNP_AE) <- gsub("AE",sample,colnames(SNP_AE))
  
  final_SNP_AE <- full_join(final_SNP_AE,SNP_AE,by=c("SNP"="SNP"))
  message('Summarized files Number:', i)
  i=i+1
  }

  final_SNP_AE[is.na(final_SNP_AE)] <- 0.5
  
  rownames(final_SNP_AE) <- final_SNP_AE$SNP
 
  write.table(final_SNP_AE[-which(colnames(final_SNP_AE)=="SNP")],output_file,sep = "\t",quote = FALSE)

  stopCluster(cl)
