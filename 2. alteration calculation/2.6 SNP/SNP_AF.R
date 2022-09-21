 #! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(doParallel))
parser <- ArgumentParser(description='Summarize SNP allele frequency from vcf.')
parser$add_argument('-i', '--input_dir', type='character', required=TRUE,
    help='input directory containing all vcf wait to be summarized(do not end with /)')
parser$add_argument('-o', '--output_dir', type='character', required=TRUE,
    help='output directory for Alt and Ref count(do not end with /)')
args <- parser$parse_args() 

  input_dir <- args$input_dir
  output_dir <- args$output_dir
  files <- dir(input_dir)
  files <- grep("rmEDIT.SNP",files,value = TRUE)
  #positive_id <- grep(positive,as.character(lapply(strsplit(files,".",fixed = 1),function(x) x[1])),value = TRUE)
  #negative_id <- grep(negative,as.character(lapply(strsplit(files,".",fixed = 1),function(x) x[1])),value = TRUE)
  
  cl <- makeCluster(16) #not to overload your computer
  registerDoParallel(cl) 

  i=1
  final_SNP_Ref <- as.data.frame(matrix(factor(0),ncol = 1,nrow = 0))
  colnames(final_SNP_Ref) <- c("SNP")
  final_SNP_Alt <- as.data.frame(matrix(factor(0),ncol = 1,nrow = 0))
  colnames(final_SNP_Alt) <- c("SNP")
  
  message('Number of files to summarize:', length(files))
  while(i<= length(files)) {
  SNP <- read.table(paste0(input_dir,"/",files[i]))
  sample <- as.character(lapply(strsplit(files[i],".",fixed = 1),function(x) x[1]))
  colnames(SNP) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","VALUE")
  SNP_passed <- SNP[SNP$FILTER=="PASS",]
  
  SNP_rowname <- paste(SNP_passed$CHROM,SNP_passed$POS,SNP_passed$REF,SNP_passed$ALT, sep = "_")
  rownames(SNP_passed) <- SNP_rowname
  
  Ref_Alt <- as.character(lapply(strsplit(as.character(SNP_passed[,"VALUE"]),":"),function(x) x[2]))
  
  Ref <- as.integer(as.character(lapply(strsplit(Ref_Alt,","),function(x) x[1])))
  Alt1 <- as.integer(as.character(lapply(strsplit(Ref_Alt,","),function(x) x[2])))
  Alt2 <- as.integer(as.character(lapply(strsplit(Ref_Alt,","),function(x) x[3])))
  Alt2[is.na(Alt2)] <- 0
  Alt <- Alt1+Alt2
  
  SNP_Ref <- data.frame("SNP"=SNP_rowname,"Ref"=Ref)
  colnames(SNP_Ref) <- c("SNP",paste0(sample,"_Ref"))
  final_SNP_Ref <- full_join(final_SNP_Ref,SNP_Ref,by=c("SNP"="SNP"))
  
  SNP_Alt <- data.frame("SNP"=SNP_rowname,"Alt"=Alt)
  colnames(SNP_Alt) <- c("SNP",paste0(sample,"_Alt"))
  final_SNP_Alt <- full_join(final_SNP_Alt,SNP_Alt,by=c("SNP"="SNP"))
  message('Summarized files Number:', i)
  i=i+1
  }

  final_SNP_Alt[is.na(final_SNP_Alt)] <- 0
  final_SNP_Ref[is.na(final_SNP_Ref)] <- 0
  
  rownames(final_SNP_Alt) <- final_SNP_Alt$SNP
  rownames(final_SNP_Ref) <- final_SNP_Ref$SNP
  
  Final <- final_SNP_Alt[,-which(colnames(final_SNP_Alt)=="SNP")]/(final_SNP_Alt[,-which(colnames(final_SNP_Alt)=="SNP")]+final_SNP_Ref[,-which(colnames(final_SNP_Ref)=="SNP")])
  
  colnames(Final) <- gsub("_Alt","",colnames(Final))
  Final[is.na(Final)] <- 0

  write.table(final_SNP_Alt[,-which(colnames(final_SNP_Alt)=="SNP")],paste0(output_dir,"/","final_SNP_Alt_count.txt"),sep = "\t",quote = FALSE)
  write.table(final_SNP_Ref[,-which(colnames(final_SNP_Ref)=="SNP")],paste0(output_dir,"/","final_SNP_Ref_count.txt"),sep = "\t",quote = FALSE)
  write.table(Final,paste0(output_dir,"/","final_SNP_AF.txt"),sep = "\t",quote = FALSE)

  stopCluster(cl)
