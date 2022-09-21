 #! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(doParallel))
parser <- ArgumentParser(description='Summarize star fusion counts, including JunctionReadCount and SpanningFragCount.')
parser$add_argument('-i', '--input_dir', type='character', required=TRUE,
    help='input directory containing all star-fusion.fusion_predictions.tsv waited to be summarized (do not end with /)')
parser$add_argument('-o', '--output_dir', type='character', required=TRUE,
    help='output directory for chimric count and spanning count (do not end with /)')
args <- parser$parse_args() 

  input_dir <- args$input_dir
  output_dir <- args$output_dir
  dir.create(output_dir)
  files <- dir(input_dir)
  
  cl <- makeCluster(16) #not to overload your computer
  registerDoParallel(cl) 

  i=1
  final_JunctionReadCount <- as.data.frame(matrix(factor(0),ncol = 1,nrow = 0))
  colnames(final_JunctionReadCount) <- c("ChimericRNA")
  final_SpanningFragCount <- as.data.frame(matrix(factor(0),ncol = 1,nrow = 0))
  colnames(final_SpanningFragCount) <- c("ChimericRNA")

  message('Number of files to summarize:', length(files))
  while(i<=length(files)) {
  ChimericCount <- read.table(paste0(input_dir,"/",files[i],"/star-fusion.fusion_predictions.tsv"),header = TRUE,check.names=FALSE,comment.char = "")
  sample <- as.character(lapply(strsplit(files[i],".",fixed = 1),function(x) x[1]))
  ChimericCount_passed <- ChimericCount[ChimericCount$FFPM > 0.1 & ChimericCount$LargeAnchorSupport == "YES_LDAS",]
  

  ChimericCount_rowname <- paste(ChimericCount_passed$`#FusionName`,ChimericCount_passed$LeftGene,ChimericCount_passed$RightGene,ChimericCount_passed$LeftBreakpoint,ChimericCount_passed$RightBreakpoint,ChimericCount_passed$SpliceType,ChimericCount_passed$annots, sep = "|")
  rownames(ChimericCount_passed) <- ChimericCount_rowname
  
  JunctionReadCount <- data.frame("ChimericRNA"=ChimericCount_rowname,"JunctionReadCount"=ChimericCount_passed$JunctionReadCount)
  colnames(JunctionReadCount) <- gsub("JunctionReadCount",sample,colnames(JunctionReadCount))
  final_JunctionReadCount <- full_join(final_JunctionReadCount,JunctionReadCount,by=c("ChimericRNA"="ChimericRNA"))

  SpanningFragCount <- data.frame("ChimericRNA"=ChimericCount_rowname,"SpanningFragCount"=ChimericCount_passed$SpanningFragCount)
  colnames(SpanningFragCount) <- gsub("SpanningFragCount",sample,colnames(SpanningFragCount))
  final_SpanningFragCount <- full_join(final_SpanningFragCount,SpanningFragCount,by=c("ChimericRNA"="ChimericRNA"))

  message('Summarized files Number:', i)
  i=i+1
  }
  
  final_JunctionReadCount[is.na(final_JunctionReadCount)] <- 0
  final_SpanningFragCount[is.na(final_SpanningFragCount)] <- 0
  rownames(final_JunctionReadCount) <- final_JunctionReadCount$ChimericRNA
  write.table(final_JunctionReadCount[-which(colnames(final_JunctionReadCount)=="ChimericRNA")],paste0(output_dir,"/JunctionReadCount.txt"),sep = "\t",quote = FALSE)

  rownames(final_SpanningFragCount) <- final_SpanningFragCount$ChimericRNA
  write.table(final_SpanningFragCount[-which(colnames(final_SpanningFragCount)=="ChimericRNA")],paste0(output_dir,"/SpanningFragCount.txt"),sep = "\t",quote = FALSE)
  stopCluster(cl)
