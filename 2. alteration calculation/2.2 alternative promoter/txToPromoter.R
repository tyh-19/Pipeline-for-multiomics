#!/usr/bin/env Rscript
message("Load required packages ...")
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(parallel))
message("Done .")

parser <- ArgumentParser(description='Group transcripts of each gene by their TSS')
parser$add_argument('-a', '--gtf', type='character', required=TRUE,help='Input gtf file')
parser$add_argument('-g','--gap-min',type='integer',default=10,help='Min gap for not collapsing two TSS to one promoter')
parser$add_argument('-o', '--output', type='character',required=TRUE,help='Output')
parser$add_argument('-p', '--threads', type='integer',default=1,help='Thread of lapply')
args <- parser$parse_args()

tssGrouper <- function(txByGene,min.gap = 10){
  ## Get unique TSS site for gene at +/- sttrand
  strandness <- as.character(runValue(strand(txByGene)))
  if(strandness=="+"){
    tss <- start(txByGene)
    tss.unique <- sort(unique(tss),decreasing = F)
  }else{
    tss <- end(txByGene)
    tss.unique <- sort(unique(tss),decreasing = T)
  }
  mcols(txByGene)[["tss.pos"]] <- paste(tss,strandness,sep=".")

  n.tss <- length(tss.unique)
  last.tss.pos <- tss.unique[1]

  if(n.tss>1){
   ## If more than 1 TSS
   tss.names <- c(paste(last.tss.pos,strandness,sep="."))
   for(i in 2:n.tss){
     if(abs(tss.unique[i] - last.tss.pos) > min.gap){
       last.tss.pos <- tss.unique[i]
     }
     tss.names <- c(tss.names,paste(last.tss.pos,strandness,sep="."))
   }
  names(tss.names) <- as.character(tss.unique)
  mcols(txByGene)[["tss.group"]] <- tss.names[as.character(tss)]
  }else{
    mcols(txByGene)[["tss.group"]] <- paste(last.tss.pos,strandness,sep=".")
  }
  
  as.data.frame(mcols(txByGene)[,c("tx_name","tss.group","tss.pos")])

}


message("Load genome annotations ...")
txdb <- makeTxDbFromGFF(file=args$gtf,format="gtf")
message("Done .")

message("Assign transcripts to TSS")
message("Use ",args$threads," cores for parallel processing")
txByGenes <- transcriptsBy(txdb, by = 'gene')
txByGenes.tss <- mclapply(txByGenes, tssGrouper, args$gap_min, mc.cores=args$threads)
txByGenes.tss.df <- dplyr::bind_rows(txByGenes.tss, .id = "gene.id")
message("Done .")

message("Write TSS annotation . ")
write.table(txByGenes.tss.df, file=args$output, sep="\t", row.names = F, quote = F)
message("Done .")
