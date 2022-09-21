
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))


parser <- ArgumentParser(description='Mutation gene fisher exact test')
parser$add_argument('-i', '--input_matrix', type='character', required=TRUE,
    help='input gene mutated site count matrix.')
parser$add_argument('-p', '--positive', type='character', required=TRUE,
    choices=c('CRC','STAD','NC','CRC|STAD'))
parser$add_argument('-n','--negative', type='character', required=TRUE,
    choices=c('CRC','STAD','NC','CRC|STAD'))
parser$add_argument('-o', '--output_dir', type='character', required=TRUE,
    help='output directory.')
args <- parser$parse_args()

  mutation_burden_matrix <- read.csv(args$input_matrix, header = TRUE, row.names = 1,sep = "\t")
  positive <- args$positive
  negative <- args$negative
  output_dir <- args$output_dir
  
  positive_id <- grep(positive,colnames(mutation_burden_matrix),value = TRUE)
  negative_id <- grep(negative,colnames(mutation_burden_matrix),value = TRUE)
  
  message('Positive sample number: ',length(positive_id))
  message('Negative sample number: ',length(negative_id))
 
  Mutated_gene_DE <- data.frame(matrix(numeric(0),ncol=5,nrow=nrow(mutation_burden_matrix)))
  colnames(Mutated_gene_DE) <- c("Gene","Positive group mutation ratio","Negative group mutation ratio","pvalue","padj")  
  i=1
  #total <- nrow(mutation_burden_matrix)
  while(i<=nrow(mutation_burden_matrix)){
  #message('Progressing: ',i,"/",total)
  
  get_zero_number <-function(x) sum(x==0)
  
  Positive_no_mutated <- apply(mutation_burden_matrix[i,positive_id],1,get_zero_number)
  Negative_no_mutated <- apply(mutation_burden_matrix[i,negative_id],1,get_zero_number)
  Positive_mutated <- ncol(mutation_burden_matrix[i,positive_id])-Positive_no_mutated
  Negative_mutated <- ncol(mutation_burden_matrix[i,negative_id])-Negative_no_mutated
  
  fisher_exact_table <- data.frame(positive=c(Positive_mutated,Positive_no_mutated),negative = c(Negative_mutated,Negative_no_mutated))
  rownames(fisher_exact_table) <- c("Mutated","No_mutation_detected")
  fisher_result <- fisher.test(fisher_exact_table)
  Mutated_gene_DE[i,"pvalue"] <- fisher_result$p.value
  
  Mutated_gene_DE[i,"Positive group mutation ratio"] <- Positive_mutated/(Positive_mutated+Positive_no_mutated)
  Mutated_gene_DE[i,"Negative group mutation ratio"] <- Negative_mutated/(Negative_mutated+Negative_no_mutated)
  
  Mutated_gene_DE[i,"Gene"] <- rownames(mutation_burden_matrix[i,])
  
  i=i+1
  }
  Mutated_gene_DE$padj <- p.adjust(Mutated_gene_DE$pvalue,method = "BH")
  rownames(Mutated_gene_DE) <- Mutated_gene_DE$Gene
  write.table(Mutated_gene_DE[,-which(colnames(Mutated_gene_DE)=="Gene")],paste0(output_dir,"/","DE_mutated_gene_fisher_exact_",positive,"vs",negative,".txt"),sep = "\t",quote = FALSE)
  
