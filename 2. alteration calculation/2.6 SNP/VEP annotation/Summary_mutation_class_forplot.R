#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))


parser <- ArgumentParser(description='Mutation classification summary')
parser$add_argument('-i', '--input_dir', type='character', required=TRUE,
    help='input directory, which contains VEP annoted SNP files(end with .VEP.txt)')
parser$add_argument('-o', '--output_dir', type='character', required=TRUE,
    help='output directory.')
args <- parser$parse_args()

  
  input_dir <- args$input_dir
  output_dir <- args$output_dir
  files <- dir(input_dir)
  files <- grep("VEP",files,value = TRUE)
  
  i=1
  mutation_burden_matrix <- as.data.frame(matrix(numeric(0),ncol=1))
  colnames(mutation_burden_matrix) <- c("Mutated_Gene") 
  mutation_burden_matrix$`Mutated_Gene` <- as.factor(mutation_burden_matrix$`Mutated_Gene`)
  while(i<=length(files)){
  vcf <- read.table(paste0(input_dir,"/",files[i]))
  sample <- as.character(lapply(strsplit(files[i],".",fixed = TRUE),function(x) x[1]))
  colnames(vcf) <- c("Uploaded_variation","Location","Allele","Gene","Feature","Feature_type","Consequence","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","IMPACT","DISTANCE","STRAND","FLAGS","VARIANT_CLASS","SYMBOL","SYMBOL_SOURCE","HGNC_ID","BIOTYPE","CANONICAL","MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO","SIFT","PolyPhen","EXON","INTRON","DOMAINS","miRNA","HGVSc","HGVSp","HGVS_OFFSET","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS")
  #filter COSMIC genes
  vcf_filtered <- vcf[grep("COSV",vcf$Existing_variation),]
  
  df <- data.frame("Mutated_Gene" = paste(vcf_filtered$Gene,vcf_filtered$SYMBOL,sep="|"),"Variant_Class" = vcf_filtered$VARIANT_CLASS, "Consequence" = vcf_filtered$Consequence)
  mutation_class <- aggregate(Variant_Class ~ Mutated_Gene, df, paste, collapse = ",")
  mutation_class$Variant_Class <- sapply(strsplit(mutation_class$Variant_Class, ","), function(x) paste(rle(x)$values, collapse=";"))
  mutation_class$Variant_Class <- paste0(mutation_class$Variant_Class,";")
  colnames(mutation_class) <- c("Mutated_Gene",sample)
  ##mutation_count <- as.data.frame(setDT(df)[,list(Count=.N),names(df)]) ###
  ##colnames(mutation_count) <- c("Mutated_Gene",sample)
  mutation_burden_matrix <- full_join(mutation_burden_matrix,mutation_class,by=c("Mutated_Gene"="Mutated_Gene"))
  
  i=i+1
  }
  
  #mutation_burden_matrix[is.na(mutation_burden_matrix)] <- 0
  rownames(mutation_burden_matrix) <- mutation_burden_matrix$`Mutated_Gene`
  mutation_burden_matrix <- mutation_burden_matrix[,-which(colnames(mutation_burden_matrix)=="Mutated_Gene")]
  nonsense_row <- which(rownames(mutation_burden_matrix)=="-|-")
  if(length(nonsense_row)!=0){
  mutation_burden_matrix <- mutation_burden_matrix[-which(rownames(mutation_burden_matrix)=="-|-"),]
  } else {
    mutation_burden_matrix <- mutation_burden_matrix
  }
  #mutation_burden_matrix <- as.data.frame(t(mutation_burden_matrix))
  #mutation_burden_matrix$Case.ID <- rownames(mutation_burden_matrix)
  #rownames(mutation_burden_matrix) <- seq(from=1,to=nrow(mutation_burden_matrix))
  write.table(mutation_burden_matrix,paste0(output_dir,"/Mutation_burden_class_forplot.txt"),sep = "\t",quote = FALSE)
