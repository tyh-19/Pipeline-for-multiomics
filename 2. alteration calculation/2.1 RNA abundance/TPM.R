#! /usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
parser <- ArgumentParser(description='TPM/CPM')
parser$add_argument('-m', '--matrix', type='character', required=TRUE,
                    help='Input count matrix. Rows are genes(TPM). Columns are samples.')
parser$add_argument('-n', '--normalization', type='character', required=TRUE,
                    choices=c("TPM","CPM"))
parser$add_argument('-I', '--ID', type='character', required=TRUE,
                    choices=c("ensembl_gene_id","hgnc_symbol","external_gene_name","EnsDb.Hsapiens.v86"))
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
                    help='output directory, do not end with /')
args <- parser$parse_args()

##TPM calculator (referred to ensembl ID longest transcript)
  TPM <- function(rawcount,ID,normalization,output,transcript_length) {
    counts <- rawcount
    rm(rawcount)
    
    message("Start preparing human biomart.")
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    message("Finish preparing biomart.")
    message("Start retriving transcript length.")
    if(ID=="ensembl_gene_id"){
    message("Using ensembl_gene_id...")
    values <- as.character(lapply(strsplit(rownames(counts),".",fixed = TRUE),function(x) x[1]))

    query <- unique(values)
    i=1
    annotations <- {}
    while(i<=(length(query)%/%5000)){
    message("Round :", i)
    annotations_temp <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"),
                         filters = "ensembl_gene_id",
                         values=query[5000*(i-1)+1:5000*i], mart= mart)#,useCache = FALSE)
    Sys.sleep(1)
    annotations <- rbind(annotations,annotations_temp)
    i=i+1
    }
    annotations_temp <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"),
                         filters = "ensembl_gene_id",
                         values=query[5000*(i-1)+1:length(query)%%5000], mart= mart)#,useCache = FALSE)
    annotations <- rbind(annotations,annotations_temp)

    message("Finish retriving transcript length.")
    longest_transcript <- arrange(annotations, ensembl_gene_id, dplyr::desc(transcript_length))
    #write.csv(longest_transcript,transcript_length)
    message("Write transcript information to: ",transcript_length)
    #merge 会改变顺序，尽量不要用
    #merged <- merge(as.data.frame(rownames(counts)),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by.x=1,by.y="hgnc_symbol",all.x=TRUE)

    #left_join 不会改变顺序
    merged <- left_join(as.data.frame(values),longest_transcript[!duplicated(longest_transcript$ensembl_gene_id),],by = c("values"="ensembl_gene_id"))
    
    } else if(ID=="hgnc_symbol") {
    message("Using hgnc_symbol...")
    #values <- rownames(counts)
    values <- as.character(lapply(strsplit(rownames(counts),".",fixed = TRUE),function(x) x[1]))

    query <- unique(values)
    i=1
    annotations <- {}
    while(i<=(length(query)%/%5000)){
    message("Round :", i)
    annotations_temp <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"),
                         filters = "hgnc_symbol",
                         values=query[5000*(i-1)+1:5000*i], mart= mart)#,useCache = FALSE)
    Sys.sleep(1)
    annotations <- rbind(annotations,annotations_temp)
    i=i+1
    }
    annotations_temp <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"),
                         filters = "hgnc_symbol",
                         values=query[5000*(i-1)+1:length(query)%%5000], mart= mart)#,useCache = FALSE)
    annotations <- rbind(annotations,annotations_temp)

    message("Finish retriving transcript length.")
    longest_transcript <- arrange(annotations, hgnc_symbol, dplyr::desc(transcript_length))
    #write.csv(longest_transcript,transcript_length)
    message("Write transcript information to: ",transcript_length)
    merged <- left_join(as.data.frame(values),longest_transcript[!duplicated(longest_transcript$hgnc_symbol),],by = c("values"="hgnc_symbol"))
    } else if(ID=="external_gene_name") {
    message("Using external_gene_name...")
    values <- rownames(counts)

    print(head(values,10))
    query <- unique(values)
    i=1
    annotations <- {}
    while(i<=(length(query)%/%5000)){
    message("Round :", i)
    annotations_temp <- getBM(attributes=c("external_gene_name","hgnc_symbol","ensembl_gene_id","transcript_length"),
                         filters = "external_gene_name",
                         values=query[5000*(i-1)+1:5000*i], mart= mart)#,useCache = FALSE)
    Sys.sleep(1)
    annotations <- rbind(annotations,annotations_temp)
    i=i+1
    }
    annotations_temp <- getBM(attributes=c("external_gene_name","hgnc_symbol","ensembl_gene_id","transcript_length"),
                         filters = "external_gene_name",
                         values=query[5000*(i-1)+1:length(query)%%5000], mart= mart)#,useCache = FALSE)
    annotations <- rbind(annotations,annotations_temp)    

    message("Finish retriving transcript length.")
    longest_transcript <- arrange(annotations, external_gene_name, dplyr::desc(transcript_length))
    #write.csv(longest_transcript,transcript_length)
    message("Write transcript information to: ",transcript_length)
    merged <- left_join(as.data.frame(values),longest_transcript[!duplicated(longest_transcript$external_gene_name),],by = c("values"="external_gene_name"))
    #write.csv(merged,transcript_length)
    } else if(ID=="EnsDb.Hsapiens.v86"){
    message("Using EnsDb.Hsapiens.v86...")
    values <- rownames(counts)
    V86_anno <- AnnotationDbi::select(EnsDb.Hsapiens.v86, key=values,columns=c("SYMBOL","GENEID"),keytype="SYMBOL")
    V86_anno <- V86_anno[grep("ENSG",V86_anno$GENEID),]
    Symbol_ensembl <- left_join(as.data.frame(values),V86_anno[!duplicated(V86_anno$SYMBOL),],by = c("values"="SYMBOL"))

    query <- unique(Symbol_ensembl$GENEID)
    i=1
    annotations <- {}
    while(i<=(length(query)%/%5000)){
    message("Round :", i)
    annotations_temp <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"),
                         filters = "ensembl_gene_id",
                         values=query[5000*(i-1)+1:5000*i], mart= mart)#,useCache = FALSE)
    Sys.sleep(1)
    annotations <- rbind(annotations,annotations_temp)
    i=i+1
    }
    annotations_temp <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","transcript_length"),
                         filters = "ensembl_gene_id",
                         values=query[5000*(i-1)+1:length(query)%%5000], mart= mart)#,useCache = FALSE)
    annotations <- rbind(annotations,annotations_temp)
    
    annotations <- left_join(Symbol_ensembl,annotations,by=c("GENEID"="ensembl_gene_id"))

    message("Finish retriving transcript length.")
    longest_transcript <- arrange(annotations, values, dplyr::desc(transcript_length))

    #write.csv(longest_transcript,transcript_length)
    message("Write transcript information to: ",transcript_length)
    merged <- left_join(as.data.frame(values),longest_transcript[!duplicated(longest_transcript$values),],by = c("values"="values"))
    colnames(merged) <- gsub("GENEID","ensembl_gene_id",colnames(merged))
    }
    
    #对于无注释长度的转录本，记为1000，即不对长度标准化
    merged$transcript_length[is.na(merged$transcript_length)] <- 1000
    
    transcript_lengths <- as.integer(merged$transcript_length)
    
    write.csv(merged,transcript_length)
    message("Finish transcript information.")
    
    if(normalization=="TPM"){
    #find gene length normalized values 
    rpk <- apply(counts, 2, function(x) x/(transcript_lengths/1000))
    write.table(rpk,paste0(output_dir,"/RPK.txt"),quote = FALSE,sep="\t")
    message("Finish RPK calculating.")
    rm(counts)
    gc()
    
    #rpk_out <- rpk
    #rownames(rpk_out) <- new_rownames$new_rownames
    #rpk_out <- as.data.frame(rpk_out)
    #rownames(rpk_out) <- gsub(".","|",fixed = TRUE,rownames(rpk_out))
    #rpk_out$gene_id <- new_rownames$new_rownames
    #rpk_out <- aggregate(. ~ gene_id, data = rpk_out, sum)
    #rownames(rpk_out) <- rpk_out$gene_id
    #rpk_out <- rpk_out[,-which(colnames(rpk_out)=="gene_id")]
    #write.table(rpk_out,paste0(output_dir,"/RPK.txt"),quote = FALSE,sep="\t")
    #rm(rpk_out)
    #gc()

    #normalize by the sample size using rpk values
    
    if(ncol(rpk) <= 100000){
    tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
    message("Finish TPM calculating.")
    rm(rpk)
    gc()
   
    tpm <- as.data.frame(tpm)
    new_rownames <- unite(merged, "new_rownames", values, ensembl_gene_id, hgnc_symbol, transcript_length, sep = "|")
    rownames(tpm) <- new_rownames$new_rownames
    write.table(tpm,output,quote = FALSE,sep="\t")
    } else {
    rpk <- as.data.frame(t(t(rpk)/colSums(rpk)))
    new_rownames <- unite(merged, "new_rownames", values, ensembl_gene_id, hgnc_symbol, transcript_length, sep = "|")
    rownames(rpk) <- new_rownames$new_rownames
    write.table(rpk,output,quote = FALSE,sep="\t")
    }
    } else if(normalization=="CPM"){
    cpm <- sapply(counts, function(x) x / sum(as.numeric(x)) * 10^6)
    message("Finish CPM calculating.")
    rm(counts)
    gc()
    cpm <- as.data.frame(cpm)
    new_rownames <- unite(merged, "new_rownames", values, ensembl_gene_id, hgnc_symbol, transcript_length, sep = "|")
    message("Finish new rownames generation.")
    print(head(new_rownames,10))
    message("nrow(cpm): ",nrow(cpm))
    message("length(new_rownames$new_rownames)",length(new_rownames$new_rownames))
    rownames(cpm) <- new_rownames$new_rownames
    write.table(cpm,output,quote = FALSE,sep="\t")
    }
    
  }

message("Input raw count matrix: ",args$matrix)
message("Identifier: ",args$ID)
message("Output directory: ",args$outdir)

message("Start reading raw count matrix.")
rawcount <- read.table(args$matrix,sep = "\t",header = TRUE,row.names = 1)
message("Finish reading raw count matrix.")
ID <- args$ID
output_dir <- args$outdir
dir.create(output_dir)
normalization <- args$normalization
if(normalization=="TPM"){
output <- paste0(output_dir,"/TPM.txt")
}else if(normalization=="CPM"){
output <- paste0(output_dir,"/CPM.txt")
}
transcript_length <- paste0(output_dir,"/Longest_transcript_length.csv")

message("Start TPM converting.")
TPM(rawcount,ID,normalization,output,transcript_length)
message("Finish TPM converting.")
