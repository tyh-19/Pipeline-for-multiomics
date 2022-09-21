#library preparation
{
  library(reshape2) #for melt
  library(edgeR) #for RNA expression and DNA methylation differential analysis
  library(ggplot2) #for plot
  library(ggpubr) #for ggarrange
  library(sunburstR) #sunburst plot
  library(htmlwidgets) #for html sunburst plot
  library(extrafont) #for plot font load
  fonts() #load fonts
  library(ggsci) #for plot color theme
  library(ggpattern) #for pattern in geom_bar
  library(ggplotify) #for complexheatmap2ggplot
  library(cowplot) #for plot arrange
  library(tidyverse) #for dataframe process
  library(FactoMineR) #plot PCA
  library(factoextra) #plot PCA
  library(scatterplot3d) #plot 3D PCA
  library(progress) #for progress bar
  library(dplyr) # for data manipulation
  library(clusterProfiler) #for enrichment analysis
  library(biomaRt) #for gene annotations
  library(multiMiR) #for miRNA target gene predict
  library(ggradar)
  library(resample)
  library(chimeraviz)
  library(ComplexHeatmap) #complexheatmap for differential alterations
  library(tidyr) #for gather (from data frame to a long table)
  
  #library biomart
  {
    library(biomaRt)
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  }
}

#function
{
  #PCA_group_distance
  PCA_group_distance <- function(for_distance,Group1,Group2){
    message("for_distance must have Dim.1, Dim.2,Group column, rownames should be sample(point) names.")
    in_group <- for_distance[for_distance$Group==Group1,]
    out_group <- for_distance[for_distance$Group==Group2,]
    i=1
    j=1
    all_in_group_distance <- {}
    while(i<=nrow(in_group)){
      while(j<=nrow(in_group)-i){
        message(rownames(in_group)[i],"_and_",rownames(in_group)[i+j])
        in_group_distance <- in_group[i,c("Dim.1","Dim.2")]-in_group[i+j,c("Dim.1","Dim.2")]
        single_in_group_distance <- data.frame("point A"=rownames(in_group)[i],"point B"=rownames(in_group)[i+j],"distance"=sqrt(rowSums(in_group_distance*in_group_distance)))
        all_in_group_distance <- rbind(all_in_group_distance,single_in_group_distance)
        j=j+1
      }
      i=i+1
    }
    rownames(all_in_group_distance) <- NULL
    mean_in_group_distance <- mean(all_in_group_distance$distance)
    
    i=1
    j=1
    all_out_group_distance <- {}
    while(i<=nrow(out_group)){
      while(j<=nrow(out_group)-i){
        message(rownames(out_group)[i],"_and_",rownames(out_group)[i+j])
        out_group_distance <- out_group[i,c("Dim.1","Dim.2")]-out_group[i+j,c("Dim.1","Dim.2")]
        single_out_group_distance <- data.frame("point A"=rownames(out_group)[i],"point B"=rownames(out_group)[i+j],"distance"=sqrt(rowSums(out_group_distance*out_group_distance)))
        all_out_group_distance <- rbind(all_out_group_distance,single_out_group_distance)
        j=j+1
      }
      i=i+1
    }
    rownames(all_out_group_distance) <- NULL
    mean_out_group_distance <- mean(all_out_group_distance$distance)
    
    i=1
    j=1
    all_between_group_distance <- {}
    while(i<=nrow(in_group)){
      while(j<=nrow(out_group)){
        message(rownames(in_group)[i],"_and_",rownames(out_group)[j])
        between_group_distance <- in_group[i,c("Dim.1","Dim.2")]-out_group[j,c("Dim.1","Dim.2")]
        single_between_group_distance <- data.frame("point A"=rownames(in_group)[i],"point B"=rownames(out_group)[j],"distance"=sqrt(rowSums(between_group_distance*between_group_distance)))
        all_between_group_distance <- rbind(all_between_group_distance,single_between_group_distance)
        j=j+1
      }
      i=i+1
    }
    rownames(all_between_group_distance) <- NULL
    mean_between_group_distance <- mean(all_between_group_distance$distance)
    
    output <- data.frame("Group"=c(Group1,Group2,paste0(Group1," vs. ",Group2)),"Distance"=c(mean_in_group_distance,mean_out_group_distance,mean_between_group_distance))
  }
  
  #sigmoid
  sigmoid = function(x) {
    1 / (1 + exp(-x))
  }
  
  ##preparation
  library(pspearman)
  spearman_CI <- function(x, y, alpha = 0.05){
    rs <- cor(x, y, method = "spearman", use = "complete.obs")
    n <- sum(complete.cases(x, y))
    sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
  }
  
  #enmuerate correlation calculation
  enmuerate_correlation <- function(forcor,cor_methods,output_dir,max = ncol(forcor)-1){
    ##initiation
    #output
    #output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation/"
    #method
    #cor_methods <- "spearman"
    #read in data matrix, row: sample, column: numeric score. 
    #forcor <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Correlation/Correlation between plasma expression and PBMC immune fraction.csv",header = TRUE, row.names = 1)
    #forcor <- forCorrlation
    dir.create(output_dir)
    message("Caculate R for: ",max,"(Number of interested subjects)")
    
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = max, clear = FALSE, width= 60)
    
    result_final <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    
    result <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    i=1
    while(i<= max){
      message(colnames(forcor[,i,drop = FALSE]))
      j=1
      while(j<=(ncol(forcor)-i)){
        message(colnames(forcor[,i,drop = FALSE])," vs ",colnames(forcor[,j+i,drop = FALSE]))
        test <- data.frame("A"=forcor[,i],"B"=forcor[,j+i])
        test[test=="x"] <- NA
        test <- na.omit(test)
        if(nrow(test)==0) {
          result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),"Not available","Not available","Not available","Not available","Not available","Not available")
          colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
          result <- rbind(result,result_tmp)
          j=j+1
        } else {
          if(cor_methods=="pearson"){
            r_twosided <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods)
            #r_twosided <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods)
            if (is.na(r_twosided$estimate)) {
              r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "two.sided")
              #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
            } else if (r_twosided$estimate>0){
              r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "greater")
              #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "greater")
            } else if (r_twosided$estimate<0) {
              r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "less")
              #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "less")
            } else {
              r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "two.sided")
              #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
            }
            result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,r$conf.int[1],r$conf.int[2],r_twosided$p.value)
          } else if(cor_methods=="spearman") {
            r_twosided <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
            #r_twosided <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
            if (is.na(r_twosided$estimate)) {
              r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
              #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
            } else if(r_twosided$estimate>0){
              r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "greater", approximation = "exact")
              #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "greater", approximation = "exact")
            } else if (r_twosided$estimate<0) {
              r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "less", approximation = "exact")
              #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "less", approximation = "exact")
            } else {
              r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
              #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
            }
            result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,spearman_CI(as.numeric(test$A),as.numeric(test$B))[1],spearman_CI(as.numeric(test$A),as.numeric(test$B))[2],r_twosided$p.value)
            #result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,spearman_CI(forcor[,i],forcor[,i+j])[1],spearman_CI(forcor[,i],forcor[,i+j])[2],r_twosided$p.value)
          }
          
          colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
          result <- rbind(result,result_tmp)
          j=j+1
        }
      }
      
      if(i%%20==0){
        result_final <- rbind(result_final,result)
        write.csv(result,paste0(output_dir,i/20,"_Result20.csv"))
        result <- as.data.frame(matrix(numeric(0),ncol=7))
        colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        i=i+1
      } else if(i==(ncol(forcor)-1)){
        result_final <- rbind(result_final,result)
        write.csv(result,paste0(output_dir,i%/%20+1,"_Result",i%%20,".csv"))
        i=i+1
      } else {
        i=i+1
      }
      pb$tick()
      Sys.sleep(1 / 100)
    }
    write.csv(result_final,paste0(output_dir,"Result_final.csv"))
    
    result_final2 <- as.data.frame(matrix(numeric(0),ncol = ncol(forcor),nrow = ))
  }
  
  #paired correlation
  paired_correlation <- function(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,pvalue_cutoff=0.1){
    message("Please make sure input column are samples ids, line are gene ids.")
    forcor1 <- t(forcor1[gene_ids,sample_ids])
    forcor2 <- t(forcor2[gene_ids,sample_ids])
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = ncol(forcor1), clear = FALSE, width= 60)
    
    result_final <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result_final) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    
    result <- as.data.frame(matrix(numeric(0),ncol=7))
    colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
    
    j=1
    while(j<=(ncol(forcor1))){
      #message(j)
      #message(colnames(forcor[,i,drop = FALSE])," vs ",colnames(forcor[,j+i,drop = FALSE]))
      test <- data.frame("A"=forcor1[,j],"B"=forcor2[,j])
      test <- test[order(test$A),]
      #test[test=="x"] <- NA
      test <- na.omit(test)
      test$ID <- factor(rownames(test),levels = rownames(test))
      
      if(nrow(test)==0) {
        result_tmp <- data.frame(paste0(colnames(forcor1)[j]," vs ",colnames(forcor2)[j]),"Not available","Not available","Not available","Not available","Not available","Not available")
        colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        result <- rbind(result,result_tmp)
        j=j+1
      } else {
        if(cor_methods=="pearson"){
          r_twosided <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods)
          #r_twosided <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods)
          if (is.na(r_twosided$estimate)) {
            r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "two.sided")
            #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
          } else if (r_twosided$estimate>0){
            r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "greater")
            #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "greater")
          } else if (r_twosided$estimate<0) {
            r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "less")
            #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "less")
          } else {
            r <- cor.test(as.numeric(test$A),as.numeric(test$B),method = cor_methods,alternative = "two.sided")
            #r <- cor.test(forcor[,i],forcor[,j+i],method = cor_methods,alternative = "two.sided")
          }
          result_tmp <- data.frame(paste0(colnames(forcor1)[j]," vs ",colnames(forcor2)[j]),r$estimate,r$p.value,r$method,r$conf.int[1],r$conf.int[2],r_twosided$p.value)
        } else if(cor_methods=="spearman") {
          r_twosided <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
          #r_twosided <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
          if (is.na(r_twosided$estimate)) {
            r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
            #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
          } else if(r_twosided$estimate>0){
            r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "greater", approximation = "exact")
            #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "greater", approximation = "exact")
          } else if (r_twosided$estimate<0) {
            r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "less", approximation = "exact")
            #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "less", approximation = "exact")
          } else {
            r <- spearman.test(as.numeric(test$A),as.numeric(test$B),alternative = "two.sided", approximation = "exact")
            #r <- spearman.test(forcor[,i],forcor[,j+i],alternative = "two.sided", approximation = "exact")
          }
          result_tmp <- data.frame(paste0(colnames(forcor1)[j]," vs ",colnames(forcor2)[j]),r$estimate,r$p.value,r$method,spearman_CI(as.numeric(test$A),as.numeric(test$B))[1],spearman_CI(as.numeric(test$A),as.numeric(test$B))[2],r_twosided$p.value)
          #result_tmp <- data.frame(paste0(colnames(forcor)[i]," vs ",colnames(forcor)[j+i]),r$estimate,r$p.value,r$method,spearman_CI(forcor[,i],forcor[,i+j])[1],spearman_CI(forcor[,i],forcor[,i+j])[2],r_twosided$p.value)
        }
        
        colnames(result_tmp) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
        result <- rbind(result,result_tmp)
        
        if(j%%5000==0){
          result_final <- rbind(result_final,result)
          write.csv(result,paste0(output_dir,j/5000,"_Result5000.csv"))
          result <- as.data.frame(matrix(numeric(0),ncol=7))
          colnames(result) <- c("DataNames","R","Pvalue","Method","Lower Confidence Interval","Upper Confidence Interval","Pvalue_twosided")
          j=j+1
        } else if(j==(ncol(forcor1)-1)){
          result_final <- rbind(result_final,result)
          write.csv(result,paste0(output_dir,j%/%5000+1,"_Result",j%%5000,".csv"))
          j=j+1
        } else {
          j=j+1
        }
        
        if(is.na(r$p.value)){next}
        if(r$p.value<=pvalue_cutoff){
          p1 <-
            ggplot(test,aes(x=ID,y=0,fill= as.numeric(as.factor(test$A))))+
            geom_tile()+
            scale_fill_gradient2(low="black",mid = "white", high = "red", midpoint = length(test$A)/2)+
            #scale_fill_manual(values=c("-"="white","early"="sky blue","late"="red"))+
            theme_bw()+
            theme(#legend.position="right",
              plot.margin = unit(x=c(0,0,0,40),units="pt"),
              legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
              legend.title = element_blank(),
              #legend.text= element_blank(),
              plot.title = element_blank(),
              axis.ticks = element_blank(),
              #axis.text.x = element_blank(),
              axis.line.y = element_blank(),
              axis.text.x= element_blank(),
              #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
              axis.text.y = element_text(face="bold", color="black", size=0, angle = 90,hjust=1,vjust = 0.5),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
          
          p2 <-
            ggplot(test,aes(x=ID,y=0,fill= as.numeric(as.factor(test$B))))+
            geom_tile()+
            scale_fill_gradient2(low="black",mid = "white", high = "red", midpoint = length(test$A)/2)+
            theme_bw()+
            theme(
              plot.margin = unit(x=c(0,0,0,40),units="pt"),
              legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
              legend.title = element_blank(),
              plot.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line.y = element_blank(),
              #axis.text.x= element_blank(),
              axis.text.x = element_text( color="black", size=15, angle = 45,hjust = 1,vjust = 1),
              axis.text.y = element_text(face="bold", color="black", size=0, angle = 90,hjust=1,vjust = 0.5),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
          
          p_final <- ggarrange(p1,p2,
                               ncol = 1,nrow = 2,align = "v",heights = c(1,2.2),legend = "none")
          plot_name <- gsub("/","-",colnames(forcor1)[j-1])
          ggsave(p_final, path = output_dir, filename = paste0(r$estimate,"_",plot_name,".pdf"), device = "pdf",width = 8.14,height = 2.90)
        }
        
      }
      
      pb$tick()
      Sys.sleep(1 / 100)
    }
    write.csv(result_final,paste0(output_dir,"Paired_result_final.csv"))
  }
  
  
}

Group_color <- c(CRC=alpha("#FCB514",alpha = 1),STAD=alpha("red",alpha = 1),HD=alpha("blue",alpha = 1))
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics")
#Figure 1
{
#cohort summary
{
  ##multi-omics sample shaozhen
  setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/01.Sample summary/Figure3/")   
  
  #sunburst
  p <- read.csv("Figure3_forplot.csv",header = F, stringsAsFactors = FALSE)
  c <- read.csv("Color.csv",header = F,stringsAsFactors=FALSE)
  c[1,1] <- "NA"
  sunburst(p,
           count=TRUE,
           sortFunction = htmlwidgets::JS(
             "function(a,b) {
                       // sort by count descending
                       // unlike the other example using data.name, value is at the top level of the object
                       return b.value - a.value}" ),
           legend=list(w=240,h=30,r=10,s=5),
           #sortFunction = htmlwidgets::JS("Health","Colorectum","Liver","Stomach","Lung","Thyroid","Esophagus"),
           legendOrder = list("NC","CRC","STAD","Female","Male","Stage I","Stage II","Stage III","Stage IV","Unknown","long RNA","DNA methylation","WGS","small RNA"),
           colors=list(range = c$V2, domain = c$V1)
  )  
  
  #Age_Stage barplot
  clinical <- read.csv("./Figure3_subtypes_forplot.csv",header = TRUE)
  View(clinical)
  clinical[clinical==""]<-NA
  
  
  clinical$Stage.summary <- factor(clinical$Stage.summary,levels = c("missing","Stage I","Stage II","Stage III","Stage IV"))
  ggplot(data=clinical[which(clinical$Type!="NC"),],aes(x=Age_10,y=number,fill=Stage.summary))+ geom_bar(stat = "identity", width=0.9, col='transparent') +
    theme_bw()+
    theme(#legend.position="bottom",
      legend.position="right",
      panel.grid.major.x = element_blank(),
      legend.title = element_text(face="bold", color="transparent",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(color="black", size=20,angle = 0,hjust = 0.5),
      axis.text.y = element_text(color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    ylab("Number of patients")+xlab("Ages(years)")+
    scale_y_continuous(expand = c(0,0),limits = c(0,35))+
    scale_fill_manual(values=c("white","#DCDCDC","#A3A3A3","#4A4A4A","#000000"))
  
  
  
  #subtype
  clinical <- read.csv("./Figure3_simple_subtype_forplot.csv",header = TRUE,stringsAsFactors = FALSE)
  color <- read.csv("Color_for_subtype.csv",header = FALSE)
  clinical[clinical==""]<-NA
  clinical[is.na(clinical)] <- "No biopsy"
  
  color$V1 <- factor(color$V1,levels = color$V1)
  
  clinical$index <- factor(clinical$index,levels = c("CRC position","MMR","STAD position","HER2"))
  clinical$Subtype <- factor(clinical$Subtype,levels = color$V1)
  ggplot(data=clinical,aes(x=number,y=number,fill=clinical$Subtype))+ geom_bar(stat = "identity", width=10, col='transparent') + facet_grid(~index)+
    #geom_text(aes(label=Subtype), position = "stack",vjust=1)+
    theme_bw()+
    theme(
      legend.position="none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.title = element_text(face="bold", color="transparent",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.text.x = element_text(color="black", size=20,angle = 0,hjust = 0.5),
      axis.text.y = element_text(color="black", size=20),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24),
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=10))+
    ylab("Number of patients")+xlab("")+scale_y_reverse()+scale_x_discrete(position = "top")+
    scale_fill_manual(values=as.character(color$V2))
  
}

#omics correlation
{
  TPM_Expression <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/PCAWG_biomarker_candidate/Expression/Expression_DNA_RNA_paired_with_MT.txt",sep="\t",header = T, row.names = 1)
  TPM_CNV <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/PCAWG_biomarker_candidate/DNA-CNV-CPM/CNV-CPM_DNA_RNA_paired.txt",sep="\t",header = T, row.names = 1)
  TPM_Methylation <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/PCAWG_biomarker_candidate/MeDIP/MeDIP_DNA_RNA_paired.txt",sep="\t",header = T, row.names = 1)
  
  { 
    colnames(TPM_Expression) <- gsub(".","-",fixed=TRUE,colnames(TPM_Expression))
    colnames(TPM_Expression) <- gsub("-pico","-expression",fixed=TRUE,colnames(TPM_Expression))
    rownames(TPM_Expression) <- as.character(lapply(strsplit(rownames(TPM_Expression),"|",fixed = TRUE), function(x) x[1]))
    
    colnames(TPM_CNV) <- gsub(".","-",fixed=TRUE,colnames(TPM_CNV))
    colnames(TPM_CNV) <- gsub("-wgs","-cnv",fixed=TRUE,colnames(TPM_CNV))
    rownames(TPM_CNV) <- as.character(lapply(strsplit(rownames(TPM_CNV),"|",fixed = TRUE), function(x) x[1]))
    
    colnames(TPM_Methylation) <- gsub(".","-",fixed=TRUE,colnames(TPM_Methylation))
    colnames(TPM_Methylation) <- gsub("-me","-methylation",fixed=TRUE,colnames(TPM_Methylation))
    rownames(TPM_Methylation) <- as.character(lapply(strsplit(rownames(TPM_Methylation),"|",fixed = TRUE), function(x) x[1]))
    
    omics_paired <- cbind(TPM_CNV,TPM_Methylation)
    omics_paired <- inner_join(rownames_to_column(omics_paired),rownames_to_column(TPM_Expression),by=c("rowname"="rowname"))
    rownames(omics_paired) <- omics_paired$rowname
    omics_paired <- omics_paired[,-which(colnames(omics_paired)=="rowname")]
    omics_paired_ensembl <- omics_paired[grep("ENSG",rownames(omics_paired)),]
    
    # cor and pvalue bubble plot
    library(corrplot)
    library(RColorBrewer)
    
    M <-cor(omics_paired_ensembl)
    library("Hmisc")
    col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                   "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                   "#4393C3", "#2166AC", "#053061")))
    # Insignificant correlation are crossed
    res2 <- rcorr(as.matrix(omics_paired_ensembl),type="pearson")
    corrplot(corr = res2$r,col = col2(200),tl.col="black",type="lower", order="original",tl.pos = "ld",tl.cex=0.7,tl.srt = 45,
             p.mat = res2$P, sig.level = 0.05,insig = "blank")  # insig = "blank"
    
    #write.csv(res2$r,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired.csv")
    #average_r <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average.csv",header = TRUE, row.names = 1)
    write.csv(res2$P,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired-20211007-p.csv")
    write.csv(res2$r,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired-20211007.csv")
    
    all_r <- res2$r
    all_p <- res2$P
    corrplot(corr = all_r[grep("cnv|methylation",rownames(all_r)),grep("cnv|methylation",colnames(all_r))],col = col2(200),tl.col="white",type="lower", order="original",tl.pos = "ld",#tl.cex=0.7,tl.srt = 45,
             p.mat = all_p[grep("cnv|methylation",rownames(all_p)),grep("cnv|methylation",colnames(all_p))], sig.level = 0.01,insig = "blank",addgrid.col = NA)  # insig = "blank"
    
    corrplot(corr = all_r[grep("cnv|expression",rownames(all_r)),grep("cnv|expression",colnames(all_r))],col = col2(200),tl.col="white",type="lower", order="original",tl.pos = "ld",#tl.cex=0.7,tl.srt = 45,
             p.mat = all_p[grep("cnv|expression",rownames(all_p)),grep("cnv|expression",colnames(all_p))], sig.level = 0.01,insig = "blank",addgrid.col = NA)  # insig = "blank"
    
    corrplot(corr = all_r[grep("methylation|expression",rownames(all_r)),grep("methylation|expression",colnames(all_r))],col = col2(200),tl.col="white",type="lower", order="original",tl.pos = "ld",#tl.cex=0.7,tl.srt = 45,
             p.mat = all_p[grep("methylation|expression",rownames(all_p)),grep("methylation|expression",colnames(all_p))], sig.level = 0.01,insig = "blank",addgrid.col = NA)  # insig = "blank"
    
    CNV_r <- all_r[grep("cnv",rownames(all_r)),grep("cnv",colnames(all_r))]
    average_CNV_r <- (sum(CNV_r)-ncol(CNV_r))/(ncol(CNV_r)*nrow(CNV_r)-ncol(CNV_r))
    
    Methy_r <- all_r[grep("methylation",rownames(all_r)),grep("methylation",colnames(all_r))]
    average_Methy_r <- (sum(Methy_r)-ncol(Methy_r))/(ncol(Methy_r)*nrow(Methy_r)-ncol(Methy_r))
    
    Expression_r <- all_r[grep("expression",rownames(all_r)),grep("expression",colnames(all_r))]
    average_Expression_r <- (sum(Expression_r)-ncol(Expression_r))/(ncol(Expression_r)*nrow(Expression_r)-ncol(Expression_r))
    
    Expression_CNV_r <- all_r[grep("expression",rownames(all_r)),grep("cnv",colnames(all_r))]
    average_Expression_CNV_r <- (sum(Expression_CNV_r))/(ncol(Expression_CNV_r)*nrow(Expression_CNV_r))
    
    Expression_Methy_r <- all_r[grep("expression",rownames(all_r)),grep("methylation",colnames(all_r))]
    average_Expression_Methy_r <- (sum(Expression_Methy_r))/(ncol(Expression_Methy_r)*nrow(Expression_Methy_r))
    
    CNV_Methy_r <- all_r[grep("cnv",rownames(all_r)),grep("methylation",colnames(all_r))]
    average_CNV_Methy_r <- (sum(CNV_Methy_r))/(ncol(CNV_Methy_r)*nrow(CNV_Methy_r))
    
    average_r <- as.data.frame(matrix(numeric(0),ncol=3,nrow=3))
    
    colnames(average_r) <- c("CNV","Methylation","Expression")
    rownames(average_r) <- c("CNV","Methylation","Expression")
    
    average_r[which(rownames(average_r)=="CNV"),which(colnames(average_r)=="CNV")] <- average_CNV_r
    average_r[which(rownames(average_r)=="Methylation"),which(colnames(average_r)=="Methylation")] <- average_Methy_r
    average_r[which(rownames(average_r)=="Expression"),which(colnames(average_r)=="Expression")] <- average_Expression_r
    average_r[which(rownames(average_r)=="Expression"),which(colnames(average_r)=="CNV")] <- average_Expression_CNV_r
    average_r[which(rownames(average_r)=="CNV"),which(colnames(average_r)=="Expression")] <- average_Expression_CNV_r
    average_r[which(rownames(average_r)=="Expression"),which(colnames(average_r)=="Methylation")] <- average_Expression_Methy_r
    average_r[which(rownames(average_r)=="Methylation"),which(colnames(average_r)=="Expression")] <- average_Expression_Methy_r
    average_r[which(rownames(average_r)=="CNV"),which(colnames(average_r)=="Methylation")] <- average_CNV_Methy_r
    average_r[which(rownames(average_r)=="Methylation"),which(colnames(average_r)=="CNV")] <- average_CNV_Methy_r
    
    
    write.csv(average_r,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average-20211007.csv")
    corrplot(corr = as.matrix(average_r),method = c("color"),col = col2(200),outline = TRUE,tl.col="black",type="lower", order="original",tl.pos = "l",tl.cex=1.5,cl.cex = 1, number.cex = as.matrix(average_r))
    
    
    {
      all_r <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired-20211007-p.csv",header = TRUE,row.names = 1)
      colnames(all_r) <- gsub(".","-",fixed = TRUE,colnames(all_r))
      Expression_r <- all_r[grep("expression",rownames(all_r)),grep("expression",colnames(all_r))]
      Methy_r <- all_r[grep("methylation",rownames(all_r)),grep("methylation",colnames(all_r))]
      CNV_r <- all_r[grep("cnv",rownames(all_r)),grep("cnv",colnames(all_r))]
      Expression_Methy_r <- all_r[grep("expression",rownames(all_r)),grep("methylation",colnames(all_r))]
      Expression_CNV_r <- all_r[grep("expression",rownames(all_r)),grep("cnv",colnames(all_r))]
      Methy_CNV_r <- all_r[grep("methylation",rownames(all_r)),grep("cnv",colnames(all_r))]
      
      ######Expression triangle
      Expression_CRC_r <- Expression_r[grep("CRC",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_r <- (sum(Expression_CRC_r))/(ncol(Expression_CRC_r)*nrow(Expression_CRC_r))
      
      Expression_STAD_r <- Expression_r[grep("STAD",rownames(Expression_r)),grep("STAD",colnames(Expression_r))]
      average_Expression_STAD_r <- (sum(Expression_STAD_r))/(ncol(Expression_STAD_r)*nrow(Expression_STAD_r))
      
      Expression_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("NC",colnames(Expression_r))]
      average_Expression_NC_r <- (sum(Expression_NC_r))/(ncol(Expression_NC_r)*nrow(Expression_NC_r))
      
      Expression_CRC_STAD_r <- Expression_r[grep("STAD",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_STAD_r <- (sum(Expression_CRC_STAD_r))/(ncol(Expression_CRC_STAD_r)*nrow(Expression_CRC_STAD_r))
      
      Expression_CRC_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_NC_r <- (sum(Expression_CRC_NC_r))/(ncol(Expression_CRC_NC_r)*nrow(Expression_CRC_NC_r))
      
      Expression_STAD_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("STAD",colnames(Expression_r))]
      average_Expression_STAD_NC_r <- (sum(Expression_STAD_NC_r))/(ncol(Expression_STAD_NC_r)*nrow(Expression_STAD_NC_r))
      
      ########Methylation triangle
      Methy_CRC_r <- Methy_r[grep("CRC",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_r <- (sum(Methy_CRC_r))/(ncol(Methy_CRC_r)*nrow(Methy_CRC_r))
      
      Methy_STAD_r <- Methy_r[grep("STAD",rownames(Methy_r)),grep("STAD",colnames(Methy_r))]
      average_Methy_STAD_r <- (sum(Methy_STAD_r))/(ncol(Methy_STAD_r)*nrow(Methy_STAD_r))
      
      Methy_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("NC",colnames(Methy_r))]
      average_Methy_NC_r <- (sum(Methy_NC_r))/(ncol(Methy_NC_r)*nrow(Methy_NC_r))
      
      Methy_CRC_STAD_r <- Methy_r[grep("STAD",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_STAD_r <- (sum(Methy_CRC_STAD_r))/(ncol(Methy_CRC_STAD_r)*nrow(Methy_CRC_STAD_r))
      
      Methy_CRC_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_NC_r <- (sum(Methy_CRC_NC_r))/(ncol(Methy_CRC_NC_r)*nrow(Methy_CRC_NC_r))
      
      Methy_STAD_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("STAD",colnames(Methy_r))]
      average_Methy_STAD_NC_r <- (sum(Methy_STAD_NC_r))/(ncol(Methy_STAD_NC_r)*nrow(Methy_STAD_NC_r))
      
      ########CNV triangle
      CNV_CRC_r <- CNV_r[grep("CRC",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_r <- (sum(CNV_CRC_r))/(ncol(CNV_CRC_r)*nrow(CNV_CRC_r))
      
      CNV_STAD_r <- CNV_r[grep("STAD",rownames(CNV_r)),grep("STAD",colnames(CNV_r))]
      average_CNV_STAD_r <- (sum(CNV_STAD_r))/(ncol(CNV_STAD_r)*nrow(CNV_STAD_r))
      
      CNV_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("NC",colnames(CNV_r))]
      average_CNV_NC_r <- (sum(CNV_NC_r))/(ncol(CNV_NC_r)*nrow(CNV_NC_r))
      
      CNV_CRC_STAD_r <- CNV_r[grep("STAD",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_STAD_r <- (sum(CNV_CRC_STAD_r))/(ncol(CNV_CRC_STAD_r)*nrow(CNV_CRC_STAD_r))
      
      CNV_CRC_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_NC_r <- (sum(CNV_CRC_NC_r))/(ncol(CNV_CRC_NC_r)*nrow(CNV_CRC_NC_r))
      
      CNV_STAD_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("STAD",colnames(CNV_r))]
      average_CNV_STAD_NC_r <- (sum(CNV_STAD_NC_r))/(ncol(CNV_STAD_NC_r)*nrow(CNV_STAD_NC_r))
      
      
      ##### RNA vs. Methylation rectangle
      Expression_CRC_Methy_CRC_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_CRC_r <- (sum(Expression_CRC_Methy_CRC_r))/(ncol(Expression_CRC_Methy_CRC_r)*nrow(Expression_CRC_Methy_CRC_r))
      
      Expression_STAD_Methy_CRC_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_CRC_r <- (sum(Expression_STAD_Methy_CRC_r))/(ncol(Expression_STAD_Methy_CRC_r)*nrow(Expression_STAD_Methy_CRC_r))
      
      Expression_NC_Methy_CRC_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_CRC_r <- (sum(Expression_NC_Methy_CRC_r))/(ncol(Expression_NC_Methy_CRC_r)*nrow(Expression_NC_Methy_CRC_r))
      
      Expression_CRC_Methy_STAD_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_STAD_r <- (sum(Expression_CRC_Methy_STAD_r))/(ncol(Expression_CRC_Methy_STAD_r)*nrow(Expression_CRC_Methy_STAD_r))
      
      Expression_STAD_Methy_STAD_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_STAD_r <- (sum(Expression_STAD_Methy_STAD_r))/(ncol(Expression_STAD_Methy_STAD_r)*nrow(Expression_STAD_Methy_STAD_r))
      
      Expression_NC_Methy_STAD_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_STAD_r <- (sum(Expression_NC_Methy_STAD_r))/(ncol(Expression_NC_Methy_STAD_r)*nrow(Expression_NC_Methy_STAD_r))
      
      Expression_CRC_Methy_NC_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_NC_r <- (sum(Expression_CRC_Methy_NC_r))/(ncol(Expression_CRC_Methy_NC_r)*nrow(Expression_CRC_Methy_NC_r))
      
      Expression_STAD_Methy_NC_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_NC_r <- (sum(Expression_STAD_Methy_NC_r))/(ncol(Expression_STAD_Methy_NC_r)*nrow(Expression_STAD_Methy_NC_r))
      
      Expression_NC_Methy_NC_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_NC_r <- (sum(Expression_NC_Methy_NC_r))/(ncol(Expression_NC_Methy_NC_r)*nrow(Expression_NC_Methy_NC_r))
      ####### RNA vs. CNV rectangle
      Expression_CRC_CNV_CRC_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_CRC_r <- (sum(Expression_CRC_CNV_CRC_r))/(ncol(Expression_CRC_CNV_CRC_r)*nrow(Expression_CRC_CNV_CRC_r))
      
      Expression_STAD_CNV_CRC_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_CRC_r <- (sum(Expression_STAD_CNV_CRC_r))/(ncol(Expression_STAD_CNV_CRC_r)*nrow(Expression_STAD_CNV_CRC_r))
      
      Expression_NC_CNV_CRC_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_CRC_r <- (sum(Expression_NC_CNV_CRC_r))/(ncol(Expression_NC_CNV_CRC_r)*nrow(Expression_NC_CNV_CRC_r))
      
      Expression_CRC_CNV_STAD_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_STAD_r <- (sum(Expression_CRC_CNV_STAD_r))/(ncol(Expression_CRC_CNV_STAD_r)*nrow(Expression_CRC_CNV_STAD_r))
      
      Expression_STAD_CNV_STAD_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_STAD_r <- (sum(Expression_STAD_CNV_STAD_r))/(ncol(Expression_STAD_CNV_STAD_r)*nrow(Expression_STAD_CNV_STAD_r))
      
      Expression_NC_CNV_STAD_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_STAD_r <- (sum(Expression_NC_CNV_STAD_r))/(ncol(Expression_NC_CNV_STAD_r)*nrow(Expression_NC_CNV_STAD_r))
      
      Expression_CRC_CNV_NC_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_NC_r <- (sum(Expression_CRC_CNV_NC_r))/(ncol(Expression_CRC_CNV_NC_r)*nrow(Expression_CRC_CNV_NC_r))
      
      Expression_STAD_CNV_NC_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_NC_r <- (sum(Expression_STAD_CNV_NC_r))/(ncol(Expression_STAD_CNV_NC_r)*nrow(Expression_STAD_CNV_NC_r))
      
      Expression_NC_CNV_NC_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_NC_r <- (sum(Expression_NC_CNV_NC_r))/(ncol(Expression_NC_CNV_NC_r)*nrow(Expression_NC_CNV_NC_r))
      
      ####### Methylation vs CNV rectangle
      Methy_CRC_CNV_CRC_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_CRC_r <- (sum(Methy_CRC_CNV_CRC_r))/(ncol(Methy_CRC_CNV_CRC_r)*nrow(Methy_CRC_CNV_CRC_r))
      
      Methy_STAD_CNV_CRC_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_CRC_r <- (sum(Methy_STAD_CNV_CRC_r))/(ncol(Methy_STAD_CNV_CRC_r)*nrow(Methy_STAD_CNV_CRC_r))
      
      Methy_NC_CNV_CRC_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_CRC_r <- (sum(Methy_NC_CNV_CRC_r))/(ncol(Methy_NC_CNV_CRC_r)*nrow(Methy_NC_CNV_CRC_r))
      
      Methy_CRC_CNV_STAD_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_STAD_r <- (sum(Methy_CRC_CNV_STAD_r))/(ncol(Methy_CRC_CNV_STAD_r)*nrow(Methy_CRC_CNV_STAD_r))
      
      Methy_STAD_CNV_STAD_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_STAD_r <- (sum(Methy_STAD_CNV_STAD_r))/(ncol(Methy_STAD_CNV_STAD_r)*nrow(Methy_STAD_CNV_STAD_r))
      
      Methy_NC_CNV_STAD_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_STAD_r <- (sum(Methy_NC_CNV_STAD_r))/(ncol(Methy_NC_CNV_STAD_r)*nrow(Methy_NC_CNV_STAD_r))
      
      Methy_CRC_CNV_NC_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_NC_r <- (sum(Methy_CRC_CNV_NC_r))/(ncol(Methy_CRC_CNV_NC_r)*nrow(Methy_CRC_CNV_NC_r))
      
      Methy_STAD_CNV_NC_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_NC_r <- (sum(Methy_STAD_CNV_NC_r))/(ncol(Methy_STAD_CNV_NC_r)*nrow(Methy_STAD_CNV_NC_r))
      
      Methy_NC_CNV_NC_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_NC_r <- (sum(Methy_NC_CNV_NC_r))/(ncol(Methy_NC_CNV_NC_r)*nrow(Methy_NC_CNV_NC_r))
      
      average_r <- as.data.frame(matrix(numeric(0),ncol=9,nrow=9))
      
      colnames(average_r) <- c("Expression-CRC","Expression-STAD","Expression-NC","Methylation-CRC","Methylation-STAD","Methylation-NC","CNV-CRC","CNV-STAD","CNV-NC")
      rownames(average_r) <- c("Expression-CRC","Expression-STAD","Expression-NC","Methylation-CRC","Methylation-STAD","Methylation-NC","CNV-CRC","CNV-STAD","CNV-NC")
      
      #Expression triangle
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_CRC_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_CRC_NC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_NC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_STAD_NC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_NC_r
      #Methylation triangle
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_CRC_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_STAD_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_CRC_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_NC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Methy_STAD_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_NC_r
      #copynumber triangle
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-STAD")] <- CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-NC")] <- CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-STAD")] <- CNV_CRC_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_STAD_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-NC")] <- CNV_CRC_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_NC_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-NC")] <- CNV_STAD_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-STAD")] <- CNV_STAD_NC_r
      #Expression vs. Methylation
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_CRC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_CRC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Expression_CRC_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_NC_r
      
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_STAD_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_STAD_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Expression_STAD_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_NC_r
      
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_NC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_NC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-NC")] <- Expression_NC_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_NC_r
      
      #Expression vs. CNV
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-CRC")] <- Expression_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-STAD")] <- Expression_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-NC")] <- Expression_CRC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-CRC")] <- Expression_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-STAD")] <- Expression_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-NC")] <- Expression_STAD_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-CRC")] <- Expression_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-STAD")] <- Expression_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-NC")] <- Expression_NC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_NC_r
      
      #Methylation vs. CNV
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-CRC")] <- Methy_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-STAD")] <- Methy_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-NC")] <- Methy_CRC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-CRC")] <- Methy_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-STAD")] <- Methy_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-NC")] <- Methy_STAD_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-CRC")] <- Methy_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-STAD")] <- Methy_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-NC")] <- Methy_NC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_NC_r
      
      write.csv(average_r,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average-20211026-p.csv")
      #library(corrplot)
      #corrplot(corr = as.matrix(average_r),
      #         method = c("circle"),col = col2(200),
      #         outline = TRUE,tl.col="black",
      #         type="full", 
      #         order="original",
      #         tl.pos = "l",tl.cex=1.5,
      #         cl.cex = 1,cl.lim= c(-0.1,1),cl.length = 12,
      #         number.cex = as.matrix(average_r),
      #         #order = "hclust",
      #         addrect = 3)
    }
    {
      all_r <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired-20211007.csv",header = TRUE,row.names = 1)
      colnames(all_r) <- gsub(".","-",fixed = TRUE,colnames(all_r))
      Expression_r <- all_r[grep("expression",rownames(all_r)),grep("expression",colnames(all_r))]
      Methy_r <- all_r[grep("methylation",rownames(all_r)),grep("methylation",colnames(all_r))]
      CNV_r <- all_r[grep("cnv",rownames(all_r)),grep("cnv",colnames(all_r))]
      Expression_Methy_r <- all_r[grep("expression",rownames(all_r)),grep("methylation",colnames(all_r))]
      Expression_CNV_r <- all_r[grep("expression",rownames(all_r)),grep("cnv",colnames(all_r))]
      Methy_CNV_r <- all_r[grep("methylation",rownames(all_r)),grep("cnv",colnames(all_r))]
      
      ######Expression triangle
      Expression_CRC_r <- Expression_r[grep("CRC",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_r <- (sum(Expression_CRC_r))/(ncol(Expression_CRC_r)*nrow(Expression_CRC_r))
      
      Expression_STAD_r <- Expression_r[grep("STAD",rownames(Expression_r)),grep("STAD",colnames(Expression_r))]
      average_Expression_STAD_r <- (sum(Expression_STAD_r))/(ncol(Expression_STAD_r)*nrow(Expression_STAD_r))
      
      Expression_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("NC",colnames(Expression_r))]
      average_Expression_NC_r <- (sum(Expression_NC_r))/(ncol(Expression_NC_r)*nrow(Expression_NC_r))
      
      Expression_CRC_STAD_r <- Expression_r[grep("STAD",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_STAD_r <- (sum(Expression_CRC_STAD_r))/(ncol(Expression_CRC_STAD_r)*nrow(Expression_CRC_STAD_r))
      
      Expression_CRC_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("CRC",colnames(Expression_r))]
      average_Expression_CRC_NC_r <- (sum(Expression_CRC_NC_r))/(ncol(Expression_CRC_NC_r)*nrow(Expression_CRC_NC_r))
      
      Expression_STAD_NC_r <- Expression_r[grep("NC",rownames(Expression_r)),grep("STAD",colnames(Expression_r))]
      average_Expression_STAD_NC_r <- (sum(Expression_STAD_NC_r))/(ncol(Expression_STAD_NC_r)*nrow(Expression_STAD_NC_r))
      
      ########Methylation triangle
      Methy_CRC_r <- Methy_r[grep("CRC",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_r <- (sum(Methy_CRC_r))/(ncol(Methy_CRC_r)*nrow(Methy_CRC_r))
      
      Methy_STAD_r <- Methy_r[grep("STAD",rownames(Methy_r)),grep("STAD",colnames(Methy_r))]
      average_Methy_STAD_r <- (sum(Methy_STAD_r))/(ncol(Methy_STAD_r)*nrow(Methy_STAD_r))
      
      Methy_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("NC",colnames(Methy_r))]
      average_Methy_NC_r <- (sum(Methy_NC_r))/(ncol(Methy_NC_r)*nrow(Methy_NC_r))
      
      Methy_CRC_STAD_r <- Methy_r[grep("STAD",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_STAD_r <- (sum(Methy_CRC_STAD_r))/(ncol(Methy_CRC_STAD_r)*nrow(Methy_CRC_STAD_r))
      
      Methy_CRC_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("CRC",colnames(Methy_r))]
      average_Methy_CRC_NC_r <- (sum(Methy_CRC_NC_r))/(ncol(Methy_CRC_NC_r)*nrow(Methy_CRC_NC_r))
      
      Methy_STAD_NC_r <- Methy_r[grep("NC",rownames(Methy_r)),grep("STAD",colnames(Methy_r))]
      average_Methy_STAD_NC_r <- (sum(Methy_STAD_NC_r))/(ncol(Methy_STAD_NC_r)*nrow(Methy_STAD_NC_r))
      
      ########CNV triangle
      CNV_CRC_r <- CNV_r[grep("CRC",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_r <- (sum(CNV_CRC_r))/(ncol(CNV_CRC_r)*nrow(CNV_CRC_r))
      
      CNV_STAD_r <- CNV_r[grep("STAD",rownames(CNV_r)),grep("STAD",colnames(CNV_r))]
      average_CNV_STAD_r <- (sum(CNV_STAD_r))/(ncol(CNV_STAD_r)*nrow(CNV_STAD_r))
      
      CNV_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("NC",colnames(CNV_r))]
      average_CNV_NC_r <- (sum(CNV_NC_r))/(ncol(CNV_NC_r)*nrow(CNV_NC_r))
      
      CNV_CRC_STAD_r <- CNV_r[grep("STAD",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_STAD_r <- (sum(CNV_CRC_STAD_r))/(ncol(CNV_CRC_STAD_r)*nrow(CNV_CRC_STAD_r))
      
      CNV_CRC_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("CRC",colnames(CNV_r))]
      average_CNV_CRC_NC_r <- (sum(CNV_CRC_NC_r))/(ncol(CNV_CRC_NC_r)*nrow(CNV_CRC_NC_r))
      
      CNV_STAD_NC_r <- CNV_r[grep("NC",rownames(CNV_r)),grep("STAD",colnames(CNV_r))]
      average_CNV_STAD_NC_r <- (sum(CNV_STAD_NC_r))/(ncol(CNV_STAD_NC_r)*nrow(CNV_STAD_NC_r))
      
      
      ##### RNA vs. Methylation rectangle
      Expression_CRC_Methy_CRC_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_CRC_r <- (sum(Expression_CRC_Methy_CRC_r))/(ncol(Expression_CRC_Methy_CRC_r)*nrow(Expression_CRC_Methy_CRC_r))
      
      Expression_STAD_Methy_CRC_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_CRC_r <- (sum(Expression_STAD_Methy_CRC_r))/(ncol(Expression_STAD_Methy_CRC_r)*nrow(Expression_STAD_Methy_CRC_r))
      
      Expression_NC_Methy_CRC_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("CRC",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_CRC_r <- (sum(Expression_NC_Methy_CRC_r))/(ncol(Expression_NC_Methy_CRC_r)*nrow(Expression_NC_Methy_CRC_r))
      
      Expression_CRC_Methy_STAD_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_STAD_r <- (sum(Expression_CRC_Methy_STAD_r))/(ncol(Expression_CRC_Methy_STAD_r)*nrow(Expression_CRC_Methy_STAD_r))
      
      Expression_STAD_Methy_STAD_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_STAD_r <- (sum(Expression_STAD_Methy_STAD_r))/(ncol(Expression_STAD_Methy_STAD_r)*nrow(Expression_STAD_Methy_STAD_r))
      
      Expression_NC_Methy_STAD_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("STAD",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_STAD_r <- (sum(Expression_NC_Methy_STAD_r))/(ncol(Expression_NC_Methy_STAD_r)*nrow(Expression_NC_Methy_STAD_r))
      
      Expression_CRC_Methy_NC_r <- Expression_Methy_r[grep("CRC",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_CRC_Methy_NC_r <- (sum(Expression_CRC_Methy_NC_r))/(ncol(Expression_CRC_Methy_NC_r)*nrow(Expression_CRC_Methy_NC_r))
      
      Expression_STAD_Methy_NC_r <- Expression_Methy_r[grep("STAD",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_STAD_Methy_NC_r <- (sum(Expression_STAD_Methy_NC_r))/(ncol(Expression_STAD_Methy_NC_r)*nrow(Expression_STAD_Methy_NC_r))
      
      Expression_NC_Methy_NC_r <- Expression_Methy_r[grep("NC",rownames(Expression_Methy_r)),grep("NC",colnames(Expression_Methy_r))]
      average_Expression_NC_Methy_NC_r <- (sum(Expression_NC_Methy_NC_r))/(ncol(Expression_NC_Methy_NC_r)*nrow(Expression_NC_Methy_NC_r))
      ####### RNA vs. CNV rectangle
      Expression_CRC_CNV_CRC_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_CRC_r <- (sum(Expression_CRC_CNV_CRC_r))/(ncol(Expression_CRC_CNV_CRC_r)*nrow(Expression_CRC_CNV_CRC_r))
      
      Expression_STAD_CNV_CRC_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_CRC_r <- (sum(Expression_STAD_CNV_CRC_r))/(ncol(Expression_STAD_CNV_CRC_r)*nrow(Expression_STAD_CNV_CRC_r))
      
      Expression_NC_CNV_CRC_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("CRC",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_CRC_r <- (sum(Expression_NC_CNV_CRC_r))/(ncol(Expression_NC_CNV_CRC_r)*nrow(Expression_NC_CNV_CRC_r))
      
      Expression_CRC_CNV_STAD_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_STAD_r <- (sum(Expression_CRC_CNV_STAD_r))/(ncol(Expression_CRC_CNV_STAD_r)*nrow(Expression_CRC_CNV_STAD_r))
      
      Expression_STAD_CNV_STAD_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_STAD_r <- (sum(Expression_STAD_CNV_STAD_r))/(ncol(Expression_STAD_CNV_STAD_r)*nrow(Expression_STAD_CNV_STAD_r))
      
      Expression_NC_CNV_STAD_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("STAD",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_STAD_r <- (sum(Expression_NC_CNV_STAD_r))/(ncol(Expression_NC_CNV_STAD_r)*nrow(Expression_NC_CNV_STAD_r))
      
      Expression_CRC_CNV_NC_r <- Expression_CNV_r[grep("CRC",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_CRC_CNV_NC_r <- (sum(Expression_CRC_CNV_NC_r))/(ncol(Expression_CRC_CNV_NC_r)*nrow(Expression_CRC_CNV_NC_r))
      
      Expression_STAD_CNV_NC_r <- Expression_CNV_r[grep("STAD",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_STAD_CNV_NC_r <- (sum(Expression_STAD_CNV_NC_r))/(ncol(Expression_STAD_CNV_NC_r)*nrow(Expression_STAD_CNV_NC_r))
      
      Expression_NC_CNV_NC_r <- Expression_CNV_r[grep("NC",rownames(Expression_CNV_r)),grep("NC",colnames(Expression_CNV_r))]
      average_Expression_NC_CNV_NC_r <- (sum(Expression_NC_CNV_NC_r))/(ncol(Expression_NC_CNV_NC_r)*nrow(Expression_NC_CNV_NC_r))
      
      ####### Methylation vs CNV rectangle
      Methy_CRC_CNV_CRC_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_CRC_r <- (sum(Methy_CRC_CNV_CRC_r))/(ncol(Methy_CRC_CNV_CRC_r)*nrow(Methy_CRC_CNV_CRC_r))
      
      Methy_STAD_CNV_CRC_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_CRC_r <- (sum(Methy_STAD_CNV_CRC_r))/(ncol(Methy_STAD_CNV_CRC_r)*nrow(Methy_STAD_CNV_CRC_r))
      
      Methy_NC_CNV_CRC_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("CRC",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_CRC_r <- (sum(Methy_NC_CNV_CRC_r))/(ncol(Methy_NC_CNV_CRC_r)*nrow(Methy_NC_CNV_CRC_r))
      
      Methy_CRC_CNV_STAD_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_STAD_r <- (sum(Methy_CRC_CNV_STAD_r))/(ncol(Methy_CRC_CNV_STAD_r)*nrow(Methy_CRC_CNV_STAD_r))
      
      Methy_STAD_CNV_STAD_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_STAD_r <- (sum(Methy_STAD_CNV_STAD_r))/(ncol(Methy_STAD_CNV_STAD_r)*nrow(Methy_STAD_CNV_STAD_r))
      
      Methy_NC_CNV_STAD_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("STAD",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_STAD_r <- (sum(Methy_NC_CNV_STAD_r))/(ncol(Methy_NC_CNV_STAD_r)*nrow(Methy_NC_CNV_STAD_r))
      
      Methy_CRC_CNV_NC_r <- Methy_CNV_r[grep("CRC",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_CRC_CNV_NC_r <- (sum(Methy_CRC_CNV_NC_r))/(ncol(Methy_CRC_CNV_NC_r)*nrow(Methy_CRC_CNV_NC_r))
      
      Methy_STAD_CNV_NC_r <- Methy_CNV_r[grep("STAD",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_STAD_CNV_NC_r <- (sum(Methy_STAD_CNV_NC_r))/(ncol(Methy_STAD_CNV_NC_r)*nrow(Methy_STAD_CNV_NC_r))
      
      Methy_NC_CNV_NC_r <- Methy_CNV_r[grep("NC",rownames(Methy_CNV_r)),grep("NC",colnames(Methy_CNV_r))]
      average_Methy_NC_CNV_NC_r <- (sum(Methy_NC_CNV_NC_r))/(ncol(Methy_NC_CNV_NC_r)*nrow(Methy_NC_CNV_NC_r))
      
      average_r <- as.data.frame(matrix(numeric(0),ncol=9,nrow=9))
      
      colnames(average_r) <- c("Expression-CRC","Expression-STAD","Expression-NC","Methylation-CRC","Methylation-STAD","Methylation-NC","CNV-CRC","CNV-STAD","CNV-NC")
      rownames(average_r) <- c("Expression-CRC","Expression-STAD","Expression-NC","Methylation-CRC","Methylation-STAD","Methylation-NC","CNV-CRC","CNV-STAD","CNV-NC")
      
      #Expression triangle
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_CRC_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_CRC_NC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_NC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_STAD_NC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_NC_r
      #Methylation triangle
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_CRC_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_STAD_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_CRC_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_NC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Methy_STAD_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_NC_r
      #copynumber triangle
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-STAD")] <- CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-NC")] <- CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-STAD")] <- CNV_CRC_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_STAD_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="CNV-NC")] <- CNV_CRC_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-CRC")] <- CNV_CRC_NC_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="CNV-NC")] <- CNV_STAD_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="CNV-STAD")] <- CNV_STAD_NC_r
      #Expression vs. Methylation
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_CRC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_CRC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Expression_CRC_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_Methy_NC_r
      
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_STAD_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_STAD_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Expression_STAD_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_Methy_NC_r
      
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Expression_NC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_CRC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Expression_NC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="Methylation-NC")] <- Expression_NC_Methy_NC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_Methy_NC_r
      
      #Expression vs. CNV
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-CRC")] <- Expression_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-STAD")] <- Expression_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-CRC"),which(colnames(average_r)=="CNV-NC")] <- Expression_CRC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-CRC")] <- Expression_CRC_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-CRC")] <- Expression_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-STAD")] <- Expression_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-STAD"),which(colnames(average_r)=="CNV-NC")] <- Expression_STAD_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-STAD")] <- Expression_STAD_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-CRC")] <- Expression_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-STAD")] <- Expression_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Expression-NC"),which(colnames(average_r)=="CNV-NC")] <- Expression_NC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Expression-NC")] <- Expression_NC_CNV_NC_r
      
      #Methylation vs. CNV
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-CRC")] <- Methy_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-STAD")] <- Methy_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-CRC"),which(colnames(average_r)=="CNV-NC")] <- Methy_CRC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-CRC")] <- Methy_CRC_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-CRC")] <- Methy_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-STAD")] <- Methy_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-STAD"),which(colnames(average_r)=="CNV-NC")] <- Methy_STAD_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-STAD")] <- Methy_STAD_CNV_NC_r
      
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-CRC")] <- Methy_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="CNV-CRC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_CRC_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-STAD")] <- Methy_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="CNV-STAD"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_STAD_r
      average_r[which(rownames(average_r)=="Methylation-NC"),which(colnames(average_r)=="CNV-NC")] <- Methy_NC_CNV_NC_r
      average_r[which(rownames(average_r)=="CNV-NC"),which(colnames(average_r)=="Methylation-NC")] <- Methy_NC_CNV_NC_r
      
      write.csv(average_r,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average-20211026.csv")
      #library(corrplot)
      #corrplot(corr = as.matrix(average_r),
      #         method = c("circle"),col = col2(200),
      #         outline = TRUE,tl.col="black",
      #         type="full", 
      #         order="original",
      #         tl.pos = "l",tl.cex=1.5,
      #         cl.cex = 1,cl.lim= c(-0.1,1),cl.length = 12,
      #         number.cex = as.matrix(average_r),
      #         #order = "hclust",
      #         addrect = 3)
    }
    
    r <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average-20211026.csv",header = TRUE,row.names = 1)
    p <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/03.Expression/Omics_correlation/omics_paired_average-20211026-p.csv",header = TRUE, row.names = 1)
    p[is.na(p)] <- 0
    r$omics <- rownames(r)
    r_forplot <- melt(r)
    r_forplot$variable <- gsub(".","-",fixed=TRUE, r_forplot$variable)
    #paste0(r_forplot$variable,r_forplot$omics)
    colnames(r_forplot) <- c("omics1","omics2","pearson")
    
    p$omics <- rownames(p)
    p_forplot <- melt(p)
    p_forplot$variable <- gsub(".","-",fixed=TRUE, p_forplot$variable)
    #paste0(r_forplot$variable,r_forplot$omics)
    colnames(p_forplot) <- c("omics3","omics4","pvalue")
    corplot <- cbind(r_forplot,p_forplot)
    corplot <- corplot[,-which(colnames(corplot) == "omics3")]
    corplot <- corplot[,-which(colnames(corplot) == "omics4")]
    
    corplot$omics1 <- factor(corplot$omics1,c("Expression-STAD","Expression-CRC","Expression-NC",
                                              "Methylation-STAD","Methylation-CRC","Methylation-NC",
                                              "CNV-STAD","CNV-CRC","CNV-NC"))
    corplot$omics2 <- factor(corplot$omics2,rev(c("Expression-STAD","Expression-CRC","Expression-NC",
                                                  "Methylation-STAD","Methylation-CRC","Methylation-NC",
                                                  "CNV-STAD","CNV-CRC","CNV-NC")))
    ggplot(corplot,aes(x=omics1,y=omics2))+
      geom_point(aes(size=as.numeric(pearson),color=as.numeric(pvalue)))+
      geom_vline(xintercept = c(3.5,6.5))+
      geom_vline(xintercept = c(1.5,2.5,4.5,5.5,7.5,8.5),colour = "grey",size = 0.1)+
      geom_hline(yintercept = c(3.5,6.5))+
      geom_hline(yintercept = c(1.5,2.5,4.5,5.5,7.5,8.5),colour = "grey",size = 0.1)+
      scale_colour_gradient2(low="#B6202E",high="#549EC9",mid="white", midpoint = 0.5)+
      #geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
      #geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
      labs(color=expression(PValue),
           size="Pearson's correlation",
           x="")+
      theme_bw()+
      theme(panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
            axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black",angle = 45,vjust=1,hjust=1),
            axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
            axis.title.y = element_blank())
  }
}
}

#Figure 2
{
  #Overview of multiomics alterations
  {
    
    #ann_plot <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/oncoplot/RNA_sample_annotation.csv", header = TRUE,)
    #Summary mutated genes
    {
      library(data.table)
      setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/gene_mutation_burden")
      input_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/gene_mutation_burden"
      output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/gene_mutation_burden"
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
      
      rownames(mutation_burden_matrix) <- mutation_burden_matrix$`Mutated_Gene`
      mutation_burden_matrix <- mutation_burden_matrix[,-which(colnames(mutation_burden_matrix)=="Mutated_Gene")]
      nonsense_row <- which(rownames(mutation_burden_matrix)=="-|-")
      if(length(nonsense_row)!=0){
        mutation_burden_matrix <- mutation_burden_matrix[-which(rownames(mutation_burden_matrix)=="-|-"),]
      } else {
        mutation_burden_matrix <- mutation_burden_matrix
      }
      mutation_burden_matrix <- as.data.frame(t(mutation_burden_matrix))
      
      mutation_burden_matrix$Case.ID <- rownames(mutation_burden_matrix)
      rownames(mutation_burden_matrix) <- seq(from=1,to=nrow(mutation_burden_matrix))
    }
    
    ##RNA mutation gene examples
    {
      mat_plot <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/SNP/Mutation_burden_class_forplot.txt",header = TRUE, row.names = 1,sep = "\t",check.names = FALSE)
      mat_plot[] <- lapply(mat_plot, as.character) 
      mat_plot[is.na(mat_plot)] <- ""
      
      ann <- ann_plot[grep("STAD|CRC",ann_plot$Samples),]
      mat <- mat_plot[grep("ENSG00000142208|ENSG00000134982|ENSG00000175054|ENSG00000103126|ENSG00000168646|ENSG00000166710|ENSG00000186174|ENSG00000157764|ENSG00000039068|ENSG00000168036|ENSG00000257923|ENSG00000104408|ENSG00000100393|ENSG00000141736|ENSG00000065361|ENSG00000178568|ENSG00000109670|ENSG00000066468|ENSG00000183454|ENSG00000100644|ENSG00000133703|ENSG00000121879|ENSG00000145675|ENSG00000101213|ENSG00000163629|ENSG00000196090|ENSG00000112531|ENSG00000164754|ENSG00000067560|ENSG00000175387|ENSG00000166949|ENSG00000141646|ENSG00000177565|ENSG00000148737|ENSG00000163513|ENSG00000141510|ENSG00000104517|ENSG00000140836",rownames(mat_plot)),]
      mat <- mat[,colnames(mat) %in% ann$Samples]
      
      ann <- ann[,-which(colnames(ann)=="Samples")]
      colours_ann <- list(
        'Gender'=c('Male'='#5CACEE','Female'='#FF0000'),
        'Stage'=c('Stage I'='#DCDCDC','Stage II'='#A3A3A3','Stage III'='#4A4A4A','Stage IV'='#000000'),
        'Location'=c('No biopsy'='#FFFFFF','Others'='#DCDCDC',
                     'Rectum'='#FFF8DC','Sigmoid colon'='#EEDD82','Descending colon'='#FEE5AC','Transverse colon'='#EDCB62','Ascending colon'='#FCB514',
                     'Fundus'='#FEF1B5','Body'='#E5BC3B','Antrum'='#CD9B1D'),
        'MMR'=c('No biopsy'='#FFFFFF','p'='#EDCB62','d'='#FCB514'),
        'HER2'=c('No biopsy'='#FFFFFF','-'='#E5BC3B','+'='#CD9B1D'),
        'Type'=c('Colorectal cancer'='#FCB514','Stomach cancer'='#CD9B1D')
      )
      
      colAnn <- HeatmapAnnotation(col = colours_ann,
                                  #cbar = anno_oncoprint_barplot(),
                                  df = ann,which = 'col')
      
      
      rownames(mat) <- lapply(strsplit(rownames(mat),"|",fixed = TRUE),function(x) x[2])
      mat[is.na(mat)] <- ""
      #mat <- mat_plot[order,]
      
      #deep sequencing
      #mat <- mat[,grep("CRC.PKU.32.pico|NC.PKU.10.pico|NC.PKU.14.pico|CRC.PKU.29.pico|NC.PKU.12.pico|CRC.PKU.30.pico|CRC.PKU.27.pico|CRC.PKU.28.pico|NC.PKU.9.pico|CRC.PKU.31.pico",colnames(mat))]
      #paired sample
      
      col <- c("SNV" = "#ff7f00", "insertion" = "#984ea3", "deletion" = "#4daf4a", "sequence_alteration" = "#003f87","indel"="#8e2323")
      alter_fun = list(
        background = alter_graphic("rect", fill = "#CCCCCC"),
        SNV = alter_graphic("rect", fill = col["SNV"]),
        insertion = alter_graphic("rect", height = 0.33, fill = col["insertion"]),
        deletion = alter_graphic("rect", height = 0.33, fill = col["deletion"],),
        sequence_alteration = alter_graphic("rect", height = 0.33, fill = col["sequence_alteration"],),
        indel = alter_graphic("rect", height = 0.33, fill = col["indel"],)
      )
      
      column_title <- "Gastrointestinal cancer plasma sample, genes in COSMIC"
      heatmap_legend_param <-
        list(
          title = "Alternations",
          at = c("SNV", "insertion", "deletion","sequence_alteration"),
          labels = c("SNV", "Insertion", "Deletion","Sequence alteration")
        )
      
      mutation_example <- oncoPrint(
        mat, alter_fun = alter_fun, col = col,
        remove_empty_rows = TRUE,
        column_order = colnames(mat),
        #top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot()),
        top_annotation = NULL,
        right_annotation = NULL,
        left_annotation = NULL,
        bottom_annotation = NULL,
        #column_order = colnames(mat),
        #row_order = unlist(lapply(strsplit(order,"|",fixed = TRUE),function(x) x[2])),
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_title = column_title,
        show_heatmap_legend = FALSE,
        heatmap_legend_param = heatmap_legend_param,
        #bottom_annotation = colAnn
        #top_annotation = HeatmapAnnotation(
        #  cbar = anno_oncoprint_barplot(),
        #  show_annotation_name = TRUE,
        #  Age = 1:ncol(mat),
        #  Gender = 1:ncol(mat),
        #  Stage = 1:ncol(mat)
        #)
      )
      mutation_example <- as.ggplot(mutation_example)
    }
    
    ##RNA chimeric gene examples
    {
      ann_plot <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/13.exact_test_for_SNP/oncoplot/RNA_sample_annotation.csv", header = TRUE,)
      mat_plot <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/JunctionReadCount.txt",header = TRUE, row.names = 1,sep = "\t",check.names = FALSE)
      #mat_plot[] <- lapply(mat_plot, as.character) 
      mat_plot[is.na(mat_plot)] <- ""
      
      ann <- ann_plot[grep("STAD|CRC",ann_plot$Samples),]
      mat <- mat_plot[grep("TCGA",rownames(mat_plot)),]
      mat <- mat[,colnames(mat) %in% ann$Samples]
      
      ann <- ann[,-which(colnames(ann)=="Samples")]
      colours_ann <- list(
        'Gender'=c('Male'='#5CACEE','Female'='#FF0000'),
        'Stage'=c('Stage I'='#DCDCDC','Stage II'='#A3A3A3','Stage III'='#4A4A4A','Stage IV'='#000000'),
        'Location'=c('No biopsy'='#FFFFFF','Others'='#DCDCDC',
                     'Rectum'='#FFF8DC','Sigmoid colon'='#EEDD82','Descending colon'='#FEE5AC','Transverse colon'='#EDCB62','Ascending colon'='#FCB514',
                     'Fundus'='#FEF1B5','Body'='#E5BC3B','Antrum'='#CD9B1D'),
        'MMR'=c('No biopsy'='#FFFFFF','p'='#EDCB62','d'='#FCB514'),
        'HER2'=c('No biopsy'='#FFFFFF','-'='#E5BC3B','+'='#CD9B1D'),
        'Type'=c('Colorectal cancer'='#FCB514','Stomach cancer'='#CD9B1D')
      )
      
      colAnn <- HeatmapAnnotation(col = colours_ann,
                                  #cbar = anno_oncoprint_barplot(),
                                  df = ann,which = 'col')
      
      
      mat$Chimeric <- as.character(lapply(strsplit(rownames(mat),"|",fixed = TRUE),function(x) x[1]))
      mat <- aggregate(mat[,-which(colnames(mat)=="Chimeric")], by = list("Chimeric"=mat$Chimeric),FUN = sum)
      rownames(mat) <- mat$Chimeric
      mat <- mat[,-which(colnames(mat)=="Chimeric")]
      mat[mat>0] <- "TCGA"
      mat[mat==0] <- ""
      #mat <- mat_plot[order,]
      
      #deep sequencing
      #mat <- mat[,grep("CRC.PKU.32.pico|NC.PKU.10.pico|NC.PKU.14.pico|CRC.PKU.29.pico|NC.PKU.12.pico|CRC.PKU.30.pico|CRC.PKU.27.pico|CRC.PKU.28.pico|NC.PKU.9.pico|CRC.PKU.31.pico",colnames(mat))]
      #paired sample
      
      
      
      
      col <- c("TCGA" = "#ff7f00")
      alter_fun = list(
        background = alter_graphic("rect", fill = "#CCCCCC"),
        TCGA = alter_graphic("rect", fill = col["TCGA"])
      )
      
      column_title <- "Gastrointestinal cancer plasma sample, genes in TCGA"
      heatmap_legend_param <-
        list(
          title = "Alternations",
          at = c("TCGA"),
          labels = c("TCGA")
        )
      
      chimeric_example <- oncoPrint(
        mat, alter_fun = alter_fun, col = col,
        column_order = colnames(mat),
        remove_empty_rows = TRUE,
        top_annotation = NULL,
        right_annotation = NULL,
        left_annotation = NULL,
        bottom_annotation = NULL,
        #column_order = colnames(mat),
        #row_order = unlist(lapply(strsplit(order,"|",fixed = TRUE),function(x) x[2])),
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_title = NULL,
        show_heatmap_legend = FALSE,
        heatmap_legend_param = heatmap_legend_param,
        #bottom_annotation = colAnn
        #top_annotation = HeatmapAnnotation(
        #  cbar = anno_oncoprint_barplot(),
        #  show_annotation_name = TRUE,
        #  Age = 1:ncol(mat),
        #  Gender = 1:ncol(mat),
        #  Stage = 1:ncol(mat)
        #)
      )
      chimeric_example <- as.ggplot(chimeric_example)
    }
    
    ann_plot <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/sample/ML_samples.csv", header = TRUE,check.names = FALSE)
    ann_plot <- ann_plot[!(ann_plot$sample_id %in% c("NC-PKU-mix3","NC-PKU-mix4","NC-PKU-mix6","NC-PKU-mix7","NC-PKU-mix8","NC-PKU-mix9")),]
    ann_plot$sample_id <- factor(ann_plot$sample_id,levels = as.character(ann_plot$sample_id))
    
    #Summary spliced genes
    {
      Splicing <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/multiomics_paired_Cancer_vs_Healthy_Inclevel_matrix.txt",sep = "\t", header = TRUE,check.names=FALSE)
      Splicing_summary <- data.frame(matrix(nrow=5,ncol = length(colnames(Splicing))))
      rownames(Splicing_summary) <- c("A3SS","A5SS","MXE","RI","SE")
      colnames(Splicing_summary) <- colnames(Splicing)
      for(sample in colnames(Splicing)){
        #print(sample)
        single_Splicing <- na.omit(Splicing[,sample,drop=FALSE])
        single_Splicing$Categroy <- as.character(lapply(strsplit(rownames(single_Splicing),"|",fixed = TRUE),function(x) x[1]))
        A3SS <- single_Splicing %>% summarise(Number=sum(.=="A3SS"))
        A5SS <- single_Splicing %>% summarise(Number=sum(.=="A5SS"))
        MXE <- single_Splicing %>% summarise(Number=sum(.=="MXE"))
        RI <- single_Splicing %>% summarise(Number=sum(.=="RI"))
        SE <- single_Splicing %>% summarise(Number=sum(.=="SE"))
        
        Splicing_summary["A3SS",sample] <- A3SS
        Splicing_summary["A5SS",sample] <- A5SS
        Splicing_summary["MXE",sample] <- MXE
        Splicing_summary["RI",sample] <- RI
        Splicing_summary["SE",sample] <- SE
      }
    dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Splicing/")
    write.csv(Splicing_summary,"/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Splicing/Splicing_summary.csv")
  
    Splicing_summary <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Splicing/Splicing_summary.csv",check.names = FALSE,header = TRUE,row.names = 1,stringsAsFactors = FALSE)
    Splicing_summary_plot <- reshape2::melt(t(Splicing_summary))
    colnames(Splicing_summary_plot) <- c("Sample","Splicing_type","Number")
    #Splicing_summary_plot <- Splicing_summary_plot[-grep("NC",Splicing_summary_plot$Sample),]
    #Splicing_summary_plot <- Splicing_summary_plot[-grep("mix..pico",Splicing_summary_plot$Sample),]
    Splicing_summary_plot <- Splicing_summary_plot[as.character(Splicing_summary_plot$Sample) %in% as.character(ann_plot$RNA_id),]
    
    Splicing_summary_plot$Sample <- gsub("-pico","",Splicing_summary_plot$Sample)
    Splicing_summary_plot$Sample <- factor(Splicing_summary_plot$Sample,levels = as.character(ann_plot$sample_id))
    p_splicing <- ggplot(Splicing_summary_plot,aes(x=Sample,y=Number,fill=Splicing_type))+
      geom_bar(stat = "identity")+
      scale_y_continuous(expand = c(0,0))+
      theme_bw()+
      theme(legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(size=0.5, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=6),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=6),
        legend.key.size = unit(x=c(1.5),units="pt"),
        plot.title = element_text(hjust = 0.5,size=12,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12))
    }
    
    #Summary chimeric genes
    {
      Chimeric <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/JunctionReadCount.txt",sep = "\t", header = TRUE,check.names=FALSE,row.names = 1,stringsAsFactors = FALSE)
      Chimeric_summary <- data.frame(matrix(nrow=2,ncol = length(colnames(Chimeric))))
      rownames(Chimeric_summary) <- c("Interchromosomal","Intrachromosomal")
      colnames(Chimeric_summary) <- colnames(Chimeric)
      for(sample in colnames(Chimeric)){
        #print(sample)
        single_Chimeric <- Chimeric[Chimeric[[sample]] > 0,sample,drop=FALSE]
        if(nrow(single_Chimeric)==0){
          Chimeric_summary["Interchromosomal",sample] <- 0
          Chimeric_summary["Intrachromosomal",sample] <- 0
          next
        }
        single_Chimeric$Chimeric_Type <- as.character(lapply(strsplit(rownames(single_Chimeric),"|",fixed = TRUE),function(x) x[length(x)]))
        single_Chimeric$Chimeric_Type <- as.character(lapply(strsplit(single_Chimeric$Chimeric_Type,"INT",fixed = TRUE),function(x) x[2]))
        single_Chimeric$Chimeric_Type <- as.character(lapply(strsplit(single_Chimeric$Chimeric_Type,"CHROMOSOMAL",fixed = TRUE),function(x) x[1]))
        single_Chimeric$Chimeric_Type <- gsub("ER","Interchromosomal",single_Chimeric$Chimeric_Type)
        single_Chimeric$Chimeric_Type <- gsub("RA","Intrachromosomal",single_Chimeric$Chimeric_Type)
      
        Interchromosomal <- single_Chimeric %>% summarise(Number=sum(.=="Interchromosomal"))
        Intrachromosomal <- single_Chimeric %>% summarise(Number=sum(.=="Intrachromosomal"))
        
        Chimeric_summary["Interchromosomal",sample] <- Interchromosomal
        Chimeric_summary["Intrachromosomal",sample] <- Intrachromosomal
      }
      dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Chimeric/")
      write.csv(Chimeric_summary,"/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Chimeric/Chimeric_summary.csv")
      
      Chimeric_summary <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Chimeric/Chimeric_summary.csv",check.names = FALSE,header = TRUE,row.names = 1,stringsAsFactors = FALSE)
      Chimeric_summary_plot <- reshape2::melt(t(Chimeric_summary))
      colnames(Chimeric_summary_plot) <- c("Sample","Chimeric_type","Number")
      #Chimeric_summary_plot <- Chimeric_summary_plot[-grep("NC",Chimeric_summary_plot$Sample),]
      #Chimeric_summary_plot <- Chimeric_summary_plot[-grep("mix..pico",Chimeric_summary_plot$Sample),]
      Chimeric_summary_plot <- Chimeric_summary_plot[as.character(Chimeric_summary_plot$Sample) %in% as.character(ann_plot$RNA_id),]
      Chimeric_summary_plot$Sample <- gsub("-pico","",Chimeric_summary_plot$Sample)
      Chimeric_summary_plot$Sample <- factor(Chimeric_summary_plot$Sample,levels = as.character(ann_plot$sample_id))
      p_chimeric <- ggplot(Chimeric_summary_plot,aes(x=Sample,y=Number,fill=Chimeric_type))+
        geom_bar(stat = "identity")+
        scale_y_continuous(expand = c(0,0))+
        theme_bw()+
        theme(legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_line(size=0.5, colour = "black"),
              legend.title = element_text(face="bold", color="black",family = "Arial", size=6),
              legend.text= element_text(face="bold", color="black",family = "Arial", size=6),
              legend.key.size = unit(x=c(1.5),units="pt"),
              plot.title = element_text(hjust = 0.5,size=12,face="bold"),
              axis.text.x = element_blank(),
              #axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
              axis.text.y = element_text(color="black", size=12),
              axis.title.x = element_text(color="black", size=12),
              axis.title.y = element_text(color="black", size=12))
    }
    
    #Summary Editing sites
    {
      Editing <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/Editing_ratio_multiomics_paired_>10read.txt",sep = "\t", header = TRUE,check.names=FALSE,row.names = 1)
      chr <- as.character(lapply(strsplit(rownames(Editing),"_"),function(x) x[1]))
      chr <- gsub("chr","",chr)
      start <- as.character(lapply(strsplit(rownames(Editing),"_"),function(x) x[2]))
      end <- as.character(lapply(strsplit(rownames(Editing),"_"),function(x) x[2]))
      #unique(as.character(lapply(strsplit(rownames(Editing),"_"),function(x) x[4])))
      #unique(as.character(lapply(strsplit(rownames(Editing),"_"),function(x) x[5])))
      
      i=1
      final_biotype <- {}
      while(i<=nrow(Editing)){
      message(paste0(i,"/",nrow(Editing)))
      values <- list(chromosome=chr[i],start=start[i],end=start[i])
      biotype <- getBM(attributes=c("ensembl_gene_id","gene_biotype"),
            filters=c("chromosome_name","start","end"), 
            values=values, mart=mart, useCache = FALSE)
      if(nrow(biotype)==0){
        biotype[1,"ensembl_gene_id"] <- NA
        biotype[1,"gene_biotype"] <- NA
        #biotype[1,"transcript_biotype"] <- NA
        biotype$ID <- rownames(Editing)[i]
        final_biotype <- rbind(final_biotype,biotype)
        i=i+1
      } else {
      biotype$ID <- rownames(Editing)[i]
      final_biotype <- rbind(final_biotype,biotype[1,])
      i=i+1
      }
      
      }
      
      
      final_biotype[is.na(final_biotype$gene_biotype),] <- "Unannotated"
      dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Editing/")
      write.csv(final_biotype,"/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Editing/Editing_summary.csv")
    
      Editing_summary <- data.frame(matrix(nrow=length(unique(final_biotype$gene_biotype)),ncol = length(colnames(Editing))))
      rownames(Editing_summary) <- c(unique(final_biotype$gene_biotype))
      colnames(Editing_summary) <- colnames(Editing)
      for(sample in colnames(Editing)){
        #print(sample)
        single_Editing <- Editing[Editing[[sample]] > 0,sample,drop=FALSE]
        single_Editing$ID <- rownames(single_Editing)
        single_Editing <- left_join(single_Editing,final_biotype,by=c('ID'='ID'))
        single_Editing[is.na(single_Editing)] <- "Unannotated"
        
        single_Editing[,-grep("ensembl_gene_id",colnames(single_Editing))]
        for(type in unique(final_biotype$gene_biotype)){
        Editing_summary[type,sample] <- single_Editing %>% summarise(Number=sum(gene_biotype==type))
        }
      }
      dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Editing/")
      write.csv(Editing_summary,"/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Editing/Editing_summary.csv")
      
      Editing_summary <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Editing/Editing_summary.csv",check.names = FALSE,header = TRUE,row.names = 1,stringsAsFactors = FALSE)
      Editing_summary_plot <- reshape2::melt(t(Editing_summary))
      colnames(Editing_summary_plot) <- c("Sample","Editing_type","Number")
      #Editing_summary_plot <- Editing_summary_plot[-grep("NC",Editing_summary_plot$Sample),]
      #Editing_summary_plot <- Editing_summary_plot[-grep("mix..pico",Editing_summary_plot$Sample),]
      Editing_summary_plot <- Editing_summary_plot[as.character(Editing_summary_plot$Sample) %in% as.character(ann_plot$RNA_id),]
      Editing_summary_plot$Sample <- gsub("-pico","",Editing_summary_plot$Sample)
      Editing_summary_plot$Sample <- factor(Editing_summary_plot$Sample,levels = as.character(ann_plot$sample_id))
      p_Editing <- ggplot(Editing_summary_plot,aes(x=Sample,y=Number,fill=Editing_type))+
        geom_bar(stat = "identity")+
        scale_y_continuous(expand = c(0,0))+
        theme_bw()+
        theme(legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_line(size=0.5, colour = "black"),
              legend.title = element_text(face="bold", color="black",family = "Arial", size=6),
              legend.text= element_text(face="bold", color="black",family = "Arial", size=6),
              legend.key.size = unit(x=c(1.5),units="pt"),
              plot.title = element_text(hjust = 0.5,size=24,face="bold"),
              axis.text.x = element_blank(),
              #axis.text.x = element_text(color="black", size=6,angle = 45,hjust = 1),
              axis.text.y = element_text(color="black", size=12),
              axis.title.x = element_text(color="black", size=12),
              axis.title.y = element_text(color="black", size=12))
    }
    
    #summary ASE
    {
      ASE <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/ASE_COSMIC_multiomics_paired_>500read.txt",sep = "\t", header = TRUE,check.names=FALSE,row.names = 1)
      
      ASE_summary <- {}
      for(sample in colnames(ASE)){
        #print(sample)
        single_ASE <- ASE[ASE[[sample]]<0.49,sample,drop=FALSE]
        single_ASE$ID <- sample
        colnames(single_ASE) <- c("AE","Sample")
        ASE_summary <- rbind(ASE_summary,single_ASE)
      }
      
      #ASE_summary_plot <- ASE_summary[-grep("NC",ASE_summary$sample),]
      #ASE_summary_plot <- ASE_summary_plot[-grep("mix..pico",ASE_summary$sample),]
      ASE_summary_plot <- ASE_summary[as.character(ASE_summary$Sample) %in% as.character(ann_plot$RNA_id),]
      ASE_summary_plot$Sample <- gsub("-pico","",ASE_summary_plot$Sample)
      ASE_summary_plot$Sample <- factor(ASE_summary_plot$Sample,levels = as.character(ann_plot$sample_id))
      p_ASE <- ggplot(ASE_summary_plot,aes(x=Sample,y=AE,fill="AE"))+
        geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.colour = "red",outlier.size=0,outlier.alpha = 0)+
        #geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 0.8),color = alpha("black",alpha = 0.2))+
        scale_y_continuous(limits = c(0,0.5),expand = c(0,0))+
        theme_bw()+
        theme(legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_line(size=0.5, colour = "black"),
              legend.title = element_text(face="bold", color="black",family = "Arial", size=12),
              legend.text= element_text(face="bold", color="black",family = "Arial", size=12),
              plot.title = element_text(hjust = 0.5,size=24,face="bold"),
              axis.text.x = element_blank(),
              #axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
              axis.text.y = element_text(color="black", size=12),
              axis.title.x = element_text(color="black", size=12),
              axis.title.y = element_text(color="black", size=12))
      p_ASE
    }
    
    #summary Alternative promoter
    {
      #too slow, run on cluster for each sample, then rbind
      {
      Altpromoter_normalized <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/AlternativePromoter_multiomics_paired_normalized.txt",sep = "\t", header = TRUE,check.names=FALSE,row.names = 1)
      Altpromoter_activity <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/TPM-by-promoter_multiomics_paired.txt",sep = "\t", header = TRUE,check.names=FALSE,row.names = 1)
      
      genes <- unique(as.character(lapply(strsplit(rownames(Altpromoter_normalized),".",fixed = TRUE), function(x) x[1])))
      
      pb <- progress_bar$new(
        format = "  Processing [:bar] :percent eta: :eta",
        total = length(colnames(Altpromoter_normalized)), clear = FALSE, width= 60)
      
      Promoter_usage_summary <- {}
      for(sample in colnames(Altpromoter_normalized)){
        single_sample_promoter_usage <- data.frame(matrix(nrow=3,ncol = length(genes)))
        colnames(single_sample_promoter_usage) <- genes
        rownames(single_sample_promoter_usage) <- paste0(c("inactivate","major","minor"),"_",sample)
        for(gene in genes[1:100]){
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
        Promoter_usage_summary <- rbind(Promoter_usage_summary,single_sample_promoter_usage)
        pb$tick()
        Sys.sleep(1 / 100)
      }
      }
      
      #summary major, minor, and inactivate promoter
      {
        Promoter_usage <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/Promoter_usage.txt",header = TRUE,sep = "\t",row.names = 1,stringsAsFactors = FALSE)
        
        Promoter_usage_summary <- rowMeans(Promoter_usage[,-grep("sample|promoter",colnames(Promoter_usage))])
        
        Promoter_usage_summary <- data.frame(sample=as.character(lapply(strsplit(names(Promoter_usage_summary),"_"),function(x) x[2])),
                                                  average_promoter_usage=Promoter_usage_summary,
                                                  group = as.character(lapply(strsplit(names(Promoter_usage_summary),"_"),function(x) x[1])))
        
        Promoter_usage_summary_plot <- Promoter_usage_summary[as.character(Promoter_usage_summary$sample) %in% as.character(ann_plot$RNA_id),]
        Promoter_usage_summary_plot$sample <- gsub("-pico","",Promoter_usage_summary_plot$sample)
        Promoter_usage_summary_plot$sample <- factor(Promoter_usage_summary_plot$sample,levels = as.character(ann_plot$sample_id))
        p_Promoter <- ggplot(Promoter_usage_summary_plot,aes(x=sample,y=average_promoter_usage,fill=group))+
          geom_bar(stat = "identity")+
          scale_y_continuous(expand = c(0,0))+
          theme_bw()+
          theme(legend.position="right",
                panel.grid=element_blank(),
                panel.border=element_blank(),
                axis.ticks.x = element_blank(),
                axis.line.x = element_blank(),
                axis.line.y = element_line(size=0.5, colour = "black"),
                legend.title = element_text(face="bold", color="black",family = "Arial", size=6),
                legend.text= element_text(face="bold", color="black",family = "Arial", size=6),
                legend.key.size = unit(x=c(1.5),units="pt"),
                plot.title = element_text(hjust = 0.5,size=24,face="bold"),
                axis.text.x = element_blank(),
                #axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
                axis.text.y = element_text(color="black", size=12),
                axis.title.x = element_text(color="black", size=12),
                axis.title.y = element_text(color="black", size=12))
        
        for(sample in unique(Promoter_usage$sample)){
          
          major <- t(Promoter_usage[Promoter_usage$sample==sample & Promoter_usage$promoter == "major",
                                    -grep("sample|promoter",colnames(Promoter_usage))])[,1]
          length(major[major > 0.5])  
          minor <- t(Promoter_usage[Promoter_usage$sample==sample & Promoter_usage$promoter == "minor",
                                    -grep("sample|promoter",colnames(Promoter_usage))])[,1]
          length(minor[minor > 0.5])  
          inactivate <- t(Promoter_usage[Promoter_usage$sample==sample & Promoter_usage$promoter == "inactivate",
                                    -grep("sample|promoter",colnames(Promoter_usage))])[,1]
          length(inactivate[inactivate > 0.5])  
        }
        
        
      }
      
    }
    
    #summary methylation
    {
      Methylation_summary <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Medip/cfMeDIP_count_matrix_allregions.txt",header = TRUE,sep = "\t",stringsAsFactors = FALSE)
      
      #Methylation_summary_plot <- Methylation_summary[-grep("rep",Methylation_summary$sample),]
      Methylation_summary_plot <- Methylation_summary[as.character(Methylation_summary$sample) %in% as.character(ann_plot$Methylation_id),]
      Methylation_summary_plot$sample <- gsub("-me","",Methylation_summary_plot$sample)
      Methylation_summary_plot$sample <- factor(Methylation_summary_plot$sample,levels = as.character(ann_plot$sample_id))
      p_Methylation <- ggplot(Methylation_summary_plot,aes(x=sample,y=assign.ratio,fill=region))+
        geom_bar(stat = "identity")+
        scale_y_continuous(expand = c(0,0))+
        theme_bw()+
        theme(legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_line(size=0.5, colour = "black"),
              legend.title = element_text(face="bold", color="black",family = "Arial", size=6),
              legend.text= element_text(face="bold", color="black",family = "Arial", size=6),
              legend.key.size = unit(x=c(1.5),units="pt"),
              plot.title = element_text(hjust = 0.5,size=24,face="bold"),
              axis.text.x = element_blank(),
              #axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
              axis.text.y = element_text(color="black", size=12),
              axis.title.x = element_text(color="black", size=12),
              axis.title.y = element_text(color="black", size=12))
    }
    
    #summary copy number
    {
      CNV_summary <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/CNV/sum-TF-10k-largeNormal.txt",header = TRUE,sep = "\t",stringsAsFactors = FALSE,row.names = NULL)
      CNV_summary <- CNV_summary[,c(1,2)]
      colnames(CNV_summary) <- c("sample","inchorCNA")

      CNV_summary_plot <- CNV_summary[as.character(CNV_summary$sample) %in% as.character(ann_plot$DNA_id),]
      CNV_summary_plot$sample <- gsub("-wgs","",CNV_summary_plot$sample)
      CNV_summary_plot$sample <- factor(CNV_summary_plot$sample,levels = as.character(ann_plot$sample_id))
      p_CNV <- ggplot(CNV_summary_plot,aes(x=sample,y=inchorCNA))+
        geom_bar(stat = "identity")+
        scale_y_continuous(expand = c(0,0))+
        theme_bw()+
        theme(legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_line(size=0.5, colour = "black"),
              legend.title = element_text(face="bold", color="black",family = "Arial", size=6),
              legend.text= element_text(face="bold", color="black",family = "Arial", size=6),
              legend.key.size = unit(x=c(1.5),units="pt"),
              plot.title = element_text(hjust = 0.5,size=24,face="bold"),
              axis.text.x = element_blank(),
              #axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
              axis.text.y = element_text(color="black", size=12),
              axis.title.x = element_text(color="black", size=12),
              axis.title.y = element_text(color="black", size=12))
    }
    
    #summary fragment size
    {
      FragSize_summary <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/DNA-FragRatio_matrix_bin100kb.correctGC.txt",row.names = 1,check.names = FALSE,header = TRUE,sep = "\t",stringsAsFactors = FALSE)
      
      FragSize_summary <- data.frame(colMeans(FragSize_summary))
      FragSize_summary$sample <- rownames(FragSize_summary)
      
      FragSize_summary_plot <- FragSize_summary[-grep("KAPA|NV|Y",FragSize_summary$sample),]
      FragSize_summary_plot <- FragSize_summary_plot[as.character(FragSize_summary_plot$sample) %in% as.character(ann_plot$DNA_id),]
      FragSize_summary_plot$sample <- gsub("-wgs","",FragSize_summary_plot$sample)
      FragSize_summary_plot$sample <- factor(FragSize_summary_plot$sample,levels = as.character(ann_plot$sample_id))
      p_FragSize <- ggplot(FragSize_summary_plot,aes(x=sample,y=colMeans.FragSize_summary.))+
        geom_point(stat = "identity")+
        scale_y_continuous(expand = c(0.001,0.001))+
        theme_bw()+
        theme(legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_line(size=0.5, colour = "black"),
              legend.title = element_text(face="bold", color="black",family = "Arial", size=6),
              legend.text= element_text(face="bold", color="black",family = "Arial", size=6),
              legend.key.size = unit(x=c(1.5),units="pt"),
              plot.title = element_text(hjust = 0.5,size=24,face="bold"),
              axis.text.x = element_blank(),
              #axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
              axis.text.y = element_text(color="black", size=12),
              axis.title.x = element_text(color="black", size=12),
              axis.title.y = element_text(color="black", size=12))
    }
    
    #RNA chimeric goem_tile
    {
      mat_plot <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/JunctionReadCount.txt",header = TRUE,row.names = 1,sep = "\t",check.names = FALSE)
      #mat_plot[] <- lapply(mat_plot, as.character) 
      mat_plot[is.na(mat_plot)] <- ""
      
      mat <- mat_plot[grep("TCGA",rownames(mat_plot)),]
      mat <- mat[,as.character(colnames(mat)) %in% as.character(ann_plot$RNA_id)]
      mat <- mat[which(rowSums(mat)!=0),]
      mat <- mat[which(rowSums(mat[,grep("NC",colnames(mat))])==0),]
      
      mat$Chimeric <- as.character(lapply(strsplit(rownames(mat),"|",fixed = TRUE),function(x) x[1]))
      mat <- aggregate(mat[,-which(colnames(mat)=="Chimeric")], by = list("Chimeric"=mat$Chimeric),FUN = sum)
      rownames(mat) <- mat$Chimeric
      #mat <- mat[grep("USP15--DPY19L2|TNFAIP8--DMXL1|SOD2--GNAS|FARS2--CDYL|RNF138--RNF125|PFKFB3--GDI2|NCK1--STAG1|IL13RA1--DOCK11|CNOT4--CALD1|ITFG1--PHKB|TNFAIP8--DMXL1|TRIP12--AGFG1|SETD2--SMARCC1|PPP2R2D--RGS10|USP4--RHOA|RAP1B--MDM1|TULP4--GTF2H5|OAZ1--SF3A2|SMAD2--CD226|SGMS1--PARG|ACER3--AP001189.4|AKT3--SDCCAG8|FER--FBXL17|R3HDM2--RAP1B|BACH1--MAP3K7CL",rownames(mat)),]
      #mat <- mat[,-which(colnames(mat)=="Chimeric")]
      #mat[mat>0] <- "TCGA"
      #mat[mat==0] <- ""
      
      dft <-gather(mat, sample, value, 2:ncol(mat))
      
      dft[dft$value>0,"value"] <- "Chimeric"
      dft[dft$value==0,"value"] <- "Not detected"
      dft$value <- factor(dft$value,levels=c("Not detected","Chimeric"))
      dft$sample <- gsub("-pico","",dft$sample)
      dft$sample <- factor(dft$sample,levels = as.character(ann_plot$sample_id))
      
      Chimeric_tile <- ggplot(dft,aes(x=sample,y=Chimeric))+
        geom_tile(aes(fill=value))+
        scale_fill_manual(values = c('white','black'))+
        theme_bw()+
        theme(legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_line(size=0.5, colour = "black"),
              legend.title = element_text(face="bold", color="black",family = "Arial", size=6),
              legend.text= element_text(face="bold", color="black",family = "Arial", size=6),
              legend.key.size = unit(x=c(1.5),units="pt"),
              plot.title = element_text(hjust = 0.5,size=24,face="bold"),
              axis.text.x = element_blank(),
              #axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
              axis.text.y = element_text(color="black", size=12,hjust = 1),
              axis.title.x = element_text(color="black", size=12),
              axis.title.y = element_text(color="black", size=12))
      
    }
    
    #RNA mutation goem_tile
    {
      mat_plot <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/SNP/Mutation_burden_class_forplot.txt",header = TRUE,row.names = 1,sep = "\t",check.names = FALSE)
      #mat_plot[] <- lapply(mat_plot, as.character) 
      mat_plot[is.na(mat_plot)] <- ""
      
      mat <- mat_plot[grep("ENSG00000142208|ENSG00000134982|ENSG00000175054|ENSG00000103126|ENSG00000168646|ENSG00000166710|ENSG00000186174|ENSG00000157764|ENSG00000039068|ENSG00000168036|ENSG00000257923|ENSG00000104408|ENSG00000100393|ENSG00000141736|ENSG00000065361|ENSG00000178568|ENSG00000109670|ENSG00000066468|ENSG00000183454|ENSG00000100644|ENSG00000133703|ENSG00000121879|ENSG00000145675|ENSG00000101213|ENSG00000163629|ENSG00000196090|ENSG00000112531|ENSG00000164754|ENSG00000067560|ENSG00000175387|ENSG00000166949|ENSG00000141646|ENSG00000177565|ENSG00000148737|ENSG00000163513|ENSG00000141510|ENSG00000104517|ENSG00000140836",rownames(mat_plot)),]
      mat <- mat_plot[grep("ENSG00000133703|ENSG00000141510|ENSG00000134982|ENSG00000076242|ENSG00000095002|ENSG00000116062|ENSG00000122512|ENSG00000141736|ENSG00000146648",rownames(mat_plot)),]
      mat <- mat[,colnames(mat) %in% as.character(ann_plot$RNA_id)]
      mat$Mutation <- as.character(lapply(strsplit(rownames(mat),"|",fixed = TRUE),function(x) x[2]))
      
      dft <-gather(mat, sample, value, 1:(ncol(mat)-1))
      
      #dft[dft$value>0,"value"] <- "Chimeric"
      #dft[dft$value==0,"value"] <- "Not detected"
      dft[is.na(dft$value),"value"] <- "Not detected"
      #dft$value <- gsub(";","",dft$value)
      dft$value <- factor(dft$value,levels=c("Not detected","SNV;","insertion;","deletion;","sequence_alteration;",
                                             "indel;SNV;","deletion;SNV;","SNV;insertion;","SNV;deletion;","SNV;sequence_alteration;",
                                             "SNV;deletion;SNV;","SNV;insertion;SNV;","insertion;SNV;"))
      dft$sample <- gsub("-pico","",dft$sample)
      dft$sample <- factor(dft$sample,levels = as.character(ann_plot$sample_id))
      
      Mutation_tile <- ggplot(dft,aes(x=sample,y=Mutation))+
        geom_tile(aes(fill=value))+
        scale_fill_manual(values = c('white','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black'))+
        theme_bw()+
        theme(legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_line(size=0.5, colour = "black"),
              legend.title = element_text(face="bold", color="black",family = "Arial", size=6),
              legend.text= element_text(face="bold", color="black",family = "Arial", size=6),
              legend.key.size = unit(x=c(1.5),units="pt"),
              plot.title = element_text(hjust = 0.5,size=24,face="bold"),
              axis.text.x = element_blank(),
              #axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
              axis.text.y = element_text(color="black", size=12,hjust = 1),
              axis.title.x = element_text(color="black", size=12),
              axis.title.y = element_text(color="black", size=12))
      
    }
    
    #clinical
    {
      clinical_Cancer <-
        ggplot(ann_plot,aes(x=sample_id,y=0,fill= as.character(ann_plot$Group)))+
        geom_tile()+
        #scale_fill_gradient(low="white",high = "black")+
        scale_fill_manual(values=c("CRC"=alpha("#FCB514",alpha = 1),"STAD"=alpha("red",alpha = 1),"HD" = alpha("blue",alpha = 1)))+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(0,-5,0,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          #legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x= element_blank(),
          #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      clinical_Gender <-
        ggplot(ann_plot,aes(x=sample_id,y=0,fill= as.character(ann_plot$Gender)))+
        geom_tile()+
        #scale_fill_gradient(low="white",high = "black")+
        scale_fill_manual(values=c("F"=alpha("#FF0000",alpha = 1),"M"=alpha("#5CACEE",alpha = 1)))+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(0,-5,0,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          #legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x= element_blank(),
          #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      clinical_Age <-
        ggplot(ann_plot,aes(x=sample_id,y=0,fill= as.numeric(as.character(ann_plot$Age))))+#`Tumor size`))+
        geom_tile()+
        scale_fill_gradient(low="white",high = "orange")+
        #scale_fill_manual(values=c("-"="white","early"="sky blue","late"="red"))+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(0,-5,0,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          #legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x= element_blank(),
          #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      clinical_P53 <-
        ggplot(ann_plot,aes(x=sample_id,y=0,fill= P53))+#`Tumor size`))+
        geom_tile()+
        #scale_fill_gradient(low="white",high = "black")+
        scale_fill_manual(values=c("NA"="grey","-"="white","+"="white","++"="white","+++"="black"))+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(0,-5,0,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          #legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x= element_blank(),
          #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      clinical_MMR <-
        ggplot(ann_plot,aes(x=sample_id,y=0,fill= `MSI/MMR`))+#`Tumor size`))+
        geom_tile()+
        #scale_fill_gradient(low="white",high = "black")+
        scale_fill_manual(values=c("NA"="grey","x"="grey","Normal"="white","Unstable"="#F05C3BFF","Unstable; lack"="#F05C3BFF","unstable; lack"="#F05C3BFF"))+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(0,-5,0,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          #legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x= element_blank(),
          #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      clinical_HER2 <-
        ggplot(ann_plot,aes(x=sample_id,y=0,fill= HER2))+#`Tumor size`))+
        geom_tile()+
        #scale_fill_gradient(low="white",high = "black")+
        scale_fill_manual(values=c("NA"="grey","-"="white","+"="white","++"="red","+++"="red"))+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(0,-5,0,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          #legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x= element_blank(),
          #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      clinical_Right <-
        ggplot(ann_plot,aes(x=sample_id,y=0,fill= ann_plot$Right_left_half))+#`Tumor size`))+
        geom_tile()+
        #scale_fill_gradient(low="white",high = "black")+
        scale_fill_manual(values=c("NA"="grey","x"="grey","L"="white","R"="dark red","L+R"="dark red"))+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(0,-5,0,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          #legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x= element_blank(),
          #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      clinical_Stage <-
        ggplot(ann_plot,aes(x=sample_id,y=0,fill= as.character(ann_plot$Stage_detailed)))+#`Tumor size`))+
        geom_tile()+
        #scale_fill_gradient(low="white",high = "black")+
        scale_fill_manual(values=c("NA"="grey","x"="grey","1"="#90E0EF","1A"="#90E0EF","1B"="#90E0EF","2A"="#00B4D8","2B"="#00B4D8","2C"="#00B4D8","3A"="#03045E","3B"="#03045E","3C"="#03045E","4"="black"))+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(0,-5,0,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          #legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x= element_blank(),
          #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
    }
    
    ggarrange(p_Promoter,p_splicing,p_chimeric,p_Editing,p_CNV,p_Methylation,mutation_example,ncol=1,align = "v")
    
    overview <- plot_grid(
      p_CNV + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                    axis.text.y = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
      p_Methylation + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                            axis.text.y = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
      p_FragSize + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                         axis.text.y = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        p_Promoter + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                           axis.text.y = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(2,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"), 
        p_splicing + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                           axis.text.y = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        p_chimeric + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                           axis.text.y = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        p_Editing + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                          axis.text.y = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        p_ASE + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                      axis.text.y = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        Chimeric_tile + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                              axis.text.y = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,2,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        Mutation_tile + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                              axis.text.y = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        clinical_P53 + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                             axis.text = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        clinical_MMR + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                             axis.text = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        clinical_HER2 + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                             axis.text = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        clinical_Gender + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                             axis.text = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        clinical_Age + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                             axis.text = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        clinical_Cancer + theme(plot.title = element_text(hjust = 1,size=6,face="bold"),axis.title.y = element_blank(),
                             axis.text = element_text(hjust = 1,size=6,face="bold"),plot.margin = unit(x=c(0,10,1,10),units="pt"),axis.title.x = element_blank(),legend.position = "none"),
        ncol = 1, align = "hv",
        rel_widths = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
        rel_heights = c(1,1,1,1,1,1,1,1,6,2,0.5,0.5,0.5,0.5,0.5,0.5)
    )
    
    legend <- plot_grid(
      get_legend(p_CNV),
      get_legend(p_Methylation),
      get_legend(p_FragSize),
      get_legend(p_Promoter),
      get_legend(p_splicing),
      get_legend(p_chimeric),
      get_legend(p_Editing),
      NULL,
      get_legend(Mutation_tile),
      get_legend(Chimeric_tile),
      get_legend(clinical_P53),
      get_legend(clinical_MMR),
      get_legend(clinical_HER2),
      get_legend(clinical_Gender),
      get_legend(clinical_Age),
      get_legend(clinical_Cancer),
      ncol = 1,
      rel_heights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
    ggsave(plot=overview,filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Overview_final_20220918.pdf",width = 10,height = 9)
    ggsave(plot=legend,filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/Overview_legend_20220918.pdf",width = 5,height = 30)
   }
  
  #Differential alteration and functional analysis
  {
    #Function
    {
      #commen function for differential analysis: wilcox test
      Wilcox_test <- function(mat_raw,des,output_matrix,output_res,norm_method='NA'){
        #norm_method <- 'NA'
        samples <- des$samples
        #samples <- gsub(".","-",samples,fixed=TRUE)
        #samples <- gsub("X","",samples,fixed=TRUE)
        group <- des$group
        positive <- des[which(des$group=="positive"),]$samples
        negative <- des[which(des$group=="negative"),]$samples
        
        if(norm_method == 'NA' ){
          message('Matrix output without normalization.')
          mat_raw[is.na(mat_raw)] <- 0
          matrix <- mat_raw
        }else{
          message('Matrix output normalized by:',norm_method)
          matrix <- cpm(mat_raw, method=norm_method)
        }
        
        colnames(matrix) <- gsub(".","-",fixed = TRUE,colnames(matrix))
        matrix <- matrix[,as.character(samples)]
        
        
        
        test_func <- function(x){
          positive_mean<-mean(x[group=="positive"])
          negative_mean<-mean(x[group=="negative"])
          if (positive_mean == negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='two.sided')$p.value
          } else if (positive_mean > negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='greater')$p.value
          } else if (positive_mean < negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='less')$p.value
          }
        }
        
        #filter events < 20% samples in minial group
        #positive_prop <- rowSums(matrix[,positive] > 0)/length(positive)
        #negative_prop <- rowSums(matrix[,negative] > 0)/length(negative)
        #matrix <- matrix[(negative_prop>=0.2)&(positive_prop>=0.2),]
        
        #write.table(matrix, output_matrix, sep='\t', quote=FALSE, row.names=TRUE)
        
        pvalues <- apply(matrix, 1, test_func)
        matrix_logcpm = log2(matrix + 1)
        logFC <- apply(matrix_logcpm[,positive], 1, mean) -
          apply(matrix_logcpm[,negative], 1, mean)
        deltaAF <- apply(matrix[,positive], 1, mean) -
          apply(matrix[,negative], 1, mean)
        
        positive_gini <- as.numeric(gini(t(matrix[,positive])))
        negative_gini <- as.numeric(gini(t(matrix[,negative])))
        
        res <- data.frame(log2FoldChange=logFC,
                          deltaAF=deltaAF,
                          pvalue=pvalues, 
                          padj=p.adjust(pvalues, method='BH'),
                          baseMean=apply(matrix, 1, mean),
                          positive_gini = positive_gini,
                          negative_gini = negative_gini)
        write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
      }
      
      #function for RNA expression: edgeR exact test
      edgeR_exact_test <- function(mat_raw,des,output_matrix,output_res){
        samples <- as.character(des$samples)
        #samples <- gsub(".","-",samples,fixed=TRUE)
        #samples <- gsub("X","",samples,fixed=TRUE)
        group <- des$group
        positive <- des[which(des$group=="positive"),]$samples
        negative <- des[which(des$group=="negative"),]$samples
        
        i=1
        mat=as.data.frame(array(dim=c(length(rownames(mat_raw)),1)))
        while (i<=length(samples)) {
          temp <- mat_raw[,which(colnames(mat_raw)==samples[i])]
          #temp <- mat_raw[,grep(samples[i],colnames(mat_raw),fixed=TRUE)]
          #print(i)
          #print(which(colnames(mat_raw)==samples[i]))
          temp <- as.data.frame(temp)
          colnames(temp) <- samples[i]
          mat <- cbind(mat,temp)
          i=i+1
        }
        
        mat <- mat[,-1]
        rownames(mat) <- rownames(mat_raw)
        #filter events < 20% samples in minial group
        positive_prop <- rowSums(mat[,positive] > 0)/length(positive)
        negative_prop <- rowSums(mat[,negative] > 0)/length(negative)
        mat <- mat[(negative_prop>=0.2)&(positive_prop>=0.2),]
        #write.table(mat,output_matrix, sep='\t', quote=FALSE, row.names=TRUE) 
        
        y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
        #keep <- filterByExpr(y,group = group)
        #filterByExpr(min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7, ...)
        #y <- y[keep,]
        y <- calcNormFactors(y, method="TMM")                      #归一化处理 TMM for RNA, CNV and Methylation
        y <- estimateDisp(y)
        test <- exactTest(y, pair = c("negative","positive"), dispersion = "auto")
        #https://support.bioconductor.org/p/64807/
        #广义线性模型可以针对多因素的情况，在这个场景下，只比较cancer与normal用exactTest就足够了
        {
          #design <- model.matrix(~0+group)   
          #y <- estimateDisp(y, design)          #广义线性模型计算离散度（common&trended&tagwise）
          #fit <- glmFit(y, design)              #差异表达分析函数                     
          #test <- glmLRT(fit, coef=2)           #差异表达分析函数 
        }
        
        res <- topTags(test, n=nrow(mat), sort.by='none')     #输出计算差异基因的结果
        res <- cbind(res$table, baseMean=2^(res$table$logCPM)) #输出差异基因结果加上baseMean一栏
        mapped_names <- colnames(res)                   #后面都是整理+输出
        for(i in 1:ncol(res)){
          if(colnames(res)[i] == 'logFC'){
            mapped_names[i] <- 'log2FoldChange'
          }else if(colnames(res)[i] == 'PValue'){
            mapped_names[i] <- 'pvalue'
          }else if(colnames(res)[i] == 'FDR') {
            mapped_names[i] <- 'padj'
          }else{
            mapped_names[i] <- colnames(res)[i]
          }
        }
        colnames(res) <- mapped_names
        positive_gini <- as.numeric(gini(t(mat[rownames(res),positive])))
        negative_gini <- as.numeric(gini(t(mat[rownames(res),negative])))
        res$positive_gini = positive_gini
        res$negative_gini = negative_gini
        write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE) #输出文件，规定名字
      }
      
      #function for RNA expression: edgeR exact test
      edgeR_exact_test_MethylationBin <- function(mat_raw,des,output_matrix,output_res){
        samples <- des$samples
        #samples <- gsub(".","-",samples,fixed=TRUE)
        #samples <- gsub("X","",samples,fixed=TRUE)
        group <- des$group
        positive <- des[which(des$group=="positive"),]$samples
        negative <- des[which(des$group=="negative"),]$samples
        
        i=1
        mat=as.data.frame(array(dim=c(length(rownames(mat_raw)),1)))
        while (i<=length(samples)) {
          temp <- mat_raw[,which(colnames(mat_raw)==samples[i])]
          #temp <- mat_raw[,grep(samples[i],colnames(mat_raw),fixed=TRUE)]
          #print(i)
          #print(which(colnames(mat_raw)==samples[i]))
          temp <- as.data.frame(temp)
          colnames(temp) <- samples[i]
          mat <- cbind(mat,temp)
          i=i+1
        }
        
        mat <- mat[,-1]
        rownames(mat) <- rownames(mat_raw)
        #filter events < 20% samples in minial group
        positive_prop <- rowSums(mat[,positive] > 0)/length(positive)
        negative_prop <- rowSums(mat[,negative] > 0)/length(negative)
        mat <- mat[(negative_prop>=0.2)&(positive_prop>=0.2),]
        #write.table(mat,output_matrix, sep='\t', quote=FALSE, row.names=TRUE) 
        
        y <- DGEList(counts=mat, samples=samples, group=group)     #规定输入格式，必须包含这几个参数
        #keep <- filterByExpr(y,group = group)
        #filterByExpr(min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7, ...)
        #y <- y[keep,]
        y <- calcNormFactors(y, method="TMM")                      #归一化处理 TMM for RNA, CNV and Methylation
        y <- estimateDisp(y)
        test <- exactTest(y, pair = c("negative","positive"), dispersion = "auto")
        #https://support.bioconductor.org/p/64807/
        #广义线性模型可以针对多因素的情况，在这个场景下，只比较cancer与normal用exactTest就足够了
        {
          #design <- model.matrix(~0+group)   
          #y <- estimateDisp(y, design)          #广义线性模型计算离散度（common&trended&tagwise）
          #fit <- glmFit(y, design)              #差异表达分析函数                     
          #test <- glmLRT(fit, coef=2)           #差异表达分析函数 
        }
        
        res <- topTags(test, n=nrow(mat), sort.by='none')     #输出计算差异基因的结果
        res <- cbind(res$table, baseMean=2^(res$table$logCPM)) #输出差异基因结果加上baseMean一栏
        mapped_names <- colnames(res)                   #后面都是整理+输出
        for(i in 1:ncol(res)){
          if(colnames(res)[i] == 'logFC'){
            mapped_names[i] <- 'log2FoldChange'
          }else if(colnames(res)[i] == 'PValue'){
            mapped_names[i] <- 'pvalue'
          }else if(colnames(res)[i] == 'FDR') {
            mapped_names[i] <- 'padj'
          }else{
            mapped_names[i] <- colnames(res)[i]
          }
        }
        colnames(res) <- mapped_names
        positive_gini <- as.numeric(gini(t(mat[rownames(res),positive])))
        negative_gini <- as.numeric(gini(t(mat[rownames(res),negative])))
        res$positive_gini = positive_gini
        res$negative_gini = negative_gini
        write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE) #输出文件，规定名字
      }
      
      #function for RNA mutation: wilcox rank sum test
      AF_wilcox_test <- function(mat_raw,des,output_matrix,output_res,norm_method='NA'){
        #norm_method <- 'NA'
        samples <- des$samples
        #samples <- gsub(".","-",samples,fixed=TRUE)
        #samples <- gsub("X","",samples,fixed=TRUE)
        group <- des$group
        positive <- des[which(des$group=="positive"),]$samples
        negative <- des[which(des$group=="negative"),]$samples
        
        if(norm_method == 'NA' ){
          message('Matrix output without normalization.')
          mat_raw[is.na(mat_raw)] <- 0
          matrix <- mat_raw
        }else{
          message('Matrix output normalized by:',norm_method)
          matrix <- cpm(mat_raw, method=norm_method)
        }
        
        colnames(matrix) <- gsub(".","-",fixed = TRUE,colnames(matrix))
        matrix <- matrix[,as.character(samples)]
        
        test_func <- function(x){
          positive_mean<-mean(x[group=="positive"])
          negative_mean<-mean(x[group=="negative"])
          if (positive_mean == negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='two.sided')$p.value
          } else if (positive_mean > negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='greater')$p.value
          } else if (positive_mean < negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='less')$p.value
          }
        }
        
        #filter events < 20% samples in minial group
        positive_prop <- rowSums(matrix[,positive] > 0)/length(positive)
        negative_prop <- rowSums(matrix[,negative] > 0)/length(negative)
        matrix <- matrix[(negative_prop>=0.2)&(positive_prop>=0.2),]
        #write.table(matrix, output_matrix, sep='\t', quote=FALSE, row.names=TRUE)
        
        positive_gini <- as.numeric(gini(t(matrix[,positive])))
        negative_gini <- as.numeric(gini(t(matrix[,negative])))
        
        pvalues <- apply(matrix, 1, test_func)
        matrix_logcpm = log2(matrix + 1)
        logFC <- apply(matrix_logcpm[,positive], 1, mean) -
          apply(matrix_logcpm[,negative], 1, mean)
        deltaAF <- apply(matrix[,positive], 1, mean) -
          apply(matrix[,negative], 1, mean)
        res <- data.frame(log2FoldChange=logFC,
                          deltaAF=deltaAF,
                          pvalue=pvalues, 
                          padj=p.adjust(pvalues, method='BH'),
                          baseMean=apply(matrix, 1, mean),
                          positive_gini = positive_gini,
                          negative_gini = negative_gini)
        write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
      }
      
      #function for RNA Editing: wilcox rank sum test
      Editing_ratio_wilcox_test <- function(mat_raw,des,output_matrix,output_res,norm_method='NA'){
        samples <- des$samples
        #samples <- gsub(".","-",samples,fixed=TRUE)
        #samples <- gsub("X","",samples,fixed=TRUE)
        group <- des$group
        positive <- des[which(des$group=="positive"),]$samples
        negative <- des[which(des$group=="negative"),]$samples
        
        if(norm_method == 'NA' ){
          message('Matrix output without normalization.')
          mat_raw[is.na(mat_raw)] <- 0
          matrix <- mat_raw
        }else{
          message('Matrix output normalized by:',norm_method)
          matrix <- cpm(mat_raw, method=norm_method)
        }
        
        colnames(matrix) <- gsub(".","-",fixed = TRUE,colnames(matrix))
        matrix <- matrix[,as.character(samples)]
        
        
        
        test_func <- function(x){
          positive_mean<-mean(x[group=="positive"])
          negative_mean<-mean(x[group=="negative"])
          if (positive_mean == negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='two.sided')$p.value
          } else if (positive_mean > negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='greater')$p.value
          } else if (positive_mean < negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='less')$p.value
          }
        }
        
        #filter Editing events < 20% samples in minial group
        positive_prop <- rowSums(matrix[,positive] > 0)/length(positive)
        negative_prop <- rowSums(matrix[,negative] > 0)/length(negative)
        matrix <- matrix[(negative_prop>=0.2)&(positive_prop>=0.2),]
        
        #write.table(matrix, output_matrix, sep='\t', quote=FALSE, row.names=TRUE)
        
        positive_gini <- as.numeric(gini(t(matrix[,positive])))
        negative_gini <- as.numeric(gini(t(matrix[,negative])))
        
        pvalues <- apply(matrix, 1, test_func)
        matrix_logcpm = log2(matrix + 1)
        logFC <- apply(matrix_logcpm[,positive], 1, mean) -
          apply(matrix_logcpm[,negative], 1, mean)
        deltaAF <- apply(matrix[,positive], 1, mean) -
          apply(matrix[,negative], 1, mean)
        res <- data.frame(log2FoldChange=logFC,
                          deltaAF=deltaAF,
                          pvalue=pvalues, 
                          padj=p.adjust(pvalues, method='BH'),
                          baseMean=apply(matrix, 1, mean),
                          positive_gini = positive_gini,
                          negative_gini = negative_gini)
        write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
      }
      
      #function for RNA ASE: wilcox rank sum test
      ASE_ratio_wilcox_test <- function(mat_raw,des,output_matrix,output_res,norm_method='NA'){
        samples <- des$samples
        #samples <- gsub(".","-",samples,fixed=TRUE)
        #samples <- gsub("X","",samples,fixed=TRUE)
        group <- des$group
        positive <- des[which(des$group=="positive"),]$samples
        negative <- des[which(des$group=="negative"),]$samples
        
        if(norm_method == 'NA' ){
          message('Matrix output without normalization.')
          mat_raw[is.na(mat_raw)] <- 0
          matrix <- mat_raw
        }else{
          message('Matrix output normalized by:',norm_method)
          matrix <- cpm(mat_raw, method=norm_method)
        }
        
        colnames(matrix) <- gsub(".","-",fixed = TRUE,colnames(matrix))
        matrix <- matrix[,as.character(samples)]
        
        
        
        test_func <- function(x){
          positive_mean<-mean(x[group=="positive"])
          negative_mean<-mean(x[group=="negative"])
          if (positive_mean == negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='two.sided')$p.value
          } else if (positive_mean > negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='greater')$p.value
          } else if (positive_mean < negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='less')$p.value
          }
        }
        
        #filter Editing events < 20% samples in minial group
        positive_prop <- rowSums(matrix[,positive] > 0)/length(positive)
        negative_prop <- rowSums(matrix[,negative] > 0)/length(negative)
        matrix <- matrix[(negative_prop>=0.2)&(positive_prop>=0.2),]
        
        #write.table(matrix, output_matrix, sep='\t', quote=FALSE, row.names=TRUE)
        
        positive_gini <- as.numeric(gini(t(matrix[,positive])))
        negative_gini <- as.numeric(gini(t(matrix[,negative])))
        
        pvalues <- apply(matrix, 1, test_func)
        matrix_logcpm = log2(matrix + 1)
        logFC <- apply(matrix_logcpm[,positive], 1, mean) -
          apply(matrix_logcpm[,negative], 1, mean)
        deltaAF <- apply(matrix[,positive], 1, mean) -
          apply(matrix[,negative], 1, mean)
        res <- data.frame(log2FoldChange=logFC,
                          deltaAF=deltaAF,
                          pvalue=pvalues, 
                          padj=p.adjust(pvalues, method='BH'),
                          baseMean=apply(matrix, 1, mean),
                          positive_gini = positive_gini,
                          negative_gini = negative_gini)
        write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
      }
      
      #function for chimeric RNA: fisher exact test 
      chimeric_fisher_exact_test <- function(mat_raw,des,output_matrix,output_res){
        norm_method <- 'NA'
        samples <- des$samples
        #samples <- gsub(".","-",samples,fixed=TRUE)
        #samples <- gsub("X","",samples,fixed=TRUE)
        group <- des$group
        positive <- des[which(des$group=="positive"),]$samples
        negative <- des[which(des$group=="negative"),]$samples
        
        if(norm_method == 'NA' ){
          message('Matrix output without normalization.')
          mat_raw[is.na(mat_raw)] <- 0
          matrix <- mat_raw
        }else{
          message('Matrix output normalized by:',norm_method)
          matrix <- cpm(mat_raw, method=norm_method)
        }
        
        colnames(matrix) <- gsub(".","-",fixed = TRUE,colnames(matrix))
        matrix <- matrix[,as.character(samples)]
        
        #filter events < 5% samples in minial group
        positive_prop <- rowSums(matrix[,positive] > 0)/length(positive)
        negative_prop <- rowSums(matrix[,negative] > 0)/length(negative)
        matrix <- matrix[(negative_prop>=0.05)&(positive_prop>=0.05),]
        
        #write.table(matrix, output_matrix, sep='\t', quote=FALSE, row.names=TRUE)
        
        positive_gini <- as.numeric(gini(t(matrix[,positive])))
        negative_gini <- as.numeric(gini(t(matrix[,negative])))
        
        get_zero_number <-function(x) sum(x==0)
        Chimeric_gene_DE <- data.frame(matrix(numeric(0),ncol=4,nrow=nrow(matrix)))
        colnames(Chimeric_gene_DE) <- c("Gene","Positive group chimeric ratio","Negative group chimeric ratio","pvalue")
        i=1
        while(i<=nrow(matrix)){
          Positive_no_chimeric <- apply(matrix[i,positive],1,get_zero_number)
          Negative_no_chimeric <- apply(matrix[i,negative],1,get_zero_number)
          Positive_chimeric <- ncol(matrix[i,positive])-Positive_no_chimeric
          Negative_chimeric <- ncol(matrix[i,negative])-Negative_no_chimeric
          
          fisher_exact_table <- data.frame(positive=c(Positive_chimeric,Positive_no_chimeric),negative = c(Negative_chimeric,Negative_no_chimeric))
          rownames(fisher_exact_table) <- c("Chimeric","No_chimeric_detected")
          fisher_result <- fisher.test(fisher_exact_table)
          Chimeric_gene_DE[i,"pvalue"] <- fisher_result$p.value
          
          Chimeric_gene_DE[i,"Positive group chimeric ratio"] <- Positive_chimeric/(Positive_chimeric+Positive_no_chimeric)
          Chimeric_gene_DE[i,"Negative group chimeric ratio"] <- Negative_chimeric/(Negative_chimeric+Negative_no_chimeric)
          Chimeric_gene_DE[i,"Gene"] <- rownames(matrix[i,])
          i=i+1
        }
        
        deltaFreq <- Chimeric_gene_DE$`Positive group chimeric ratio` -
          Chimeric_gene_DE$`Negative group chimeric ratio`
        res <- data.frame(deltaFreq=deltaFreq,
                          pvalue=Chimeric_gene_DE$pvalue, 
                          padj=p.adjust(Chimeric_gene_DE$pvalue, method='BH'),
                          positive_gini = positive_gini,
                          negative_gini = negative_gini)
        rownames(res) <- Chimeric_gene_DE$Gene
        write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
      }
      
      #function for RNA alternative promoter: wilcox rank sum test
      AlernativePromoter_wilcox_test <- function(mat_raw,des,output_matrix,output_res){
        norm_method <- 'NA'
        samples <- des$samples
        #samples <- gsub(".","-",samples,fixed=TRUE)
        #samples <- gsub("X","",samples,fixed=TRUE)
        group <- des$group
        positive <- des[which(des$group=="positive"),]$samples
        negative <- des[which(des$group=="negative"),]$samples
        
        if(norm_method == 'NA' ){
          message('Matrix output without normalization.')
          mat_raw[is.na(mat_raw)] <- 0
          matrix <- mat_raw
        }else{
          message('Matrix output normalized by:',norm_method)
          matrix <- cpm(mat_raw, method=norm_method)
        }
        
        colnames(matrix) <- gsub(".","-",fixed = TRUE,colnames(matrix))
        matrix <- matrix[,as.character(samples)]
        
        
        
        test_func <- function(x){
          positive_mean<-mean(x[group=="positive"])
          negative_mean<-mean(x[group=="negative"])
          if (positive_mean == negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='two.sided')$p.value
          } else if (positive_mean > negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='greater')$p.value
          } else if (positive_mean < negative_mean) {
            wilcox.test(x=x[group=="positive"], y=x[group=="negative"], alternative='less')$p.value
          }
        }
        
        #filter events < 20% samples in minial group
        positive_prop <- rowSums(matrix[,positive] > 0)/length(positive)
        negative_prop <- rowSums(matrix[,negative] > 0)/length(negative)
        matrix <- matrix[(negative_prop>=0.2)&(positive_prop>=0.2),]
        
        #write.table(matrix, output_matrix, sep='\t', quote=FALSE, row.names=TRUE)
        
        positive_gini <- as.numeric(gini(t(matrix[,positive])))
        negative_gini <- as.numeric(gini(t(matrix[,negative])))
        
        pvalues <- apply(matrix, 1, test_func)
        matrix_logcpm = log2(matrix + 1)
        logFC <- apply(matrix_logcpm[,positive], 1, mean) -
          apply(matrix_logcpm[,negative], 1, mean)
        deltaAF <- apply(matrix[,positive], 1, mean) -
          apply(matrix[,negative], 1, mean)
        res <- data.frame(log2FoldChange=logFC,
                          deltaAF=deltaAF,
                          pvalue=pvalues, 
                          padj=p.adjust(pvalues, method='BH'),
                          baseMean=apply(matrix, 1, mean),
                          positive_gini = positive_gini,
                          negative_gini = negative_gini)
        write.table(res, output_res, sep='\t', quote=FALSE, row.names=TRUE)
      }
      
      #function for enrichment analysis for RNA expression: clusterprofiler
      KEGG_GO_Expression <- function(Differential_result,
                                     output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                     pvalue_cutoff,log2Foldchange_cutoff){
        #all_gene <- row.names(Differential_result)
        #strsplit is highly depend on the colnames of alteration matrix
        #all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
        #background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
        #                    filters = "ensembl_gene_id",
        #                    values=all_gene, mart= mart,useCache = FALSE)
        #background <- background[-which(background$entrezgene_id=="")]
        #background <- background$entrezgene_id
        
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange < -log2Foldchange_cutoff),])
        filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
        forenrich_down <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                                filters = "ensembl_gene_id",
                                values=filtered_down, mart= mart,useCache = FALSE)
        if(length(which(forenrich_down$entrezgene_id==""))==0){
          forenrich_down <- forenrich_down
        } else {
          forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
        }
        forenrich_down <- forenrich_down$entrezgene_id
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange > log2Foldchange_cutoff),])
        filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
        forenrich_up <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                              filters = "ensembl_gene_id",
                              values=filtered_up, mart= mart,useCache = FALSE)
        if(length(which(forenrich_up$entrezgene_id==""))==0){
          forenrich_up <- forenrich_up
        } else {
          forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
        }
        forenrich_up <- forenrich_up$entrezgene_id
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      
      #function for enrichment analysis for miRNA: clusterprofiler
      KEGG_GO_miRNA_target <- function(Differential_result,
                                output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                pvalue_cutoff,log2Foldchange_cutoff){
        #all_tx <- row.names(Differential_result)
        #strsplit is highly depend on the colnames of alteration matrix
        #all_tx <- as.character(lapply(strsplit(all_tx,"\\."),function(x) x[1]))
        #all_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
        #                         filters = "ensembl_transcript_id",
        #                         values=all_tx, mart= mart,useCache = FALSE)
        #all_mir <- all_mir[-which(all_mir$mirbase_id=="")]
        #background <- get_multimir(mirna = all_mir$mirbase_id, summary = TRUE)
        #background <- unique(background@data$target_entrez)
        
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange < -log2Foldchange_cutoff),])
        filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
        down_tx <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
        down_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
                          filters = "ensembl_transcript_id",
                          values=down_tx, mart= mart,useCache = FALSE)
        if(length(which(down_mir$entrezgene_id==""))==0){
          down_mir <- down_mir
        } else {
          down_mir <- down_mir[-which(down_mir$mirbase_id==""),]
        }
        forenrich_down <- get_multimir(mirna = down_mir$mirbase_id, summary = TRUE)
        forenrich_down <- unique(forenrich_down@data$target_entrez)
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange > log2Foldchange_cutoff),])
        filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
        up_tx <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
        up_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
                        filters = "ensembl_transcript_id",
                        values=up_tx, mart= mart,useCache = FALSE)
        if(length(which(up_mir$entrezgene_id==""))==0){
          up_mir <- up_mir
        } else {
          up_mir <- up_mir[-which(up_mir$mirbase_id==""),]
        }
        forenrich_up <- get_multimir(mirna = up_mir$mirbase_id, summary = TRUE)
        forenrich_up <- unique(forenrich_up@data$target_entrez)
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      
      KEGG_GO_miRNA <- function(Differential_result,
                                output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                pvalue_cutoff,log2Foldchange_cutoff){
        #all_tx <- row.names(Differential_result)
        #strsplit is highly depend on the colnames of alteration matrix
        #all_tx <- as.character(lapply(strsplit(all_tx,"\\."),function(x) x[1]))
        #all_mir <- getBM(attributes=c("ensembl_transcript_id", "mirbase_id"),
        #                         filters = "ensembl_transcript_id",
        #                         values=all_tx, mart= mart,useCache = FALSE)
        #all_mir <- all_mir[-which(all_mir$mirbase_id=="")]
        #background <- get_multimir(mirna = all_mir$mirbase_id, summary = TRUE)
        #background <- unique(background@data$target_entrez)
        
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange < -log2Foldchange_cutoff),])
        filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
        down_tx <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
        down_mir <- getBM(attributes=c("ensembl_transcript_id", "entrezgene_id"),
                          filters = "ensembl_transcript_id",
                          values=down_tx, mart= mart,useCache = FALSE)
        forenrich_down <- down_mir$entrezgene_id
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange > log2Foldchange_cutoff),])
        filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
        up_tx <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
        up_mir <- getBM(attributes=c("ensembl_transcript_id", "entrezgene_id"),
                        filters = "ensembl_transcript_id",
                        values=up_tx, mart= mart,useCache = FALSE)
        forenrich_up <- up_mir$entrezgene_id
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      
      #function for enrichment analysis for Splicing: clusterprofiler
      KEGG_GO_Splicing <- function(Differential_result,
                                   output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                   PValue_cutoff,IncLevelDifference_cutoff){
        #all_gene <- row.names(Differential_result)
        #strsplit is highly depend on the colnames of alteration matrix
        #all_gene <- as.character(lapply(strsplit(all_gene,"|",fixed = TRUE),function(x) x[2]))
        #all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
        #background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
        #                    filters = "ensembl_gene_id",
        #                    values=all_gene, mart= mart,useCache = FALSE)
        #background <- background[-which(background$entrezgene_id=="")]
        #background <- background$entrezgene_id
        
        filtered_down <- rownames(Differential_result[(Differential_result$PValue<PValue_cutoff)&(Differential_result$IncLevelDifference < -IncLevelDifference_cutoff),])
        filtered_down <- as.character(lapply(strsplit(filtered_down,"|",fixed = TRUE),function(x) x[2]))
        filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
        forenrich_down <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                                filters = "ensembl_gene_id",
                                values=filtered_down, mart= mart,useCache = FALSE)
        if(length(which(forenrich_down$entrezgene_id==""))==0){
          forenrich_down <- forenrich_down
        } else {
          forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
        }
        forenrich_down <- forenrich_down$entrezgene_id
        
        filtered_up <- rownames(Differential_result[(Differential_result$PValue<PValue_cutoff)&(Differential_result$IncLevelDifference > IncLevelDifference_cutoff),])
        filtered_up <- as.character(lapply(strsplit(filtered_up,"|",fixed = TRUE),function(x) x[2]))
        filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
        forenrich_up <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                              filters = "ensembl_gene_id",
                              values=filtered_up, mart= mart,useCache = FALSE)
        if(length(which(forenrich_up$entrezgene_id==""))==0){
          forenrich_up <- forenrich_up
        } else {
          forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
        }
        forenrich_up <- forenrich_up$entrezgene_id
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      
      #function for enrichment analysis for Mutation site: clusterprofiler
      KEGG_GO_Mutation <- function(Differential_result,
                                   output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                   pvalue_cutoff,deltaAF_cutoff){
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF < -deltaAF_cutoff),])
        chr <- unlist(lapply(strsplit(filtered_down,"_"), function(x) x[1]))
        chr <- gsub("chr","",chr)
        start <- as.integer(unlist(lapply(strsplit(filtered_down,"_"), function(x) x[2])))
        if(length(which(is.na(start)))==0){
          chr <- chr
          start <- start
        } else {
          chr <- chr[-which(is.na(start))]
          start <- start[-which(is.na(start))]
        }
        
        i=1
        forenrich_down <- data.frame()
        while(i<=length(chr)){
          message("Processing: ",i,"/",length(chr))
          attributes <- c("ensembl_gene_id","start_position","end_position","entrezgene_id")
          filters <- c("chromosome_name","start","end")
          values <- list(chromosome=chr[i],start=start[i],end=start[i])
          forenrich_down_tmp <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id"),
                                      filters=c("chromosome_name","start","end"), 
                                      values=values, mart=mart, useCache = FALSE)
          forenrich_down <- rbind(forenrich_down,forenrich_down_tmp)
          i=i+1
        }
        
        forenrich_down <- na.omit(forenrich_down)
        if(length(which(forenrich_down$entrezgene_id==""))==0){
          forenrich_down <- forenrich_down
        } else {
          forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
        }
        forenrich_down <- forenrich_down$entrezgene_id
        forenrich_down <- unique(forenrich_down)
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF > deltaAF_cutoff),])
        chr <- unlist(lapply(strsplit(filtered_up,"_"), function(x) x[1]))
        chr <- gsub("chr","",chr)
        start <- as.integer(unlist(lapply(strsplit(filtered_up,"_"), function(x) x[2])))
        if(length(which(is.na(start)))==0){
          chr <- chr
          start <- start
        } else {
          chr <- chr[-which(is.na(start))]
          start <- start[-which(is.na(start))]
        }
        
        i=1
        forenrich_up <- data.frame()
        while(i<=length(chr)){
          message("Processing: ",i,"/",length(chr))
          attributes <- c("ensembl_gene_id","start_position","end_position","entrezgene_id")
          filters <- c("chromosome_name","start","end")
          values <- list(chromosome=chr[i],start=start[i],end=start[i])
          forenrich_up_tmp <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id"),
                                    filters=c("chromosome_name","start","end"), 
                                    values=values, mart=mart, useCache = FALSE)
          forenrich_up <- rbind(forenrich_up,forenrich_up_tmp)
          i=i+1
        }
        
        forenrich_up <- na.omit(forenrich_up)
        if(length(which(forenrich_up$entrezgene_id==""))==0){
          forenrich_up <- forenrich_up
        } else {
          forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
        }
        forenrich_up <- forenrich_up$entrezgene_id
        forenrich_up <- unique(forenrich_up)
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      
      #function for enrichment analysis for Editing site: clusterprofiler
      KEGG_GO_Editing <- function(Differential_result,
                                  output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                  pvalue_cutoff,deltaAF_cutoff){
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF < -deltaAF_cutoff),])
        chr <- unlist(lapply(strsplit(filtered_down,"_"), function(x) x[1]))
        chr <- gsub("chr","",chr)
        start <- as.integer(unlist(lapply(strsplit(filtered_down,"_"), function(x) x[2])))
        if(length(which(is.na(start)))==0){
          chr <- chr
          start <- start
        } else {
          chr <- chr[-which(is.na(start))]
          start <- start[-which(is.na(start))]
        }
        
        i=1
        forenrich_down <- data.frame()
        while(i<=length(chr)){
          message("Processing: ",i,"/",length(chr))
          attributes <- c("ensembl_gene_id","start_position","end_position","entrezgene_id")
          filters <- c("chromosome_name","start","end")
          values <- list(chromosome=chr[i],start=start[i],end=start[i])
          forenrich_down_tmp <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id"),
                                      filters=c("chromosome_name","start","end"), 
                                      values=values, mart=mart, useCache = FALSE)
          forenrich_down <- rbind(forenrich_down,forenrich_down_tmp)
          i=i+1
        }
        
        forenrich_down <- na.omit(forenrich_down)
        if(length(which(forenrich_down$entrezgene_id==""))==0){
          forenrich_down <- forenrich_down
        } else {
          forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
        }
        forenrich_down <- forenrich_down$entrezgene_id
        forenrich_down <- unique(forenrich_down)
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF > deltaAF_cutoff),])
        chr <- unlist(lapply(strsplit(filtered_up,"_"), function(x) x[1]))
        chr <- gsub("chr","",chr)
        start <- as.integer(unlist(lapply(strsplit(filtered_up,"_"), function(x) x[2])))
        if(length(which(is.na(start)))==0){
          chr <- chr
          start <- start
        } else {
          chr <- chr[-which(is.na(start))]
          start <- start[-which(is.na(start))]
        }
        
        i=1
        forenrich_up <- data.frame()
        while(i<=length(chr)){
          message("Processing: ",i,"/",length(chr))
          attributes <- c("ensembl_gene_id","start_position","end_position","entrezgene_id")
          filters <- c("chromosome_name","start","end")
          values <- list(chromosome=chr[i],start=start[i],end=start[i])
          forenrich_up_tmp <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id"),
                                    filters=c("chromosome_name","start","end"), 
                                    values=values, mart=mart, useCache = FALSE)
          forenrich_up <- rbind(forenrich_up,forenrich_up_tmp)
          i=i+1
        }
        
        forenrich_up <- na.omit(forenrich_up)
        if(length(which(forenrich_up$entrezgene_id==""))==0){
          forenrich_up <- forenrich_up
        } else {
          forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
        }
        forenrich_up <- forenrich_up$entrezgene_id
        forenrich_up <- unique(forenrich_up)
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      #function for enrichment analysis for ASE site: clusterprofiler
      KEGG_GO_ASE <- function(Differential_result,
                              output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                              pvalue_cutoff,deltaAF_cutoff){
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF < -deltaAF_cutoff),])
        chr <- unlist(lapply(strsplit(filtered_down,"_",fixed=TRUE), function(x) x[1]))
        chr <- gsub("chr","",chr)
        start <- as.integer(unlist(lapply(strsplit(filtered_down,"_",fixed=TRUE), function(x) x[2])))
        if(length(which(is.na(start)))==0){
          chr <- chr
          start <- start
        } else {
          chr <- chr[-which(is.na(start))]
          start <- start[-which(is.na(start))]
        }
        
        i=1
        forenrich_down <- data.frame()
        while(i<=length(chr)){
          message("Processing: ",i,"/",length(chr))
          attributes <- c("ensembl_gene_id","start_position","end_position","entrezgene_id")
          filters <- c("chromosome_name","start","end")
          values <- list(chromosome=chr[i],start=start[i],end=start[i])
          forenrich_down_tmp <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id"),
                                      filters=c("chromosome_name","start","end"), 
                                      values=values, mart=mart, useCache = FALSE)
          forenrich_down <- rbind(forenrich_down,forenrich_down_tmp)
          i=i+1
        }
        
        forenrich_down <- na.omit(forenrich_down)
        if(length(which(forenrich_down$entrezgene_id==""))==0){
          forenrich_down <- forenrich_down
        } else {
          forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
        }
        forenrich_down <- forenrich_down$entrezgene_id
        forenrich_down <- unique(forenrich_down)
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF > deltaAF_cutoff),])
        chr <- unlist(lapply(strsplit(filtered_up,"|",fixed=TRUE), function(x) x[1]))
        chr <- gsub("chr","",chr)
        start <- as.integer(unlist(lapply(strsplit(filtered_up,"|",fixed=TRUE), function(x) x[2])))
        if(length(which(is.na(start)))==0){
          chr <- chr
          start <- start
        } else {
          chr <- chr[-which(is.na(start))]
          start <- start[-which(is.na(start))]
        }
        
        i=1
        forenrich_up <- data.frame()
        while(i<=length(chr)){
          message("Processing: ",i,"/",length(chr))
          attributes <- c("ensembl_gene_id","start_position","end_position","entrezgene_id")
          filters <- c("chromosome_name","start","end")
          values <- list(chromosome=chr[i],start=start[i],end=start[i])
          forenrich_up_tmp <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","entrezgene_id"),
                                    filters=c("chromosome_name","start","end"), 
                                    values=values, mart=mart, useCache = FALSE)
          forenrich_up <- rbind(forenrich_up,forenrich_up_tmp)
          i=i+1
        }
        
        forenrich_up <- na.omit(forenrich_up)
        if(length(which(forenrich_up$entrezgene_id==""))==0){
          forenrich_up <- forenrich_up
        } else {
          forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
        }
        forenrich_up <- forenrich_up$entrezgene_id
        forenrich_up <- unique(forenrich_up)
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      
      #function for enrichment analysis for chimeric RNA: clusterprofiler
      KEGG_GO_chimeric <- function(Differential_result,
                                   output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                   pvalue_cutoff,deltaFreq_cutoff){
        #all_gene <- row.names(Differential_result)
        #strsplit is highly depend on the colnames of alteration matrix
        #all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
        #background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
        #                    filters = "ensembl_gene_id",
        #                    values=all_gene, mart= mart,useCache = FALSE)
        #background <- background[-which(background$entrezgene_id=="")]
        #background <- background$entrezgene_id
        
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaFreq < -deltaFreq_cutoff),])
        filtered_down_head <- as.character(lapply(strsplit(filtered_down,"\\|"),function(x) x[2]))
        filtered_down_head <- as.character(lapply(strsplit(filtered_down_head,"^"),function(x) x[2]))
        filtered_down_head <- as.character(lapply(strsplit(filtered_down_head,".",fixed = TRUE),function(x) x[1]))
        filtered_down_tail <- as.character(lapply(strsplit(filtered_down,"\\|"),function(x) x[3]))
        filtered_down_tail <- as.character(lapply(strsplit(filtered_down_tail,"^"),function(x) x[2]))
        filtered_down_tail <- as.character(lapply(strsplit(filtered_down_tail,".",fixed = TRUE),function(x) x[1]))
        filtered_down <- fct_c(as.factor(filtered_down_head),as.factor(filtered_down_tail))
        forenrich_down <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                                filters = "ensembl_gene_id",
                                values=filtered_down, mart= mart,useCache = FALSE)
        if(length(which(forenrich_down$entrezgene_id==""))==0){
          forenrich_down <- forenrich_down
        } else {
          forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
        }
        forenrich_down <- forenrich_down$entrezgene_id
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaFreq > deltaFreq_cutoff),])
        filtered_up_head <- as.character(lapply(strsplit(filtered_up,"\\|"),function(x) x[2]))
        filtered_up_head <- as.character(lapply(strsplit(filtered_up_head,"^"),function(x) x[2]))
        filtered_up_head <- as.character(lapply(strsplit(filtered_up_head,".",fixed = TRUE),function(x) x[1]))
        filtered_up_tail <- as.character(lapply(strsplit(filtered_up,"\\|"),function(x) x[3]))
        filtered_up_tail <- as.character(lapply(strsplit(filtered_up_tail,"^"),function(x) x[2]))
        filtered_up_tail <- as.character(lapply(strsplit(filtered_up_tail,".",fixed = TRUE),function(x) x[1]))
        filtered_up <- fct_c(as.factor(filtered_up_head),as.factor(filtered_up_tail))
        forenrich_up <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                              filters = "ensembl_gene_id",
                              values=filtered_up, mart= mart,useCache = FALSE)
        if(length(which(forenrich_up$entrezgene_id==""))==0){
          forenrich_up <- forenrich_up
        } else {
          forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
        }
        forenrich_up <- forenrich_up$entrezgene_id
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      
      #function for enrichment analysis for RNA expression: clusterprofiler
      KEGG_GO_AlternativePromoter <- function(Differential_result,
                                              output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                              pvalue_cutoff,deltaAF_cutoff){
        #all_gene <- row.names(Differential_result)
        #strsplit is highly depend on the colnames of alteration matrix
        #all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
        #background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
        #                    filters = "ensembl_gene_id",
        #                    values=all_gene, mart= mart,useCache = FALSE)
        #background <- background[-which(background$entrezgene_id=="")]
        #background <- background$entrezgene_id
        
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF < -deltaAF_cutoff),])
        filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
        forenrich_down <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                                filters = "ensembl_gene_id",
                                values=filtered_down, mart= mart,useCache = FALSE)
        if(length(which(forenrich_down$entrezgene_id==""))==0){
          forenrich_down <- as.character(na.omit(forenrich_down$entrezgene_id))
        } else {
          forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
          forenrich_down <- as.character(na.omit(forenrich_down))
        }
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF > deltaAF_cutoff),])
        filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
        forenrich_up <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                              filters = "ensembl_gene_id",
                              values=filtered_up, mart= mart,useCache = FALSE)
        if(length(which(forenrich_up$entrezgene_id==""))==0){
          forenrich_up <- as.character(na.omit(forenrich_up$entrezgene_id))
        } else {
          forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
          forenrich_up <- as.character(na.omit(forenrich_up))
        }
        
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      
      #function for enrichment analysis for DNA methylation: clusterprofiler
      KEGG_GO_Methylation <- function(Differential_result,
                                      output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                      pvalue_cutoff,log2Foldchange_cutoff){
        #all_gene <- row.names(Differential_result)
        #strsplit is highly depend on the colnames of alteration matrix
        #all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
        #background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
        #                    filters = "ensembl_gene_id",
        #                    values=all_gene, mart= mart,useCache = FALSE)
        #background <- background[-which(background$entrezgene_id=="")]
        #background <- background$entrezgene_id
        
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange < -log2Foldchange_cutoff),])
        filtered_down <- gsub("promoter_","",filtered_down)
        filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
        forenrich_down <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                                filters = "ensembl_gene_id",
                                values=filtered_down, mart= mart,useCache = FALSE)
        if(length(which(forenrich_down$entrezgene_id==""))==0){
          forenrich_down <- forenrich_down
        } else {
          forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
        }
        forenrich_down <- forenrich_down$entrezgene_id
        forenrich_down <- as.character(na.omit(forenrich_down))
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange > log2Foldchange_cutoff),])
        filtered_up <- gsub("promoter_","",filtered_up)
        filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
        forenrich_up <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                              filters = "ensembl_gene_id",
                              values=filtered_up, mart= mart,useCache = FALSE)
        if(length(which(forenrich_up$entrezgene_id==""))==0){
          forenrich_up <- forenrich_up
        } else {
          forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
        }
        forenrich_up <- forenrich_up$entrezgene_id
        forenrich_up <- as.character(na.omit(forenrich_up))
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      
      #function for enrichment analysis for DNA WPS: clusterprofiler
      KEGG_GO_WPS <- function(Differential_result,
                                     output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                     pvalue_cutoff,deltaAF_cutoff){
        #all_gene <- row.names(Differential_result)
        #strsplit is highly depend on the colnames of alteration matrix
        #all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
        #background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
        #                    filters = "ensembl_gene_id",
        #                    values=all_gene, mart= mart,useCache = FALSE)
        #background <- background[-which(background$entrezgene_id=="")]
        #background <- background$entrezgene_id
        
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF < -deltaAF_cutoff),])
        #filtered_down <- gsub("NOVsum_","",filtered_down)
        #filtered_down <- gsub("promoter300100exon1end_","",filtered_down)
        filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
        forenrich_down <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                                filters = "ensembl_gene_id",
                                values=filtered_down, mart= mart,useCache = FALSE)
        if(length(which(forenrich_down$entrezgene_id==""))==0){
          forenrich_down <- forenrich_down
        } else {
          forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
        }
        forenrich_down <- forenrich_down$entrezgene_id
        forenrich_down <- as.character(na.omit(forenrich_down))
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$deltaAF > deltaAF_cutoff),])

        filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
        forenrich_up <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                              filters = "ensembl_gene_id",
                              values=filtered_up, mart= mart,useCache = FALSE)
        if(length(which(forenrich_up$entrezgene_id==""))==0){
          forenrich_up <- forenrich_up
        } else {
          forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
        }
        forenrich_up <- forenrich_up$entrezgene_id
        forenrich_up <- as.character(na.omit(forenrich_up))
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
      
      #function for enrichment analysis for DNA copy number: clusterprofiler
      KEGG_GO_CNV <- function(Differential_result,
                              output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                              pvalue_cutoff,log2Foldchange_cutoff){
        #all_gene <- row.names(Differential_result)
        #strsplit is highly depend on the colnames of alteration matrix
        #all_gene <- as.character(lapply(strsplit(all_gene,"\\."),function(x) x[1]))
        #background <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
        #                    filters = "ensembl_gene_id",
        #                    values=all_gene, mart= mart,useCache = FALSE)
        #background <- background[-which(background$entrezgene_id=="")]
        #background <- background$entrezgene_id
        
        filtered_down <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange < -log2Foldchange_cutoff),])
        filtered_down <- as.character(lapply(strsplit(filtered_down,"\\."),function(x) x[1]))
        forenrich_down <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                                filters = "ensembl_gene_id",
                                values=filtered_down, mart= mart,useCache = FALSE)
        if(length(which(forenrich_down$entrezgene_id==""))==0){
          forenrich_down <- forenrich_down
        } else {
          forenrich_down <- forenrich_down[-which(forenrich_down$entrezgene_id==""),]
        }
        forenrich_down <- forenrich_down$entrezgene_id
        forenrich_down <- as.character(na.omit(forenrich_down))
        
        filtered_up <- rownames(Differential_result[(Differential_result$pvalue<pvalue_cutoff)&(Differential_result$log2FoldChange > log2Foldchange_cutoff),])
        filtered_up <- gsub("NOVsum_","",filtered_up)
        filtered_up <- as.character(lapply(strsplit(filtered_up,"\\."),function(x) x[1]))
        forenrich_up <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                              filters = "ensembl_gene_id",
                              values=filtered_up, mart= mart,useCache = FALSE)
        if(length(which(forenrich_up$entrezgene_id==""))==0){
          forenrich_up <- forenrich_up
        } else {
          forenrich_up <- forenrich_up[-which(forenrich_up$entrezgene_id==""),]
        }
        forenrich_up <- forenrich_up$entrezgene_id
        forenrich_up <- as.character(na.omit(forenrich_up))
        
        #KEGG
        {
          KEGG_res_down <- enrichKEGG(
            forenrich_down,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_down <- KEGG_res_down@result
          KEGG_output_down$GeneEnrichedIn <- "Down regulated"
          
          KEGG_res_up <- enrichKEGG(
            forenrich_up,
            organism = "hsa",
            keyType = "kegg",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1,
            use_internal_data = FALSE)
          KEGG_output_up <- KEGG_res_up@result
          KEGG_output_up$GeneEnrichedIn <- "Up regulated"
          
          KEGG_output <- rbind(KEGG_output_up,KEGG_output_down)
          write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
        }
        #GO_BP
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "BP",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_BP,quote = FALSE,sep = "\t")
        }
        #GO_CC
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "CC",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_CC,quote = FALSE,sep = "\t")
        }
        #GO_MF
        {
          library(clusterProfiler)
          GO_res_down <- enrichGO(
            forenrich_down,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_down <- GO_res_down@result
          GO_output_down$GeneEnrichedIn <- "Down regulated"
          
          GO_res_up <- enrichGO(
            forenrich_up,
            'org.Hs.eg.db',
            ont = "MF",
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            #universe=as.character(background),
            minGSSize = 0,
            maxGSSize = 500,
            qvalueCutoff = 1)
          GO_output_up <- GO_res_up@result
          GO_output_up$GeneEnrichedIn <- "Up regulated"
          
          GO_output <- rbind(GO_output_up,GO_output_down)
          
          write.table(GO_output,output_GO_MF,quote = FALSE,sep = "\t")
        }
      }
    }
    
    #Differential analysis and enrichment analysis (wilcox)
    {
      setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 3")
      ##Expression
      {
      #RNA expression matrix: /data/taoyuhuan/projects/exOmics_RNA/multiomics_paired/output/Intron-spanning/count/count_matrix/featurecount_intron_spanning.txt
      #Differential expression 
      {
        #mat_raw <- read.table("./Expression/matrix/featurecount_intron_spanning.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #in molecular cancer paper
        mat_raw <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/Plasma_featurecount_intron_spanning.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
        #mat_raw <- mat_raw[grep("ENSG",rownames(mat_raw)),]
        
        des <- read.csv("./Expression/group/des_Expression_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Expression/matrix/Expression_for_GIvsNC"
        output_res <- "./Expression/output/Expression_GIvsNC_edger_exact.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Expression/group/des_Expression_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Expression/matrix/Expression_for_CRCvsNC"
        output_res <- "./Expression/output/Expression_CRCvsNC_edger_exact.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Expression/group/des_Expression_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Expression/matrix/Expression_for_STADvsNC"
        output_res <- "./Expression/output/Expression_STADvsNC_edger_exact.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Expression/group/des_Expression_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Expression/matrix/Expression_for_STADvsCRC"
        output_res <- "./Expression/output/Expression_STADvsCRC_edger_exact.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
      }
      
      #Functional analysis
      {
        #Enrichment analysis: gastrointestinal cancer vs healthy donors
        {
          Differential_result <- read.csv("./Expression/output/Expression_GIvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Expression/output/Expression_GIvsNC_KEGG.txt"
          output_GO_BP <- "./Expression/output/Expression_GIvsNC_GO_BP.txt"
          output_GO_CC <- "./Expression/output/Expression_GIvsNC_GO_CC.txt"
          output_GO_MF <- "./Expression/output/Expression_GIvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Expression(Differential_result,
                             output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                             pvalue_cutoff,log2Foldchange_cutoff)
        }
        #Enrichment analysis: STAD vs CRC
        {
          Differential_result <- read.csv("./Expression/output/Expression_STADvsCRC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Expression/output/Expression_STADvsCRC_KEGG.txt"
          output_GO_BP <- "./Expression/output/Expression_STADvsCRC_GO_BP.txt"
          output_GO_CC <- "./Expression/output/Expression_STADvsCRC_GO_CC.txt"
          output_GO_MF <- "./Expression/output/Expression_STADvsCRC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Expression(Differential_result,
                             output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                             pvalue_cutoff,log2Foldchange_cutoff)
        }
        #Enrichment analysis: STAD vs NC
        {
          Differential_result <- read.csv("./Expression/output/Expression_STADvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Expression/output/Expression_STADvsNC_KEGG.txt"
          output_GO_BP <- "./Expression/output/Expression_STADvsNC_GO_BP.txt"
          output_GO_CC <- "./Expression/output/Expression_STADvsNC_GO_CC.txt"
          output_GO_MF <- "./Expression/output/Expression_STADvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Expression(Differential_result,
                             output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                             pvalue_cutoff,log2Foldchange_cutoff)
        }
        #Enrichment analysis: CRC vs NC
        {
          Differential_result <- read.csv("./Expression/output/Expression_CRCvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Expression/output/Expression_CRCvsNC_KEGG.txt"
          output_GO_BP <- "./Expression/output/Expression_CRCvsNC_GO_BP.txt"
          output_GO_CC <- "./Expression/output/Expression_CRCvsNC_GO_CC.txt"
          output_GO_MF <- "./Expression/output/Expression_CRCvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Expression(Differential_result,
                             output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                             pvalue_cutoff,log2Foldchange_cutoff)
        }
      }
      }
      
      ##miRNA
      {
      #Differential Analysis: gastrointestinal cancer vs healthy donors
      {
        mat_raw <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/miRNA.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
        
        des <- read.csv("./miRNA/group/des_miRNA_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./miRNA/matrix/miRNA_for_GIvsNC.txt"
        output_res <- "./miRNA/output/miRNA_GIvsNC_edger_exact.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./miRNA/group/des_miRNA_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./miRNA/matrix/miRNA_for_CRCvsNC.txt"
        output_res <- "./miRNA/output/miRNA_CRCvsNC_edger_exact.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./miRNA/group/des_miRNA_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./miRNA/matrix/miRNA_for_STADvsNC.txt"
        output_res <- "./miRNA/output/miRNA_STADvsNC_edger_exact.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./miRNA/group/des_miRNA_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./miRNA/matrix/miRNA_for_STADvsCRC.txt"
        output_res <- "./miRNA/output/miRNA_STADvsCRC_edger_exact.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
      }
      #Functional analysis
      {
        #Enrichment analysis: gastrointestinal cancer vs healthy donors
        {
          Differential_result <- read.csv("./miRNA/output/miRNA_GIvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./miRNA/output/miRNA_GIvsNC_KEGG.txt"
          output_GO_BP <- "./miRNA/output/miRNA_GIvsNC_GO_BP.txt"
          output_GO_CC <- "./miRNA/output/miRNA_GIvsNC_GO_CC.txt"
          output_GO_MF <- "./miRNA/output/miRNA_GIvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_miRNA(Differential_result,
                        output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                        pvalue_cutoff,log2Foldchange_cutoff)
        }
        #Enrichment analysis: STAD vs CRC
        {
          Differential_result <- read.csv("./miRNA/output/miRNA_STADvsCRC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./miRNA/output/miRNA_STADvsCRC_KEGG.txt"
          output_GO_BP <- "./miRNA/output/miRNA_STADvsCRC_GO_BP.txt"
          output_GO_CC <- "./miRNA/output/miRNA_STADvsCRC_GO_CC.txt"
          output_GO_MF <- "./miRNA/output/miRNA_STADvsCRC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_miRNA(Differential_result,
                        output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                        pvalue_cutoff,log2Foldchange_cutoff)
        }
        #Enrichment analysis: STAD vs NC
        {
          Differential_result <- read.csv("./miRNA/output/miRNA_STADvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./miRNA/output/miRNA_STADvsNC_KEGG.txt"
          output_GO_BP <- "./miRNA/output/miRNA_STADvsNC_GO_BP.txt"
          output_GO_CC <- "./miRNA/output/miRNA_STADvsNC_GO_CC.txt"
          output_GO_MF <- "./miRNA/output/miRNA_STADvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_miRNA(Differential_result,
                        output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                        pvalue_cutoff,log2Foldchange_cutoff)
        }
        #Enrichment analysis: CRC vs NC
        {
          Differential_result <- read.csv("./miRNA/output/miRNA_CRCvsNC_edger_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./miRNA/output/miRNA_CRCvsNC_KEGG.txt"
          output_GO_BP <- "./miRNA/output/miRNA_CRCvsNC_GO_BP.txt"
          output_GO_CC <- "./miRNA/output/miRNA_CRCvsNC_GO_CC.txt"
          output_GO_MF <- "./miRNA/output/miRNA_CRCvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_miRNA(Differential_result,
                        output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                        pvalue_cutoff,log2Foldchange_cutoff)
        }
      }
      }
    
      ##Splicing
      {
      #splicing output: /data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/multiomics_paired/allpassed-20211020
      #example inclevel matrix: /data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/multiomics_paired/allpassed-20211020/CRCvsNC/summary/Splicing_CRCvsNC_allpassed.txt
      #Differential splicing done by rMATs
      #example differential splicing: /data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/multiomics_paired/allpassed-20211020/CRCvsNC/summary/Splicing_differential_CRCvsNC_allpassed.txt
      #readin stats and inclevel(for example, Rscript is located in /data/taoyuhuan/projects/exOmics_RNA/level_3_splicing)
      {
        #setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/")
        #file_tosum <- "./Splicing/output"
        #files <- dir(file_tosum)
        #stats <- grep("stats",files,value = TRUE)
        #inclevels <- grep("inc_level",files,value = TRUE)
        #output_stat <- "./Splicing/output/Splicing_GIvsNC_rMATs.txt"
        #output_inclevel <- "./Splicing/matrix/Inclevel_for_GIvsNC.txt"
        
        #stat_all <- data.frame()
        #for(i in stats){
        #  stat <- read.csv(paste0(file_tosum,"/",i),sep = "\t", header = TRUE, row.names = 1)
        #  stat_all <- rbind(stat_all,stat)
        #}
        #write.table(stat_all,output_stat,quote = FALSE, sep = "\t")
        
        #inclevel_matrix <- data.frame()
        #for(i in inclevels){
        #  inclevel <- read.csv(paste0(file_tosum,"/",i),sep = "\t", header = TRUE, row.names = 1)
        #  inclevel_matrix <- rbind(inclevel_matrix,inclevel)
        #}
        #write.table(inclevel_matrix,output_inclevel,quote = FALSE, sep = "\t")
      }
      
      #Functional
      {
        #Enrichment analysis: gastrointestinal cancer vs healthy donors
        {
          Differential_result <- read.csv("./Splicing/output/multiomics_paired_Cancer_vs_Healthy.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Splicing/output/Splicing_GIvsNC_KEGG.txt"
          output_GO_BP <- "./Splicing/output/Splicing_GIvsNC_GO_BP.txt"
          output_GO_CC <- "./Splicing/output/Splicing_GIvsNC_GO_CC.txt"
          output_GO_MF <- "./Splicing/output/Splicing_GIvsNC_GO_MF.txt"
          PValue_cutoff <- 0.05
          IncLevelDifference_cutoff <- 0.05
          
          KEGG_GO_Splicing(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           PValue_cutoff,IncLevelDifference_cutoff)
        }
        #Enrichment analysis: CRC vs NC
        {
          Differential_result <- read.csv("./Splicing/output/multiomics_paired_CRC_vs_Healthy.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Splicing/output/Splicing_CRCvsNC_KEGG.txt"
          output_GO_BP <- "./Splicing/output/Splicing_CRCvsNC_GO_BP.txt"
          output_GO_CC <- "./Splicing/output/Splicing_CRCvsNC_GO_CC.txt"
          output_GO_MF <- "./Splicing/output/Splicing_CRCvsNC_GO_MF.txt"
          PValue_cutoff <- 0.05
          IncLevelDifference_cutoff <- 0.05
          
          KEGG_GO_Splicing(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           PValue_cutoff,IncLevelDifference_cutoff)
        }
        #Enrichment analysis: STAD vs NC
        {
          Differential_result <- read.csv("./Splicing/output/multiomics_paired_STAD_vs_Healthy.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Splicing/output/Splicing_STADvsNC_KEGG.txt"
          output_GO_BP <- "./Splicing/output/Splicing_STADvsNC_GO_BP.txt"
          output_GO_CC <- "./Splicing/output/Splicing_STADvsNC_GO_CC.txt"
          output_GO_MF <- "./Splicing/output/Splicing_STADvsNC_GO_MF.txt"
          PValue_cutoff <- 0.05
          IncLevelDifference_cutoff <- 0.05
          
          KEGG_GO_Splicing(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           PValue_cutoff,IncLevelDifference_cutoff)
        }
        #Enrichment analysis: STAD vs CRC
        {
          Differential_result <- read.csv("./Splicing/output/multiomics_paired_STAD_vs_CRC.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Splicing/output/Splicing_STADvsCRC_KEGG.txt"
          output_GO_BP <- "./Splicing/output/Splicing_STADvsCRC_GO_BP.txt"
          output_GO_CC <- "./Splicing/output/Splicing_STADvsCRC_GO_CC.txt"
          output_GO_MF <- "./Splicing/output/Splicing_STADvsCRC_GO_MF.txt"
          PValue_cutoff <- 0.05
          IncLevelDifference_cutoff <- 0.05
          
          KEGG_GO_Splicing(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           PValue_cutoff,IncLevelDifference_cutoff)
        }
      }
      }
      
      ##mutation
      {
      #mutation output (variant allele fraction per site):/data/taoyuhuan/projects/Machine_learning/20210722_multiomics/SNP_20211011/Matrix/final_SNP_AF.txt
      #Differential mutation analysis by allele fraction per site by wilcox rank sum test
      {
        mat_raw <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/final_SNP_AF.txt",sep = "\t",header = TRUE, row.names = 1)
        des <- read.csv("./Mutation/group/des_Mutation_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Mutation/matrix/Mutation_site_for_GIvsNC.txt"
        output_res <- "./Mutation/output/Mutation_site_GIvsNC_wilcox.txt"
        AF_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Mutation/group/des_Mutation_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Mutation/matrix/Mutation_site_for_CRCvsNC.txt"
        output_res <- "./Mutation/output/Mutation_site_CRCvsNC_wilcox.txt"
        AF_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Mutation/group/des_Mutation_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Mutation/matrix/Mutation_site_for_STADvsNC.txt"
        output_res <- "./Mutation/output/Mutation_site_STADvsNC_wilcox.txt"
        AF_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Mutation/group/des_Mutation_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Mutation/matrix/Mutation_site_for_STADvsCRC.txt"
        output_res <- "./Mutation/output/Mutation_site_STADvsCRC_wilcox.txt"
        AF_wilcox_test(mat_raw,des,output_matrix,output_res)
      }
      
      #Functional
      {
        #Enrichment analysis: gastrointestinal cancer vs healthy donors
        {
          Differential_result <- read.csv("./Mutation/output/Mutation_site_GIvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Mutation/output/Mutation_GIvsNC_KEGG.txt"
          output_GO_BP <- "./Mutation/output/Mutation_GIvsNC_GO_BP.txt"
          output_GO_CC <- "./Mutation/output/Mutation_GIvsNC_GO_CC.txt"
          output_GO_MF <- "./Mutation/output/Mutation_GIvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_Mutation(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaAF_cutoff)
        }
        #Enrichment analysis: CRC vs NC
        {
          Differential_result <- read.csv("./Mutation/output/Mutation_site_CRCvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Mutation/output/Mutation_CRCvsNC_KEGG.txt"
          output_GO_BP <- "./Mutation/output/Mutation_CRCvsNC_GO_BP.txt"
          output_GO_CC <- "./Mutation/output/Mutation_CRCvsNC_GO_CC.txt"
          output_GO_MF <- "./Mutation/output/Mutation_CRCvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_Mutation(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaAF_cutoff)
        }
        #Enrichment analysis: STAD vs NC
        {
          Differential_result <- read.csv("./Mutation/output/Mutation_site_STADvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Mutation/output/Mutation_STADvsNC_KEGG.txt"
          output_GO_BP <- "./Mutation/output/Mutation_STADvsNC_GO_BP.txt"
          output_GO_CC <- "./Mutation/output/Mutation_STADvsNC_GO_CC.txt"
          output_GO_MF <- "./Mutation/output/Mutation_STADvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_Mutation(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaAF_cutoff)
        }
        #Enrichment analysis: STAD vs CRC
        {
          Differential_result <- read.csv("./Mutation/output/Mutation_site_STADvsCRC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Mutation/output/Mutation_STADvsCRC_KEGG.txt"
          output_GO_BP <- "./Mutation/output/Mutation_STADvsCRC_GO_BP.txt"
          output_GO_CC <- "./Mutation/output/Mutation_STADvsCRC_GO_CC.txt"
          output_GO_MF <- "./Mutation/output/Mutation_STADvsCRC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_Mutation(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaAF_cutoff)
        }
      }
      }
      
      ##Editing
      {
      #Differential editing analysis by editing ratio per site by wilcox rank sum test
      {
        mat_raw <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/Editing_ratio.txt",sep = "\t",header = TRUE, row.names = 1)
        des <- read.csv("./Editing/group/des_Editing_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Editing/matrix/Editing_site_for_GIvsNC.txt"
        output_res <- "./Editing/output/Editing_site_GIvsNC_wilcox.txt"
        Editing_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Editing/group/des_Editing_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Editing/matrix/Editing_site_for_CRCvsNC.txt"
        output_res <- "./Editing/output/Editing_site_CRCvsNC_wilcox.txt"
        Editing_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Editing/group/des_Editing_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Editing/matrix/Editing_site_for_STADvsNC.txt"
        output_res <- "./Editing/output/Editing_site_STADvsNC_wilcox.txt"
        Editing_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Editing/group/des_Editing_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Editing/matrix/Editing_site_for_STADvsCRC.txt"
        output_res <- "./Editing/output/Editing_site_STADvsCRC_wilcox.txt"
        Editing_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
      }
      
      #Functional
      {
        #Enrichment analysis: gastrointestinal cancer vs healthy donors
        {
          Differential_result <- read.csv("./Editing/output/Editing_site_GIvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Editing/output/Editing_GIvsNC_KEGG.txt"
          output_GO_BP <- "./Editing/output/Editing_GIvsNC_GO_BP.txt"
          output_GO_CC <- "./Editing/output/Editing_GIvsNC_GO_CC.txt"
          output_GO_MF <- "./Editing/output/Editing_GIvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_Editing(Differential_result,
                          output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                          pvalue_cutoff,deltaAF_cutoff)
        }
        #Enrichment analysis: CRC vs NC
        {
          Differential_result <- read.csv("./Editing/output/Editing_site_CRCvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Editing/output/Editing_CRCvsNC_KEGG.txt"
          output_GO_BP <- "./Editing/output/Editing_CRCvsNC_GO_BP.txt"
          output_GO_CC <- "./Editing/output/Editing_CRCvsNC_GO_CC.txt"
          output_GO_MF <- "./Editing/output/Editing_CRCvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_Editing(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaAF_cutoff)
        }
        #Enrichment analysis: STAD vs NC
        {
          Differential_result <- read.csv("./Editing/output/Editing_site_STADvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Editing/output/Editing_STADvsNC_KEGG.txt"
          output_GO_BP <- "./Editing/output/Editing_STADvsNC_GO_BP.txt"
          output_GO_CC <- "./Editing/output/Editing_STADvsNC_GO_CC.txt"
          output_GO_MF <- "./Editing/output/Editing_STADvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_Editing(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaAF_cutoff)
        }
        #Enrichment analysis: STAD vs CRC
        {
          Differential_result <- read.csv("./Editing/output/Editing_site_STADvsCRC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Editing/output/Editing_STADvsCRC_KEGG.txt"
          output_GO_BP <- "./Editing/output/Editing_STADvsCRC_GO_BP.txt"
          output_GO_CC <- "./Editing/output/Editing_STADvsCRC_GO_CC.txt"
          output_GO_MF <- "./Editing/output/Editing_STADvsCRC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_Editing(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaAF_cutoff)
        }
      }
      }
      
      ##ASE
      {
      #Differential editing analysis by editing ratio per site by wilcox rank sum test
      {
        mat_raw <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/ASE_COSMIC_multiomics_paired.txt",sep = "\t",header = TRUE, row.names = 1)
        des <- read.csv("./ASE/group/des_ASE_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./ASE/matrix/ASE_site_for_GIvsNC.txt"
        output_res <- "./ASE/output/ASE_site_GIvsNC_wilcox.txt"
        ASE_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./ASE/group/des_ASE_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./ASE/matrix/ASE_site_for_CRCvsNC.txt"
        output_res <- "./ASE/output/ASE_site_CRCvsNC_wilcox.txt"
        ASE_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./ASE/group/des_ASE_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./ASE/matrix/ASE_site_for_STADvsNC.txt"
        output_res <- "./ASE/output/ASE_site_STADvsNC_wilcox.txt"
        ASE_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./ASE/group/des_ASE_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./ASE/matrix/ASE_site_for_STADvsCRC.txt"
        output_res <- "./ASE/output/ASE_site_STADvsCRC_wilcox.txt"
        ASE_ratio_wilcox_test(mat_raw,des,output_matrix,output_res)
      }
      
      #Functional
      {
        #Enrichment analysis: gastrointestinal cancer vs healthy donors
        {
          Differential_result <- read.csv("./ASE/output/ASE_site_GIvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./ASE/output/ASE_GIvsNC_KEGG.txt"
          output_GO_BP <- "./ASE/output/ASE_GIvsNC_GO_BP.txt"
          output_GO_CC <- "./ASE/output/ASE_GIvsNC_GO_CC.txt"
          output_GO_MF <- "./ASE/output/ASE_GIvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.01
          
          KEGG_GO_ASE(Differential_result,
                      output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                      pvalue_cutoff,deltaAF_cutoff)
        }
        
        #Enrichment analysis: STAD vs NC
        {
          Differential_result <- read.csv("./ASE/output/ASE_site_STADvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./ASE/output/ASE_STADvsNC_KEGG.txt"
          output_GO_BP <- "./ASE/output/ASE_STADvsNC_GO_BP.txt"
          output_GO_CC <- "./ASE/output/ASE_STADvsNC_GO_CC.txt"
          output_GO_MF <- "./ASE/output/ASE_STADvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.01
          
          KEGG_GO_ASE(Differential_result,
                      output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                      pvalue_cutoff,deltaAF_cutoff)
        }
        
        #Enrichment analysis: CRC vs NC
        {
          Differential_result <- read.csv("./ASE/output/ASE_site_CRCvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./ASE/output/ASE_CRCvsNC_KEGG.txt"
          output_GO_BP <- "./ASE/output/ASE_CRCvsNC_GO_BP.txt"
          output_GO_CC <- "./ASE/output/ASE_CRCvsNC_GO_CC.txt"
          output_GO_MF <- "./ASE/output/ASE_CRCvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.01
          
          KEGG_GO_ASE(Differential_result,
                      output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                      pvalue_cutoff,deltaAF_cutoff)
        }
        
        #Enrichment analysis: STAD vs CRC
        {
          Differential_result <- read.csv("./ASE/output/ASE_site_STADvsCRC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./ASE/output/ASE_STADvsCRC_KEGG.txt"
          output_GO_BP <- "./ASE/output/ASE_STADvsCRC_GO_BP.txt"
          output_GO_CC <- "./ASE/output/ASE_STADvsCRC_GO_CC.txt"
          output_GO_MF <- "./ASE/output/ASE_STADvsCRC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.01
          
          KEGG_GO_ASE(Differential_result,
                      output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                      pvalue_cutoff,deltaAF_cutoff)
        }
      }
      }
      
      ##chimeric
      {
      #Differential analysis
      {
        mat_raw <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/JunctionReadCount.txt",sep = "\t",header = TRUE, row.names = 1,check.names = FALSE)
        des <- read.csv("./chimeric/group/des_chimeric_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./chimeric/matrix/chimeric_for_GIvsNC.txt"
        output_res <- "./chimeric/output/chimeric_GIvsNC_fisher_exact.txt"
        chimeric_fisher_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./chimeric/group/des_chimeric_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./chimeric/matrix/chimeric_for_CRCvsNC.txt"
        output_res <- "./chimeric/output/chimeric_CRCvsNC_fisher_exact.txt"
        chimeric_fisher_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./chimeric/group/des_chimeric_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./chimeric/matrix/chimeric_for_STADvsNC.txt"
        output_res <- "./chimeric/output/chimeric_STADvsNC_fisher_exact.txt"
        chimeric_fisher_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./chimeric/group/des_chimeric_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./chimeric/matrix/chimeric_for_STADvsCRC.txt"
        output_res <- "./chimeric/output/chimeric_STADvsCRC_fisher_exact.txt"
        chimeric_fisher_exact_test(mat_raw,des,output_matrix,output_res)
      }
      
      #functional
      {
        #Enrichment analysis: gastrointestinal cancer vs healthy donors
        {
          Differential_result <- read.csv("./chimeric/output/chimeric_GIvsNC_fisher_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./chimeric/output/chimeric_GIvsNC_KEGG.txt"
          output_GO_BP <- "./chimeric/output/chimeric_GIvsNC_GO_BP.txt"
          output_GO_CC <- "./chimeric/output/chimeric_GIvsNC_GO_CC.txt"
          output_GO_MF <- "./chimeric/output/chimeric_GIvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaFreq_cutoff <- 0
          
          KEGG_GO_chimeric(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaFreq_cutoff)
        }
        
        #Enrichment analysis: CRC vs NC
        {
          Differential_result <- read.csv("./chimeric/output/chimeric_CRCvsNC_fisher_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./chimeric/output/chimeric_CRCvsNC_KEGG.txt"
          output_GO_BP <- "./chimeric/output/chimeric_CRCvsNC_GO_BP.txt"
          output_GO_CC <- "./chimeric/output/chimeric_CRCvsNC_GO_CC.txt"
          output_GO_MF <- "./chimeric/output/chimeric_CRCvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaFreq_cutoff <- 0
          
          KEGG_GO_chimeric(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaFreq_cutoff)
        }
        
        #Enrichment analysis: STAD vs NC
        {
          Differential_result <- read.csv("./chimeric/output/chimeric_STADvsNC_fisher_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./chimeric/output/chimeric_STADvsNC_KEGG.txt"
          output_GO_BP <- "./chimeric/output/chimeric_STADvsNC_GO_BP.txt"
          output_GO_CC <- "./chimeric/output/chimeric_STADvsNC_GO_CC.txt"
          output_GO_MF <- "./chimeric/output/chimeric_STADvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaFreq_cutoff <- 0
          
          KEGG_GO_chimeric(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaFreq_cutoff)
        }
        
        #Enrichment analysis: STAD vs CRC
        {
          Differential_result <- read.csv("./chimeric/output/chimeric_STADvsCRC_fisher_exact.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./chimeric/output/chimeric_STADvsCRC_KEGG.txt"
          output_GO_BP <- "./chimeric/output/chimeric_STADvsCRC_GO_BP.txt"
          output_GO_CC <- "./chimeric/output/chimeric_STADvsCRC_GO_CC.txt"
          output_GO_MF <- "./chimeric/output/chimeric_STADvsCRC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaFreq_cutoff <- 0
          
          KEGG_GO_chimeric(Differential_result,
                           output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                           pvalue_cutoff,deltaFreq_cutoff)
        }
      }
      }
    
      ##Alternative promoter
      #Alt normalized output: /Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Alternative promoter/matrix/Altpromoter_normalized_forDE.txt
      {
      #Differential analysis
      {
        mat_raw <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/AlternativePromoter_multiomics_paired_normalized.txt",sep = "\t",header = TRUE, row.names = 1,check.names = FALSE)
        #write.table(mat_raw,"./Alternative promoter/matrix/Altpromoter_normalized_forDE2.txt",sep = "\t",quote = FALSE)
        des <- read.csv("./Alternative promoter/group/des_AlternativePromoter_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Alternative promoter/matrix/AlternativePromoter_for_GIvsNC.txt"
        output_res <- "./Alternative promoter/output/AlternativePromoter_GIvsNC_wilcox.txt"
        AlernativePromoter_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Alternative promoter/group/des_AlternativePromoter_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Alternative promoter/matrix/AlternativePromoter_for_CRCvsNC.txt"
        output_res <- "./Alternative promoter/output/AlternativePromoter_CRCvsNC_wilcox.txt"
        AlernativePromoter_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Alternative promoter/group/des_AlternativePromoter_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Alternative promoter/matrix/AlternativePromoter_for_STADvsNC.txt"
        output_res <- "./Alternative promoter/output/AlternativePromoter_STADvsNC_wilcox.txt"
        AlernativePromoter_wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Alternative promoter/group/des_AlternativePromoter_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Alternative promoter/matrix/AlternativePromoter_for_STADvsCRC.txt"
        output_res <- "./Alternative promoter/output/AlternativePromoter_STADvsCRC_wilcox.txt"
        AlernativePromoter_wilcox_test(mat_raw,des,output_matrix,output_res)
        
      }
      
      #Functional
      {
        #Enrichment analysis: gastrointestinal cancer vs healthy donors
        {
          Differential_result <- read.csv("./Alternative promoter/output/AlternativePromoter_GIvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Alternative promoter/output/AlternativePromoter_GIvsNC_KEGG.txt"
          output_GO_BP <- "./Alternative promoter/output/AlternativePromoter_GIvsNC_GO_BP.txt"
          output_GO_CC <- "./Alternative promoter/output/AlternativePromoter_GIvsNC_GO_CC.txt"
          output_GO_MF <- "./Alternative promoter/output/AlternativePromoter_GIvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_AlternativePromoter(Differential_result,
                                      output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                      pvalue_cutoff,deltaAF_cutoff)
        }
        
        #Enrichment analysis: CRC vs NC
        {
          Differential_result <- read.csv("./Alternative promoter/output/AlternativePromoter_CRCvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Alternative promoter/output/AlternativePromoter_CRCvsNC_KEGG.txt"
          output_GO_BP <- "./Alternative promoter/output/AlternativePromoter_CRCvsNC_GO_BP.txt"
          output_GO_CC <- "./Alternative promoter/output/AlternativePromoter_CRCvsNC_GO_CC.txt"
          output_GO_MF <- "./Alternative promoter/output/AlternativePromoter_CRCvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_AlternativePromoter(Differential_result,
                                      output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                      pvalue_cutoff,deltaAF_cutoff)
        }
        
        #Enrichment analysis: STAD vs NC
        {
          Differential_result <- read.csv("./Alternative promoter/output/AlternativePromoter_STADvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Alternative promoter/output/AlternativePromoter_STADvsNC_KEGG.txt"
          output_GO_BP <- "./Alternative promoter/output/AlternativePromoter_STADvsNC_GO_BP.txt"
          output_GO_CC <- "./Alternative promoter/output/AlternativePromoter_STADvsNC_GO_CC.txt"
          output_GO_MF <- "./Alternative promoter/output/AlternativePromoter_STADvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_AlternativePromoter(Differential_result,
                                      output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                      pvalue_cutoff,deltaAF_cutoff)
        }
        
        #Enrichment analysis: STAD vs CRC
        {
          Differential_result <- read.csv("./Alternative promoter/output/AlternativePromoter_STADvsCRC_wilcox.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Alternative promoter/output/AlternativePromoter_STADvsCRC_KEGG.txt"
          output_GO_BP <- "./Alternative promoter/output/AlternativePromoter_STADvsCRC_GO_BP.txt"
          output_GO_CC <- "./Alternative promoter/output/AlternativePromoter_STADvsCRC_GO_CC.txt"
          output_GO_MF <- "./Alternative promoter/output/AlternativePromoter_STADvsCRC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          deltaAF_cutoff <- 0.2
          
          KEGG_GO_AlternativePromoter(Differential_result,
                                      output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                                      pvalue_cutoff,deltaAF_cutoff)
        }
      }
      }
    
      ##Methylation
      {
      #Differential
      {
        mat_raw <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/CPM-TMM_matrix_promoter.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
        des <- read.csv("./Medip/group/des_Medip_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Medip/matrix/Medip_for_GIvsNC_promoter.txt"
        output_res <- "./Medip/output/Medip_GIvsNC_edger_exact_promoter.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Medip/group/des_Medip_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Medip/matrix/Medip_for_CRCvsNC_promoter.txt"
        output_res <- "./Medip/output/Medip_CRCvsNC_edger_exact_promoter.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Medip/group/des_Medip_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Medip/matrix/Medip_for_STADvsNC_promoter.txt"
        output_res <- "./Medip/output/Medip_STADvsNC_edger_exact_promoter.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Medip/group/des_Medip_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Medip/matrix/Medip_for_STADvsCRC_promoter.txt"
        output_res <- "./Medip/output/Medip_STADvsCRC_edger_exact_promoter.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
      }
      
      #Functional
      {
        #Enrichment analysis: gastrointestinal cancer vs healthy donors
        {
          Differential_result <- read.csv("./Medip/output/Medip_GIvsNC_edger_exact_promoter.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Medip/output/Methylation_GIvsNC_KEGG.txt"
          output_GO_BP <- "./Medip/output/Methylation_GIvsNC_GO_BP.txt"
          output_GO_CC <- "./Medip/output/Methylation_GIvsNC_GO_CC.txt"
          output_GO_MF <- "./Medip/output/Methylation_GIvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Methylation(Differential_result,
                              output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                              pvalue_cutoff,log2Foldchange_cutoff)
        }
        
        #Enrichment analysis: CRC vs healthy donors
        {
          Differential_result <- read.csv("./Medip/output/Medip_CRCvsNC_edger_exact_promoter.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Medip/output/Methylation_CRCvsNC_KEGG.txt"
          output_GO_BP <- "./Medip/output/Methylation_CRCvsNC_GO_BP.txt"
          output_GO_CC <- "./Medip/output/Methylation_CRCvsNC_GO_CC.txt"
          output_GO_MF <- "./Medip/output/Methylation_CRCvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Methylation(Differential_result,
                              output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                              pvalue_cutoff,log2Foldchange_cutoff)
        }
        
        #Enrichment analysis: STAD vs healthy donors
        {
          Differential_result <- read.csv("./Medip/output/Medip_STADvsNC_edger_exact_promoter.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Medip/output/Methylation_STADvsNC_KEGG.txt"
          output_GO_BP <- "./Medip/output/Methylation_STADvsNC_GO_BP.txt"
          output_GO_CC <- "./Medip/output/Methylation_STADvsNC_GO_CC.txt"
          output_GO_MF <- "./Medip/output/Methylation_STADvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Methylation(Differential_result,
                              output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                              pvalue_cutoff,log2Foldchange_cutoff)
        }
        
        #Enrichment analysis: STAD vs CRC
        {
          Differential_result <- read.csv("./Medip/output/Medip_STADvsCRC_edger_exact_promoter.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./Medip/output/Methylation_STADvsCRC_KEGG.txt"
          output_GO_BP <- "./Medip/output/Methylation_STADvsCRC_GO_BP.txt"
          output_GO_CC <- "./Medip/output/Methylation_STADvsCRC_GO_CC.txt"
          output_GO_MF <- "./Medip/output/Methylation_STADvsCRC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Methylation(Differential_result,
                              output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                              pvalue_cutoff,log2Foldchange_cutoff)
        }
      }
      }
    
      ##Methylation bin
      {
      #Differential
      {
        mat_raw <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/bin-all-cpgislands-shores-shelves-allsamples.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
        des <- read.csv("./Medip_bin/group/des_Medip_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Medip_bin/matrix/Medip_for_CRCvsNC_bin.txt"
        output_res <- "./Medip_bin/output/Medip_CRCvsNC_edger_exact_bin.txt"
        edgeR_exact_test_MethylationBin(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Medip_bin/group/des_Medip_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Medip_bin/matrix/Medip_for_STADvsNC_bin.txt"
        output_res <- "./Medip_bin/output/Medip_STADvsNC_edger_exact_bin.txt"
        edgeR_exact_test_MethylationBin(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Medip_bin/group/des_Medip_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Medip_bin/matrix/Medip_for_STADvsCRC_bin.txt"
        output_res <- "./Medip_bin/output/Medip_STADvsCRC_edger_exact_bin.txt"
        edgeR_exact_test_MethylationBin(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./Medip_bin/group/des_Medip_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./Medip_bin/matrix/Medip_for_GIvsNC_bin.txt"
        output_res <- "./Medip_bin/output/Medip_GIvsNC_edger_exact_bin.txt"
        edgeR_exact_test_MethylationBin(mat_raw,des,output_matrix,output_res)
      }
      }
    
      ##CNV
      {
      #Differential
      {
        mat_raw <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/CPM-TMM_matrix_gene.correctGC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
        des <- read.csv("./CNV/group/des_CNV_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./CNV/matrix/CNV_for_GIvsNC_gene.txt"
        output_res <- "./CNV/output/CNV_GIvsNC_edger_exact_gene.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./CNV/group/des_CNV_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./CNV/matrix/CNV_for_CRCvsNC_gene.txt"
        output_res <- "./CNV/output/CNV_CRCvsNC_edger_exact_gene.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./CNV/group/des_CNV_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./CNV/matrix/CNV_for_STADvsNC_gene.txt"
        output_res <- "./CNV/output/CNV_STADvsNC_edger_exact_gene.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./CNV/group/des_CNV_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./CNV/matrix/CNV_for_STADvsCRC_gene.txt"
        output_res <- "./CNV/output/CNV_STADvsCRC_edger_exact_gene.txt"
        edgeR_exact_test(mat_raw,des,output_matrix,output_res)
      }
      
      #Functional
      {
        #GIvsNC
        {
          Differential_result <- read.csv("./CNV/output/CNV_GIvsNC_edger_exact_gene.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./CNV/output/CNV_GIvsNC_KEGG.txt"
          output_GO_BP <- "./CNV/output/CNV_GIvsNC_GO_BP.txt"
          output_GO_CC <- "./CNV/output/CNV_GIvsNC_GO_CC.txt"
          output_GO_MF <- "./CNV/output/CNV_GIvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Nucleosome(Differential_result,
                             output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                             pvalue_cutoff,log2Foldchange_cutoff)
        }
        #CRCvsNC
        {
          Differential_result <- read.csv("./CNV/output/CNV_CRCvsNC_edger_exact_gene.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./CNV/output/CNV_CRCvsNC_KEGG.txt"
          output_GO_BP <- "./CNV/output/CNV_CRCvsNC_GO_BP.txt"
          output_GO_CC <- "./CNV/output/CNV_CRCvsNC_GO_CC.txt"
          output_GO_MF <- "./CNV/output/CNV_CRCvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Nucleosome(Differential_result,
                             output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                             pvalue_cutoff,log2Foldchange_cutoff)
        }
        #STADvsNC
        {
          Differential_result <- read.csv("./CNV/output/CNV_STADvsNC_edger_exact_gene.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./CNV/output/CNV_STADvsNC_KEGG.txt"
          output_GO_BP <- "./CNV/output/CNV_STADvsNC_GO_BP.txt"
          output_GO_CC <- "./CNV/output/CNV_STADvsNC_GO_CC.txt"
          output_GO_MF <- "./CNV/output/CNV_STADvsNC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Nucleosome(Differential_result,
                             output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                             pvalue_cutoff,log2Foldchange_cutoff)
        }
        #STADvsCRC
        {
          Differential_result <- read.csv("./CNV/output/CNV_STADvsCRC_edger_exact_gene.txt",header = T, row.names = 1,sep = "\t")
          output_KEGG <- "./CNV/output/CNV_STADvsCRC_KEGG.txt"
          output_GO_BP <- "./CNV/output/CNV_STADvsCRC_GO_BP.txt"
          output_GO_CC <- "./CNV/output/CNV_STADvsCRC_GO_CC.txt"
          output_GO_MF <- "./CNV/output/CNV_STADvsCRC_GO_MF.txt"
          pvalue_cutoff <- 0.05
          log2Foldchange_cutoff <- 0.59
          
          KEGG_GO_Nucleosome(Differential_result,
                             output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                             pvalue_cutoff,log2Foldchange_cutoff)
        }
      }
      }
      
      ##FragSize
      {
      #Differential
      {
        mat_raw <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/DNA-FragRatio_matrix_bin100kb.correctGC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
        des <- read.csv("./FragSize/group/des_FragSize_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./FragSize/matrix/FragSize_for_GIvsNC_bin100kb.correctGC.txt"
        output_res <- "./FragSize/output/FragSize_GIvsNC_wilcox_bin100kb.correctGC.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./FragSize/group/des_FragSize_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./FragSize/matrix/FragSize_for_CRCvsNC_bin100kb.correctGC.txt"
        output_res <- "./FragSize/output/FragSize_CRCvsNC_wilcox_bin100kb.correctGC.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./FragSize/group/des_FragSize_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./FragSize/matrix/FragSize_for_STADvsNC_bin100kb.correctGC.txt"
        output_res <- "./FragSize/output/FragSize_STADvsNC_wilcox_bin100kb.correctGC.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./FragSize/group/des_FragSize_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./FragSize/matrix/FragSize_for_STADvsCRC_bin100kb.correctGC.txt"
        output_res <- "./FragSize/output/FragSize_STADvsCRC_wilcox_bin100kb.correctGC.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
      }
      }
      
      ##EndMotif
      {
      #Differential
      {
        mat_raw <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/DNA-EndMotifRatio_matrix_4mers.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
        des <- read.csv("./EndMotif/group/des_EndMotif_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./EndMotif/matrix/EndMotif_for_GIvsNC_4mers.txt"
        output_res <- "./EndMotif/output/EndMotif_GIvsNC_wilcox_4mers.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./EndMotif/group/des_EndMotif_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./EndMotif/matrix/EndMotif_for_CRCvsNC_4mers.txt"
        output_res <- "./EndMotif/output/EndMotif_CRCvsNC_wilcox_4mers.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./EndMotif/group/des_EndMotif_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./EndMotif/matrix/EndMotif_for_STADvsNC_4mers.txt"
        output_res <- "./EndMotif/output/EndMotif_STADvsNC_wilcox_4mers.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./EndMotif/group/des_EndMotif_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./EndMotif/matrix/EndMotif_for_STADvsCRC_4mers.txt"
        output_res <- "./EndMotif/output/EndMotif_STADvsCRC_wilcox_4mers.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
      }
      }
      
      ##WPS liyu
      {
      #Differential
      {
        mat_raw <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/WPS-divide-bgCOV_matrix_150TSS50.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
        des <- read.csv("./WPS/group/des_WPS_CRCvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./WPS/matrix/WPS_for_CRCvsNC.txt"
        output_res <- "./WPS/output/WPS_CRCvsNC_wilcox.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./WPS/group/des_WPS_STADvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./WPS/matrix/WPS_for_STADvsNC.txt"
        output_res <- "./WPS/output/WPS_STADvsNC_wilcox.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./WPS/group/des_WPS_STADvsCRC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./WPS/matrix/WPS_for_STADvsCRC.txt"
        output_res <- "./WPS/output/WPS_STADvsCRC_wilcox.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
        
        des <- read.csv("./WPS/group/des_WPS_GIvsNC.csv", header = TRUE, check.names=FALSE, sep=',')
        output_matrix <- "./WPS/matrix/WPS_for_GIvsNC.txt"
        output_res <- "./WPS/output/WPS_GIvsNC_wilcox.txt"
        Wilcox_test(mat_raw,des,output_matrix,output_res)
      }
        
      #Functional
      {
          #GIvsNC
          {
            Differential_result <- read.csv("./WPS/output/WPS_GIvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
            output_KEGG <- "./WPS/output/WPS_GIvsNC_KEGG.txt"
            output_GO_BP <- "./WPS/output/WPS_GIvsNC_GO_BP.txt"
            output_GO_CC <- "./WPS/output/WPS_GIvsNC_GO_CC.txt"
            output_GO_MF <- "./WPS/output/WPS_GIvsNC_GO_MF.txt"
            pvalue_cutoff <- 0.05
            deltaAF_cutoff <- 0.5
            
            KEGG_GO_WPS(Differential_result,
                               output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                               pvalue_cutoff,deltaAF_cutoff)
          }
          #CRCvsNC
          {
            Differential_result <- read.csv("./WPS/output/WPS_CRCvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
            output_KEGG <- "./WPS/output/WPS_CRCvsNC_KEGG.txt"
            output_GO_BP <- "./WPS/output/WPS_CRCvsNC_GO_BP.txt"
            output_GO_CC <- "./WPS/output/WPS_CRCvsNC_GO_CC.txt"
            output_GO_MF <- "./WPS/output/WPS_CRCvsNC_GO_MF.txt"
            pvalue_cutoff <- 0.05
            deltaAF_cutoff <- 0.5
            
            KEGG_GO_WPS(Differential_result,
                               output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                               pvalue_cutoff,deltaAF_cutoff)
          }
          #STADvsNC
          {
            Differential_result <- read.csv("./WPS/output/WPS_STADvsNC_wilcox.txt",header = T, row.names = 1,sep = "\t")
            output_KEGG <- "./WPS/output/WPS_STADvsNC_KEGG.txt"
            output_GO_BP <- "./WPS/output/WPS_STADvsNC_GO_BP.txt"
            output_GO_CC <- "./WPS/output/WPS_STADvsNC_GO_CC.txt"
            output_GO_MF <- "./WPS/output/WPS_STADvsNC_GO_MF.txt"
            pvalue_cutoff <- 0.05
            deltaAF_cutoff <- 0.5
            
            KEGG_GO_WPS(Differential_result,
                               output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                               pvalue_cutoff,deltaAF_cutoff)
          }
          #STADvsCRC
          {
            Differential_result <- read.csv("./WPS/output/WPS_STADvsCRC_wilcox.txt",header = T, row.names = 1,sep = "\t")
            output_KEGG <- "./WPS/output/WPS_STADvsCRC_KEGG.txt"
            output_GO_BP <- "./WPS/output/WPS_STADvsCRC_GO_BP.txt"
            output_GO_CC <- "./WPS/output/WPS_STADvsCRC_GO_CC.txt"
            output_GO_MF <- "./WPS/output/WPS_STADvsCRC_GO_MF.txt"
            pvalue_cutoff <- 0.05
            deltaAF_cutoff <- 0.5
            
            KEGG_GO_WPS(Differential_result,
                               output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC,
                               pvalue_cutoff,deltaAF_cutoff)
          }
        }
      }
    }
  }
  
  #Differential event number
  {
    ML_samples <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/sample/ML_samples.csv",header = TRUE, row.names = 1, check.names = FALSE)
    ##Expression
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./Expression/output/Expression_STADvsCRC_edger_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      output_res <- read.table("./Expression/output/Expression_STADvsNC_edger_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      output_res <- read.table("./Expression/output/Expression_CRCvsNC_edger_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./Expression/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##miRNA
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      
      output_res <- read.table("./miRNA/output/miRNA_STADvsCRC_edger_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      output_res <- read.table("./miRNA/output/miRNA_STADvsNC_edger_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      output_res <- read.table("./miRNA/output/miRNA_CRCvsNC_edger_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./miRNA/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
      
    }
    
    ##Alternative promoter
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./Alternative promoter/output/AlternativePromoter_STADvsCRC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.2,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.2,])
      
      output_res <- read.table("./Alternative promoter/output/AlternativePromoter_STADvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.2,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.2,])
      
      output_res <- read.table("./Alternative promoter/output/AlternativePromoter_CRCvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.2,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.2,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./Alternative promoter/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##ASE
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./ASE/output/ASE_site_STADvsCRC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.1,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.1,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.1,])
      
      output_res <- read.table("./ASE/output/ASE_site_STADvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.1,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.1,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.1,])
      
      output_res <- read.table("./ASE/output/ASE_site_CRCvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.1,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.1,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.1,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./ASE/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##chimeric
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./chimeric/output/chimeric_STADvsCRC_fisher_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaFreq*output_res$deltaFreq) > 0.1,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaFreq > 0.1,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaFreq < -0.1,])
      
      output_res <- read.table("./chimeric/output/chimeric_STADvsNC_fisher_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaFreq*output_res$deltaFreq) > 0.1,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaFreq > 0.1,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaFreq < -0.1,])
      
      output_res <- read.table("./chimeric/output/chimeric_CRCvsNC_fisher_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaFreq*output_res$deltaFreq) > 0.1,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaFreq > 0.1,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaFreq < -0.1,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./chimeric/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##Editing
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./Editing/output/Editing_site_STADvsCRC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.2,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.2,])
      
      output_res <- read.table("./Editing/output/Editing_site_STADvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.2,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.2,])
      
      output_res <- read.table("./Editing/output/Editing_site_CRCvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.2,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.2,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./Editing/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##Mutation
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./Mutation/output/Mutation_site_STADvsCRC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.2,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.2,])
      
      output_res <- read.table("./Mutation/output/Mutation_site_STADvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.2,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.2,])
      
      output_res <- read.table("./Mutation/output/Mutation_site_CRCvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.2,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.2,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./Mutation/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##Splicing
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./Splicing/output/multiomics_paired_STAD_vs_CRC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$PValue<0.05 & sqrt(output_res$IncLevelDifference*output_res$IncLevelDifference) > 0.05,])
      STADvsCRC_features_up <- rownames(output_res[output_res$PValue<0.05 & output_res$IncLevelDifference > 0.05,])
      STADvsCRC_features_down <- rownames(output_res[output_res$PValue<0.05 & output_res$IncLevelDifference < -0.05,])
      
      output_res <- read.table("./Splicing/output/multiomics_paired_STAD_vs_Healthy.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$PValue<0.05 & sqrt(output_res$IncLevelDifference*output_res$IncLevelDifference) > 0.05,])
      STADvsNC_features_up <- rownames(output_res[output_res$PValue<0.05 & output_res$IncLevelDifference > 0.05,])
      STADvsNC_features_down <- rownames(output_res[output_res$PValue<0.05 & output_res$IncLevelDifference < -0.05,])
  
      output_res <- read.table("./Splicing/output/multiomics_paired_CRC_vs_Healthy.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$PValue<0.05 & sqrt(output_res$IncLevelDifference*output_res$IncLevelDifference) > 0.05,])
      CRCvsNC_features_up <- rownames(output_res[output_res$PValue<0.05 & output_res$IncLevelDifference > 0.05,])
      CRCvsNC_features_down <- rownames(output_res[output_res$PValue<0.05 & output_res$IncLevelDifference < -0.05,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./Splicing/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##Medip
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./Medip/output/Medip_STADvsCRC_edger_exact_promoter.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      output_res <- read.table("./Medip/output/Medip_STADvsNC_edger_exact_promoter.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      output_res <- read.table("./Medip/output/Medip_CRCvsNC_edger_exact_promoter.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./Medip/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##Medip_bin
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./Medip_bin/output/Medip_STADvsCRC_edger_exact_bin.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      output_res <- read.table("./Medip_bin/output/Medip_STADvsNC_edger_exact_bin.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      output_res <- read.table("./Medip_bin/output/Medip_CRCvsNC_edger_exact_bin.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./Medip_bin/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##CNV
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./CNV/output/CNV_STADvsCRC_edger_exact_gene.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      output_res <- read.table("./CNV/output/CNV_STADvsNC_edger_exact_gene.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      output_res <- read.table("./CNV/output/CNV_CRCvsNC_edger_exact_gene.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange > 0.59,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$log2FoldChange < -0.59,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./CNV/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##EndMotif
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./EndMotif/output/EndMotif_STADvsCRC_wilcox_4mers.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < 0,])
      
      output_res <- read.table("./EndMotif/output/EndMotif_STADvsNC_wilcox_4mers.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < 0,])
      
      output_res <- read.table("./EndMotif/output/EndMotif_CRCvsNC_wilcox_4mers.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < 0,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./EndMotif/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##FragSize
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./FragSize/output/FragSize_STADvsCRC_wilcox_bin100kb.correctGC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < 0,])
      
      output_res <- read.table("./FragSize/output/FragSize_STADvsNC_wilcox_bin100kb.correctGC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < 0,])
      
      output_res <- read.table("./FragSize/output/FragSize_CRCvsNC_wilcox_bin100kb.correctGC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < 0,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./FragSize/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    ##WPS
    {
      differential_event_number <- as.data.frame(array(,dim=c(3,3)))
      colnames(differential_event_number) <- c("Total","Up-regulated","Down-regulated")
      rownames(differential_event_number) <- c("STADvsCRC","STADvsHD","CRCvsHD")
      
      output_res <- read.table("./WPS/output/WPS_STADvsCRC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsCRC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.5,])
      STADvsCRC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.5,])
      STADvsCRC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.5,])
      
      output_res <- read.table("./WPS/output/WPS_STADvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.5,])
      STADvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.5,])
      STADvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.5,])
      
      output_res <- read.table("./WPS/output/WPS_CRCvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.5,])
      CRCvsNC_features_up <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF > 0.5,])
      CRCvsNC_features_down <- rownames(output_res[output_res$pvalue<0.05 & output_res$deltaAF < -0.5,])
      
      differential_event_number["STADvsCRC","Total"] <- length(STADvsCRC_features)
      differential_event_number["STADvsCRC","Up-regulated"] <- length(STADvsCRC_features_up)
      differential_event_number["STADvsCRC","Down-regulated"] <- length(STADvsCRC_features_down)
      differential_event_number["STADvsHD","Total"] <- length(STADvsNC_features)
      differential_event_number["STADvsHD","Up-regulated"] <- length(STADvsNC_features_up)
      differential_event_number["STADvsHD","Down-regulated"] <- length(STADvsNC_features_down)
      differential_event_number["CRCvsHD","Total"] <- length(CRCvsNC_features)
      differential_event_number["CRCvsHD","Up-regulated"] <- length(CRCvsNC_features_up)
      differential_event_number["CRCvsHD","Down-regulated"] <- length(CRCvsNC_features_down)
      
      write.table(differential_event_number,"./WPS/output/Differential_event_number.txt",sep = "\t",quote=FALSE)
    }
    
    Alterations <- c("Expression","Alternative promoter","ASE","chimeric",
                     "Editing","Mutation","Splicing","miRNA",
                     "Medip","Medip_bin",
                     "CNV","EndMotif","FragSize","WPS")
    
    i=1
    Statistics_all <- {}
    #Statistics_all <- as.data.frame(array(,dim=c(0,4)))
    #colnames(Statistics_all) <- c("Total","Up-regulated","Down-regulated","Alteration")
    while(i<=length(Alterations)){
      Statistics_i <- read.table(paste0("./",Alterations[i],"/output/Differential_event_number.txt"),sep = "\t", header = TRUE,row.names = NULL)
      Statistics_i$Alteration <- Alterations[i]
      Statistics_all <- rbind(Statistics_all,Statistics_i)
      
      i=i+1
    }
    
    
    #RNA_color <- c("Human gene abundance"="#A5435C","cfRNA alternative promoter usage"="#A5435C","cfRNA allele specific expression"="#A5435C","Chimeric cfRNA"="#A5435C",
    #               "cfRNA editing"="#A5435C","Microbe genus abundance"="#A5435C","cfRNA SNV"="#A5435C","cfRNA splicing"="#A5435C","cfRNA transposable elements abudance"="#A5435C",
    #               "cf-miRNA abundance"=alpha("#A5435C",alpha=0.5),
    #               "cfDNA methylation"="#44A5F9",
    #               "cfDNA copy number"="#FFE234","cfDNA nucleosome occpuancy"="#FFE234","cfDNA end motif usage"="#FFE234","cfDNA fragment size"="#FFE234"
    #               )
    
    RNA_color <- c("cfRNA abundance"="#5CACEE","cfRNA alternative promoter usage"="#5CACEE","cfRNA allele specific expression"="#5CACEE","Chimeric cfRNA"="#5CACEE",
                   "cfRNA editing"="#5CACEE","cfRNA abundance (microbe genus abundance)"="#5CACEE","cfRNA SNV"="#5CACEE","cfRNA splicing"="#5CACEE","cfRNA abudance (transposable elements)"="#5CACEE",
                   "cf-miRNA abundance"=alpha("#0a14a6",alpha=1),
                   "cfDNA methylation (promoter)"="#666666","cfDNA methylation (CpG island)"="#666666",
                   "cfDNA copy number"="#EEB4B4","cfDNA window protection score"="#EEB4B4","cfDNA end motif usage"="#EEB4B4","cfDNA fragment size"="#EEB4B4"
    )
    #scale_fill_manual(values=c(alpha("#5CACEE",alpha = 1),alpha("#EEB4B4",alpha = 0.45),alpha("#666666",alpha = 0.1),alpha("#CD0000",alpha = 0.8)))
    
    {
    comparison <- "CRCvsHD"
    Differential_events_plot <- Statistics_all[which(Statistics_all$row.names==comparison),]
    Differential_events_plot$Alteration <- gsub("Expression","cfRNA abundance",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("Alternative promoter","cfRNA alternative promoter usage",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("ASE","cfRNA allele specific expression",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("chimeric","Chimeric cfRNA",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("Editing","cfRNA editing",fixed=TRUE,Differential_events_plot$Alteration)
    #Differential_events_plot$Alteration <- gsub("Microbe","cfRNA abundance (microbe genus abundance)",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("Mutation","cfRNA SNV",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("Splicing","cfRNA splicing",fixed=TRUE,Differential_events_plot$Alteration)
    #Differential_events_plot$Alteration <- gsub("TE","cfRNA abudance (transposable elements)",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("miRNA","cf-miRNA abundance",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("Medip_bin","cfDNA methylation (CpG island)",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("Medip","cfDNA methylation (promoter)",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("CNV","cfDNA copy number",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("WPS","cfDNA window protection score",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("EndMotif","cfDNA end motif usage",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- gsub("FragSize","cfDNA fragment size",fixed=TRUE,Differential_events_plot$Alteration)
    Differential_events_plot$Alteration <- factor(Differential_events_plot$Alteration,
                                                  levels = rev(rev(as.character(Differential_events_plot[order(Differential_events_plot$Total,decreasing = TRUE),]$Alteration))))
    }
    
    p <- ggplot(Differential_events_plot)+
      geom_bar(aes(x=Alteration,y=-log10(Down.regulated+1),fill=Alteration),stat = "identity",colour = "black")+
      geom_text(aes(x=Alteration,y=-log10(20000),label=Down.regulated),size = 4,angle = 0, hjust = 0)+
      geom_bar(aes(x=Alteration,y=log10(Up.regulated+1),fill=Alteration),stat = "identity",colour = "black")+
      geom_text(aes(x=Alteration,y=log10(20000),label=Up.regulated),size = 4,angle = 0, hjust = 1)+
      scale_fill_manual(values= alpha(RNA_color,alpha = 1))+
      geom_hline(aes(yintercept = 0))+
      coord_flip()+
      xlab("")+
      ylab("")+
      labs(title=comparison)+
      scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),labels = c("1000","100","10","0","10","100","1000"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(15,5,10,15),units="pt"),
        legend.position='none',
        panel.grid=element_blank(),
        panel.border=element_rect(color = "black",size = 0),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = -0.6,size=20,face="bold"),
        #axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_blank(),
        #axis.line.y = element_line(color = "black"),
        axis.line = element_blank(),
        axis.text.x = element_text(color="black", size=16, angle = 45,hjust = 1,vjust = 1),
        #axis.text.y = element_blank(),
        axis.text.y = element_text( color="black", size=16, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))
    ggsave(plot = p,filename = paste0("./",comparison,"_event_number.pdf"),width=8.52,height=5.47,device = "pdf")
  }
  
  #gini index density and PCA
  {
    Group_color <- c(CRC=alpha("#FCB514",alpha = 0.5),STAD=alpha("red",alpha = 0.5),HD=alpha("blue",alpha = 0.5))
    comparison_type <- "CRC_STAD_vs_HD"
    ML_samples <- read.csv("./Alterations/sample/ML_samples.csv",header = TRUE, row.names = 1, check.names = FALSE)
    #Group_color <- c(CRC=alpha("#FCB514",alpha = 0.5),STAD=alpha("red",alpha = 0.5),HD=alpha("blue",alpha = 0.5),GIC="#FCB514")
    #Splicing
    {
      dir.create("./Figure 3/Splicing/gini")
      dir.create("./Figure 3/Splicing/PCA")
      Splicing <- read.table("./Alterations/matrix_raw/RNA/multiomics_paired_Cancer_vs_Healthy_Inclevel_matrix.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- Splicing
      Splicing <- na.omit(Splicing)
      #Splicing[is.na(Splicing)] <- 1 
      
      output_res <- read.table("./Figure 3/Splicing/output/multiomics_paired_STAD_vs_Healthy.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$PValue<0.05 & sqrt(output_res$IncLevelDifference*output_res$IncLevelDifference) > 0.05,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/Splicing/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/Splicing/output/multiomics_paired_CRC_vs_Healthy.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$PValue<0.05 & sqrt(output_res$IncLevelDifference*output_res$IncLevelDifference) > 0.05,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/Splicing/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
      Splicing <- Splicing[top_features,]
      
      STAD_gini <- data.frame(gini=gini(t(Splicing[,grep("STAD",colnames(Splicing))])),group="STAD")
      CRC_gini <- data.frame(gini=gini(t(Splicing[,grep("CRC",colnames(Splicing))])),group="CRC")
      NC_gini <- data.frame(gini=gini(t(Splicing[,grep("NC",colnames(Splicing))])),group="HD")
      gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
      
      STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
      CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
      STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
      
      ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                            "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                            "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                            "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
      
      write.csv(ks_test,"./Figure 3/Splicing/gini/ks_test.csv")
      
      p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
        geom_density(color="black")+
        #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
        scale_fill_manual(values=Group_color)+
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
        scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
        xlab("gini index")+
        ylab("Density")+
        theme_bw()+
        theme(
          #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))
      ggsave("./Figure 3/Splicing/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(ML_samples$RNA_id)])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/Splicing/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/Splicing/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #Expression
    {
      dir.create("./Figure 3/Expression/gini")
      dir.create("./Figure 3/Expression/PCA")
      Expression <- read.table("./Alterations/matrix_raw/RNA/Plasma_TPM.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- Expression
      Expression <- na.omit(Expression)
      #Expression[is.na(Expression)] <- 1 
      
      output_res <- read.table("./Figure 3/Expression/output/Expression_STADvsNC_edger_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/Expression/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/Expression/output/Expression_CRCvsNC_edger_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/Expression/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        Expression <- Expression[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(Expression[,grep("STAD",colnames(Expression))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(Expression[,grep("CRC",colnames(Expression))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(Expression[,grep("NC",colnames(Expression))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/Expression/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/Expression/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(ML_samples$RNA_id)])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/Expression/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/Expression/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #miRNA
    {
      dir.create("./Figure 3/miRNA/gini")
      dir.create("./Figure 3/miRNA/PCA")
      miRNA <- read.table("./Alterations/matrix_raw/RNA/miRNA.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- miRNA
      miRNA <- na.omit(miRNA)
      #miRNA[is.na(miRNA)] <- 1 
      
      output_res <- read.table("./Figure 3/miRNA/output/miRNA_STADvsNC_edger_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/miRNA/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/miRNA/output/miRNA_CRCvsNC_edger_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/miRNA/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        miRNA <- miRNA[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(miRNA[,grep("STAD",colnames(miRNA))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(miRNA[,grep("CRC",colnames(miRNA))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(miRNA[,grep("NC",colnames(miRNA))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/miRNA/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/miRNA/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$miRNA_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/miRNA/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/miRNA/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #Editing
    {
      dir.create("./Figure 3/Editing/gini")
      dir.create("./Figure 3/Editing/PCA")
      Editing <- read.table("./Alterations/matrix_raw/RNA/Editing_ratio.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- Editing
      Editing <- na.omit(Editing)
      #Editing[is.na(Editing)] <- 1 
      
      output_res <- read.table("./Figure 3/Editing/output/Editing_site_STADvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/Editing/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/Editing/output/Editing_site_CRCvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/Editing/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        Editing <- Editing[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(Editing[,grep("STAD",colnames(Editing))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(Editing[,grep("CRC",colnames(Editing))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(Editing[,grep("NC",colnames(Editing))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/Editing/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/Editing/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$RNA_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/Editing/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/Editing/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #ASE
    {
      dir.create("./Figure 3/ASE/gini")
      dir.create("./Figure 3/ASE/PCA")
      ASE <- read.table("./Alterations/matrix_raw/RNA/ASE_COSMIC_multiomics_paired.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- ASE
      ASE <- na.omit(ASE)
      #ASE[is.na(ASE)] <- 1 
      
      output_res <- read.table("./Figure 3/ASE/output/ASE_site_STADvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.1,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/ASE/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/ASE/output/ASE_site_CRCvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.1,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/ASE/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        ASE <- ASE[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(ASE[,grep("STAD",colnames(ASE))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(ASE[,grep("CRC",colnames(ASE))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(ASE[,grep("NC",colnames(ASE))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/ASE/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/ASE/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$RNA_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/ASE/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/ASE/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #Mutation
    {
      dir.create("./Figure 3/Mutation/gini")
      dir.create("./Figure 3/Mutation/PCA")
      Mutation <- read.table("./Alterations/matrix_raw/RNA/final_SNP_AF.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- Mutation
      Mutation <- na.omit(Mutation)
      #Mutation[is.na(Mutation)] <- 1 
      
      output_res <- read.table("./Figure 3/Mutation/output/Mutation_site_STADvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/Mutation/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/Mutation/output/Mutation_site_CRCvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/Mutation/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        Mutation <- Mutation[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(Mutation[,grep("STAD",colnames(Mutation))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(Mutation[,grep("CRC",colnames(Mutation))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(Mutation[,grep("NC",colnames(Mutation))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/Mutation/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/Mutation/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$RNA_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/Mutation/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/Mutation/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #Altpromoter
    {
      dir.create("./Figure 3/Alternative promoter/gini")
      dir.create("./Figure 3/Alternative promoter/PCA")
      AlternativePromoter <- read.table("./Alterations/matrix_raw/RNA/AlternativePromoter_multiomics_paired_normalized.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- AlternativePromoter
      AlternativePromoter <- na.omit(AlternativePromoter)
      #AlternativePromoter[is.na(AlternativePromoter)] <- 1 
      
      output_res <- read.table("./Figure 3/Alternative promoter/output/AlternativePromoter_STADvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/Alternative promoter/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/Alternative promoter/output/AlternativePromoter_CRCvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.2,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/Alternative promoter/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        AlternativePromoter <- AlternativePromoter[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(AlternativePromoter[,grep("STAD",colnames(AlternativePromoter))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(AlternativePromoter[,grep("CRC",colnames(AlternativePromoter))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(AlternativePromoter[,grep("NC",colnames(AlternativePromoter))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/Alternative promoter/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/Alternative promoter/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$RNA_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/Alternative promoter/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/Alternative promoter/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #chimeric (do not have significant differential chimeric alterations, so use all events for gini and PCA)
    {
      dir.create("./Figure 3/chimeric/gini")
      dir.create("./Figure 3/chimeric/PCA")
      chimeric <- read.table("./Alterations/matrix_raw/RNA/JunctionReadCount.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- chimeric
      chimeric <- na.omit(chimeric)
      #chimeric[is.na(chimeric)] <- 1 
      
      output_res <- read.table("./Figure 3/chimeric/output/chimeric_STADvsNC_fisher_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaFreq*output_res$deltaFreq) > 0.1,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/chimeric/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/chimeric/output/chimeric_CRCvsNC_fisher_exact.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaFreq*output_res$deltaFreq) > 0.1,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/chimeric/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        chimeric <- chimeric#[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(chimeric[,grep("STAD",colnames(chimeric))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(chimeric[,grep("CRC",colnames(chimeric))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(chimeric[,grep("NC",colnames(chimeric))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STADgini <- STAD_gini$gini
        STADgini <- STADgini[!is.nan(STADgini)]
        CRCgini <- CRC_gini$gini
        CRCgini <- CRCgini[!is.nan(CRCgini)]
        NCgini <- NC_gini$gini
        NCgini <- NCgini[!is.nan(NCgini)]
        
        STAD_NC <- ks.test(STADgini,NCgini)
        CRC_NC <- ks.test(CRCgini,NCgini)
        STAD_CRC <- ks.test(STADgini,CRCgini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/chimeric/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/chimeric/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[,as.character(na.omit(ML_samples$RNA_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/chimeric/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/chimeric/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #CNV
    {
      dir.create("./Figure 3/CNV/gini")
      dir.create("./Figure 3/CNV/PCA")
      CNV <- read.table("./Alterations/matrix_raw/DNA/CPM-TMM_matrix_gene.correctGC.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- CNV
      CNV <- na.omit(CNV)
      #CNV[is.na(CNV)] <- 1 
      
      output_res <- read.table("./Figure 3/CNV/output/CNV_STADvsNC_edger_exact_gene.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/CNV/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/CNV/output/CNV_CRCvsNC_edger_exact_gene.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/CNV/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        CNV <- CNV[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(CNV[,grep("STAD",colnames(CNV))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(CNV[,grep("CRC",colnames(CNV))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(CNV[,grep("NC",colnames(CNV))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/CNV/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/CNV/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$DNA_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/CNV/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/CNV/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #WPS
    {
      dir.create("./Figure 3/WPS/gini")
      dir.create("./Figure 3/WPS/PCA")
      WPS <- read.table("./Alterations/matrix_raw/DNA/WPS-divide-bgCOV_matrix_150TSS50.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- WPS
      WPS <- na.omit(WPS)
      #WPS[is.na(WPS)] <- 1 
      
      output_res <- read.table("./Figure 3/WPS/output/WPS_STADvsCRC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.5,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/WPS/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/WPS/output/WPS_CRCvsNC_wilcox.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0.5,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/WPS/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        WPS <- WPS[top_features,]
        
        STAD_gini <- data.frame(gini=gini(-t(WPS[,grep("STAD",colnames(WPS))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(-t(WPS[,grep("CRC",colnames(WPS))])),group="CRC")
        NC_gini <- data.frame(gini=gini(-t(WPS[,grep("NC",colnames(WPS))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/WPS/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/WPS/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$DNA_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/WPS/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/WPS/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #EndMotif
    {
      dir.create("./Figure 3/EndMotif/gini")
      dir.create("./Figure 3/EndMotif/PCA")
      EndMotif <- read.table("./Alterations/matrix_raw/DNA/DNA-EndMotifRatio_matrix_4mers.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- EndMotif
      EndMotif <- na.omit(EndMotif)
      #EndMotif[is.na(EndMotif)] <- 1 
      
      output_res <- read.table("./Figure 3/EndMotif/output/EndMotif_STADvsNC_wilcox_4mers.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/EndMotif/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/EndMotif/output/EndMotif_CRCvsNC_wilcox_4mers.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/EndMotif/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        EndMotif <- EndMotif[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(EndMotif[,grep("STAD",colnames(EndMotif))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(EndMotif[,grep("CRC",colnames(EndMotif))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(EndMotif[,grep("NC",colnames(EndMotif))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/EndMotif/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/EndMotif/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$DNA_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/EndMotif/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/EndMotif/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #FragSize
    {
      dir.create("./Figure 3/FragSize/gini")
      dir.create("./Figure 3/FragSize/PCA")
      FragSize <- read.table("./Alterations/matrix_raw/DNA/DNA-FragRatio_matrix_bin100kb.correctGC.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- FragSize
      FragSize <- na.omit(FragSize)
      #FragSize[is.na(FragSize)] <- 1 
      
      output_res <- read.table("./Figure 3/FragSize/output/FragSize_STADvsNC_wilcox_bin100kb.correctGC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/FragSize/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/FragSize/output/FragSize_CRCvsNC_wilcox_bin100kb.correctGC.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$deltaAF*output_res$deltaAF) > 0,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/FragSize/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        FragSize <- FragSize[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(FragSize[,grep("STAD",colnames(FragSize))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(FragSize[,grep("CRC",colnames(FragSize))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(FragSize[,grep("NC",colnames(FragSize))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/FragSize/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/FragSize/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$DNA_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/FragSize/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/FragSize/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #Medip(promoter)
    {
      dir.create("./Figure 3/Medip/gini")
      dir.create("./Figure 3/Medip/PCA")
      Medip <- read.table("./Alterations/matrix_raw/DNA/CPM-TMM_matrix_promoter.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- Medip
      Medip <- na.omit(Medip)
      #Medip[is.na(Medip)] <- 1 
      
      output_res <- read.table("./Figure 3/Medip/output/Medip_STADvsNC_edger_exact_promoter.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/Medip/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/Medip/output/Medip_CRCvsNC_edger_exact_promoter.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/Medip/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        Medip <- Medip[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(Medip[,grep("STAD",colnames(Medip))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(Medip[,grep("CRC",colnames(Medip))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(Medip[,grep("NC",colnames(Medip))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/Medip/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/Medip/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$Methylation_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/Medip/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/Medip/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
    
    #Medip(CpG island)
    {
      dir.create("./Figure 3/Medip_bin/gini")
      dir.create("./Figure 3/Medip_bin/PCA")
      Medip_bin <- read.table("./Alterations/matrix_raw/DNA/bin-all-cpgislands-shores-shelves-allsamples.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      full_matrix <- Medip_bin
      Medip_bin <- na.omit(Medip_bin)
      #Medip_bin[is.na(Medip_bin)] <- 1 
      
      output_res <- read.table("./Figure 3/Medip_bin/output/Medip_STADvsNC_edger_exact_bin.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      STADvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      write.table(output_res[STADvsNC_features,],"./Figure 3/Medip_bin/PCA/STADvsNC_features.txt",sep = "\t",quote = FALSE)
      output_res <- read.table("./Figure 3/Medip_bin/output/Medip_STADvsNC_edger_exact_bin.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
      CRCvsNC_features <- rownames(output_res[output_res$pvalue<0.05 & sqrt(output_res$log2FoldChange*output_res$log2FoldChange) > 0.59,])
      write.table(output_res[CRCvsNC_features,],"./Figure 3/Medip_bin/PCA/CRCvsNC_features.txt",sep = "\t",quote = FALSE)
      top_features <-  unique(c(STADvsNC_features,CRCvsNC_features))
      
      #gini
      {
        Medip_bin <- Medip_bin[top_features,]
        
        STAD_gini <- data.frame(gini=gini(t(Medip_bin[,grep("STAD",colnames(Medip_bin))])),group="STAD")
        CRC_gini <- data.frame(gini=gini(t(Medip_bin[,grep("CRC",colnames(Medip_bin))])),group="CRC")
        NC_gini <- data.frame(gini=gini(t(Medip_bin[,grep("NC",colnames(Medip_bin))])),group="HD")
        gini_forplot <- rbind(NC_gini,STAD_gini,CRC_gini)
        
        STAD_NC <- ks.test(STAD_gini$gini,NC_gini$gini)
        CRC_NC <- ks.test(CRC_gini$gini,NC_gini$gini)
        STAD_CRC <- ks.test(STAD_gini$gini,CRC_gini$gini)
        
        ks_test <- data.frame("group"=c("STAD vs. HD","CRC vs. HD","CRC vs. STAD"),
                              "pvalue"=c(STAD_NC$p.value,CRC_NC$p.value,STAD_CRC$p.value),
                              "method"=c(STAD_NC$method,CRC_NC$method,STAD_CRC$method),
                              "alternative"=c(STAD_NC$alternative,CRC_NC$alternative,STAD_CRC$alternative))
        
        write.csv(ks_test,"./Figure 3/Medip_bin/gini/ks_test.csv")
        
        p <- ggplot(gini_forplot, aes(x=gini,..scaled.., fill=group))+
          geom_density(color="black")+
          #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
          scale_fill_manual(values=Group_color)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c("0","0.2","0.4","0.6","0.8","1"),expand = c(0,0),limits = c(0,1))+
          xlab("gini index")+
          ylab("Density")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        ggsave("./Figure 3/Medip_bin/gini/gini_density.pdf",p,height = 2.46,width = 4.3, device = "pdf")
      }
      
      #PCA
      {
        matrix_forPCA <- as.matrix(t(na.omit(full_matrix[top_features,as.character(na.omit(ML_samples$Methylation_id))])))
        matrix_forPCA <- t(scale(t(matrix_forPCA)))
        
        group <- as.character(lapply(strsplit(as.character(rownames(matrix_forPCA)),"-",fixed = TRUE),function(x) x[1]))
        group <- gsub("NC","HD",group)
        
        PCA.res <- PCA(matrix_forPCA, scale.unit = TRUE, ncp = 3, graph = FALSE)
        
        for_distance <- as.data.frame(PCA.res$ind$coord)
        for_distance$Group <- as.character(lapply(strsplit(rownames(for_distance),"-"),function(x) x[1]))
        Group1 <- "STAD"
        Group2 <- "NC"
        STAD_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "CRC"
        Group2 <- "NC"
        CRC_NC_out <- PCA_group_distance(for_distance,Group1,Group2)
        Group1 <- "STAD"
        Group2 <- "CRC"
        STAD_CRC_out <- PCA_group_distance(for_distance,Group1,Group2)
        PCA_distance <- unique(rbind(STAD_NC_out,CRC_NC_out,STAD_CRC_out))
        write.csv(PCA_distance,"./Figure 3/Medip_bin/PCA/PCA_distance.csv")
        
        pdf(file=paste0("./Figure 3/Medip_bin/PCA/",comparison_type,"_",date(),".pdf"), height = 4 ,width = 4)
        fviz_pca_ind(PCA.res,
                     geom.ind = "point",
                     palette = Group_color,#c(CRC="#FCB514",STAD="red",HD="blue"), #"jco",
                     col.ind = group,
                     addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.999999999999999,axes.linetype=NA)+
          #labs(title = paste0("cfRNA expression"," | ",comparison_type))+
          labs(title = "")+
          theme_bw()+
          theme(
            #plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=20,face="bold",vjust=0),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            #axis.line = element_line(color = "black"),
            axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=16),
            axis.title.y = element_text(face="bold",color="black", size=16))
        dev.off()
      }
    }
  }
  
  #PCA summary
  {
    PCA_summary <- read.csv("Figure 3/PCA_distance_summary.csv")
    PCA_summary$Average.Distance.between.Cancer.and.HD <- (PCA_summary$Distance.between.CRC.and.HD+PCA_summary$Distance.between.STAD.and.HD)/2
    PCA_summary_plot <- reshape2::melt(PCA_summary[,c("Group","Average.Distance.between.Cancer.and.HD","Distance.between.STAD.and.CRC")])
    colnames(PCA_summary_plot) <- c("Alteration","Group","Distance")
    PCA_summary_plot$Group <- gsub("Average.Distance.between.Cancer.and.HD","Cancer vs. HD", PCA_summary_plot$Group)
    PCA_summary_plot$Group <- gsub("Distance.between.STAD.and.CRC","STAD vs. CRC", PCA_summary_plot$Group)
    
    #PCA_summary_plot <- reshape2::melt(PCA_summary[,c("Group","Distance.between.STAD.and.HD","Distance.between.CRC.and.HD","Distance.between.STAD.and.CRC")])
    #colnames(PCA_summary_plot) <- c("Alteration","Group","Distance")
    #PCA_summary_plot$Group <- gsub("Distance.between.STAD.and.HD","STAD vs. HD", PCA_summary_plot$Group)
    #PCA_summary_plot$Group <- gsub("Distance.between.CRC.and.HD","CRC vs. HD", PCA_summary_plot$Group)
    #PCA_summary_plot$Group <- gsub("Distance.between.STAD.and.CRC","STAD vs. CRC", PCA_summary_plot$Group)
    
    PCA_summary_plot$distant <- NA
    PCA_summary_plot[PCA_summary_plot$Distance>=1.2,]$distant <- "Yes"
    PCA_summary_plot[PCA_summary_plot$Distance<1.2,]$distant <- "No"
    
    #rank <- as.character(PCA_summary_plot[order(PCA_summary_plot[PCA_summary_plot$Group == "CRC vs. HD",]$Distance,decreasing = TRUE),]$Alteration)
    rank <- c("Medip (promoter)","Medip (CpG island)","WPS","EndMotif","FragSize","CNV","miRNA","Expression","Splicing","chimeric","Alternative promoter","Editing","ASE","Mutation")
    PCA_summary_plot$Alteration <- factor(PCA_summary_plot$Alteration, level = rank)
    PCA_summary_plot$Group <- factor(PCA_summary_plot$Group, level = c("Cancer vs. HD","CRC vs. HD","STAD vs. HD","STAD vs. CRC"))
    p <- ggplot(PCA_summary_plot, aes(x=Alteration,y=Distance))+
      geom_bar(aes(fill = distant),stat = "identity",position = "dodge")+
      facet_wrap(~Group,nrow = 2,strip.position = "right")+
      geom_hline(yintercept = c(1.2),linetype = 3)+
      #geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
      scale_fill_manual(values=c("Yes"="red","No"="light grey"))+
      scale_y_continuous(limits = c(0,3.5),expand = c(0,0))+
      xlab("")+
      ylab("Out-class/In-class distance")+
      theme_bw()+
      theme(
        plot.margin = unit(x=c(5,5,5,30),units="pt"),
        legend.position="none",
        strip.background = element_rect(fill = "white",color = "white"),
        strip.text = element_text(hjust = 0.5,size=12),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_text(color="black", size=20),
        axis.title.y = element_text(color="black", size=16))
    p
    ggsave("./Figure 3/Distance_merged.pdf",p,height = 4.46,width = 5.54, device = "pdf")
  }

  
  #illustration of different features
  {
    #chimeric
    {
    devtools::install_github("stianlagstad/chimeraviz")

    # Load chimeraviz
    library(chimeraviz)
    
    starfusionData <- read.table("/Users/yuhuan/Downloads/star-fusion.fusion_predictions.abridged.tsv",comment.char = "",sep = "\t",header = TRUE,check.names = FALSE)
    
    
    starFusion <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/15.examples/chimeric/star-fusion.fusion_predictions.abridged.tsv"
    fusions <- import_starfusion(starFusion, "hg38", 10)
    #create_fusion_report(fusions, "output.html")
    plot_circle(fusions)
    #?RCircos.Reset.Plot.Parameters
    
    
    library(Rsamtools)
    sortBam("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/15.examples/chimeric/CRC-PKU-12-picoAligned.out.reheader.bam",destination = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/15.examples/chimeric/CRC-PKU-12-picoAligned.out.sorted")
    indexBam("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/15.examples/chimeric/CRC-PKU-12-picoAligned.out.sorted.bam")
    starBam <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/15.examples/chimeric/CRC-PKU-12-picoAligned.out.sorted.bam"     
    
    
    
    #if (!require("BiocManager", quietly = TRUE))
    #  install.packages("BiocManager")
    #BiocManager::install("EnsDb.Hsapiens.v86")
    library(EnsDb.Hsapiens.v86)
    
    ## Making a "short cut"
    edb <- EnsDb.Hsapiens.v86
    ## print some informations for this package
    edb
    
    plot_circle(fusions)
    fusion <- get_fusion_by_id(fusions, 2)
    plot_fusion(fusion = fusion,
                ylim = c(0, 20),
                bamfile = starBam,
                edb = edb)
    
    plot_fusion_transcript(
      fusion = fusion,
      edb = edb,
      bamfile = starBam)
    
    seqnamesTabix(starBam)
    }
    
    #Mutation
    {
      Coding_consequence <- as.data.frame(t(read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 2/SNP/Variant_class.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)))
      Coding_consequence <- Coding_consequence[-grep("mix.-pico",rownames(Coding_consequence)),]
      rownames(Coding_consequence) <- gsub("NC","HD",rownames(Coding_consequence))
      Coding_consequence$Group <- as.character(lapply(strsplit(rownames(Coding_consequence),"-"), function(x) x[1]))
      
      Mean <- aggregate(Coding_consequence[,-which(colnames(Coding_consequence)=="Group")], by = list("Group"=Coding_consequence$Group),FUN = mean)
      SD <- aggregate(Coding_consequence[,-which(colnames(Coding_consequence)=="Group")], by = list("Group"=Coding_consequence$Group),FUN = sd)
    
      Mean_plot <- melt(Mean)
      colnames(Mean_plot) <- gsub("value","Mean",colnames(Mean_plot))
      SD_plot <- melt(SD)
      
      Mean_plot$SD <- SD_plot$value
      Mean_plot$Group <- factor(Mean_plot$Group,levels=c("HD","CRC","STAD"))
      ggplot(Mean_plot,aes(x=variable,y=Mean,fill=Group))+
        #facet_wrap(~metric,nrow = 4)+
        geom_bar(stat = "identity",colour = "black",position = "dodge")+
        #geom_text(aes(x=variable,y=Mean+1,label=round(Mean,digits = 0)),size = 4,angle = 0,position = "dodge")+
        geom_errorbar(aes(x=variable, ymin=Mean-SD, ymax=Mean+SD), width=0.9, colour="black", alpha=0.5, size=1,position = "dodge")+
        xlab("")+
        labs(title = "")+
        #scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("0","0.2","0.4","0.6","0.8","1.0"),expand = c(0,0),limits = c(0,1))+
        #scale_fill_aaas()+
        scale_fill_manual(values = Group_color)+
        #geom_vline(xintercept=1.5,linetype = 2)+
        #geom_vline(xintercept=3.5,linetype = 2)+
        xlab("")+
        ylab("Number")+
        theme_bw()+
        theme(
          plot.margin = unit(x=c(15,5,-10,70),units="pt"),
          legend.position="none",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face="bold", color="black", size=18),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))
    }
  }
}

#Figure 3
{
  #Differential pathway number
  {
    comparison <- "STADvsCRC"
    {
      Differential_expression <- read.csv(paste0("Expression/output/Expression_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      Expression_up_gene <- length(rownames(Differential_expression[(Differential_expression$GeneEnrichedIn=="Up regulated")&(Differential_expression$p.adjust < 0.05),]))
      Expression_down_gene <- length(rownames(Differential_expression[(Differential_expression$GeneEnrichedIn=="Down regulated")&(Differential_expression$p.adjust < 0.05),]))
      
      Differential_miRNA <- read.csv(paste0("miRNA/output/miRNA_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      miRNA_up_gene <- length(rownames(Differential_miRNA[(Differential_miRNA$GeneEnrichedIn=="Up regulated")&(Differential_miRNA$p.adjust < 0.05),]))
      miRNA_down_gene <- length(rownames(Differential_miRNA[(Differential_miRNA$GeneEnrichedIn=="Down regulated")&(Differential_miRNA$p.adjust < 0.05),]))
      
      Differential_Splicing <- read.csv(paste0("Splicing/output/Splicing_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      Splicing_up_gene <- length(rownames(Differential_Splicing[(Differential_Splicing$GeneEnrichedIn=="Up regulated")&(Differential_Splicing$p.adjust < 0.05),]))
      Splicing_down_gene <- length(rownames(Differential_Splicing[(Differential_Splicing$GeneEnrichedIn=="Down regulated")&(Differential_Splicing$p.adjust < 0.05),]))
      
      Differential_Editing <- read.csv(paste0("Editing/output/Editing_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      Editing_up_site <- length(rownames(Differential_Editing[(Differential_Editing$GeneEnrichedIn=="Up regulated")&(Differential_Editing$p.adjust < 0.05),]))
      Editing_down_site <- length(rownames(Differential_Editing[(Differential_Editing$GeneEnrichedIn=="Down regulated")&(Differential_Editing$p.adjust < 0.05),]))
      
      Differential_Mutation <- read.csv(paste0("Mutation/output/Mutation_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      Mutation_up_site <- length(rownames(Differential_Mutation[(Differential_Mutation$GeneEnrichedIn=="Up regulated")&(Differential_Mutation$p.adjust < 0.05),]))
      Mutation_down_site <- length(rownames(Differential_Mutation[(Differential_Mutation$GeneEnrichedIn=="Down regulated")&(Differential_Mutation$p.adjust < 0.05),]))
      
      Differential_AlternativePromoter <- read.csv(paste0("Alternative promoter/output/AlternativePromoter_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      AlternativePromoter_up_gene <- length(rownames(Differential_AlternativePromoter[(Differential_AlternativePromoter$GeneEnrichedIn=="Up regulated")&(Differential_AlternativePromoter$p.adjust < 0.05),]))
      AlternativePromoter_down_gene <- length(rownames(Differential_AlternativePromoter[(Differential_AlternativePromoter$GeneEnrichedIn=="Down regulated")&(Differential_AlternativePromoter$p.adjust < 0.05),]))
      
      #Differential_ASE <- read.csv(paste0("ASE/output/ASE_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      #ASE_up_site <- length(rownames(Differential_ASE[(Differential_ASE$GeneEnrichedIn=="Up regulated")&(Differential_ASE$p.adjust < 0.05),]))
      #ASE_down_site <- length(rownames(Differential_ASE[(Differential_ASE$GeneEnrichedIn=="Down regulated")&(Differential_ASE$p.adjust < 0.05),]))
      ASE_up_site <- 0
      ASE_down_site <- 0
      
      #Differential_chimeric <- read.csv(paste0("chimeric/output/chimeric_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      #chimeric_up_gene <- length(rownames(Differential_chimeric[(Differential_chimeric$GeneEnrichedIn=="Up regulated")&(Differential_chimeric$p.adjust < 0.05),]))
      #chimeric_down_gene <- length(rownames(Differential_chimeric[(Differential_chimeric$GeneEnrichedIn=="Down regulated")&(Differential_chimeric$p.adjust < 0.05),]))
      chimeric_up_gene <- 0
      chimeric_down_gene <- 0
      
      Differential_CNV <- read.csv(paste0("CNV/output/CNV_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE)
      CNV_up_gene <- length(rownames(Differential_CNV[(Differential_CNV$GeneEnrichedIn=="Up regulated")&(Differential_CNV$p.adjust < 0.05),]))
      CNV_down_gene <- length(rownames(Differential_CNV[(Differential_CNV$GeneEnrichedIn=="Down regulated")&(Differential_CNV$p.adjust < 0.05),]))
      #CNV_up_gene <- 0
      #CNV_down_gene <- 0
      
      Differential_Methylation <- read.csv(paste0("Medip/output/Methylation_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      Methylation_up_gene <- length(rownames(Differential_Methylation[(Differential_Methylation$GeneEnrichedIn=="Up regulated")&(Differential_Methylation$p.adjust < 0.05),]))
      Methylation_down_gene <- length(rownames(Differential_Methylation[(Differential_Methylation$GeneEnrichedIn=="Down regulated")&(Differential_Methylation$p.adjust < 0.05),]))
      
      Differential_WPS <- read.csv(paste0("WPS/output/WPS_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      WPS_up_gene <- length(rownames(Differential_WPS[(Differential_WPS$GeneEnrichedIn=="Up regulated")&(Differential_WPS$p.adjust < 0.05),]))
      WPS_down_gene <- length(rownames(Differential_WPS[(Differential_WPS$GeneEnrichedIn=="Down regulated")&(Differential_WPS$p.adjust < 0.05),]))
      #Nucleosome_up_gene <- 0
      #Nucleosome_down_gene <- 0
      
      #Differential_FragSize <- read.csv(paste0("FragSize/output/FragSize_",comparison,"_KEGG.txt"),sep = "\t", header = TRUE, row.names = 1)
      #FragSize_up_gene <- length(rownames(Differential_FragSize[(Differential_Nucleosome$GeneEnrichedIn=="Up regulated")&(Differential_Nucleosome$p.adjust < 0.05),]))
      #FragSize_down_gene <- length(rownames(Differential_FragSize[(Differential_Nucleosome$GeneEnrichedIn=="Down regulated")&(Differential_Nucleosome$p.adjust < 0.05),]))
      FragSize_up_gene <- 0
      FragSize_down_gene <- 0
      
      
      up_events <- data.frame("Alteration"=c("cfRNA abundance","cf-miRNA abundance","cfRNA splicing","cfRNA editing","cfRNA SNV","cfRNA alternative promoter usage","cfRNA allele specific expression","Chimeric cfRNA","cfDNA copy number","cfDNA methylation (promoter)","cfDNA window protection score","cfDNA fragment size"),"Number"=c(Expression_up_gene,miRNA_up_gene,Splicing_up_gene,Editing_up_site,Mutation_up_site,AlternativePromoter_up_gene,ASE_up_site,chimeric_up_gene,CNV_up_gene,Methylation_up_gene,WPS_up_gene,FragSize_up_gene))
      down_events <- data.frame("Alteration"=c("cfRNA abundance","cf-miRNA abundance","cfRNA splicing","cfRNA editing","cfRNA SNV","cfRNA alternative promoter usage","cfRNA allele specific expression","Chimeric cfRNA","cfDNA copy number","cfDNA methylation (promoter)","cfDNA window protection score","cfDNA fragment size"),"Number"=c(-Expression_down_gene,-miRNA_down_gene,-Splicing_down_gene,-Editing_down_site,-Mutation_down_site,-AlternativePromoter_down_gene,-ASE_down_site,-chimeric_down_gene,-CNV_down_gene,-Methylation_down_gene,-WPS_down_gene,-FragSize_down_gene))
    }
    RNA_color <- c("cfRNA abundance"="#5CACEE","cfRNA alternative promoter usage"="#5CACEE","cfRNA allele specific expression"="#5CACEE","Chimeric cfRNA"="#5CACEE",
                   "cfRNA editing"="#5CACEE","cfRNA abundance (microbe genus abundance)"="#5CACEE","cfRNA SNV"="#5CACEE","cfRNA splicing"="#5CACEE","cfRNA abudance (transposable elements)"="#5CACEE",
                   "cf-miRNA abundance"=alpha("#0a14a6",alpha=1),
                   "cf-miRNA abundance"=alpha("#0a14a6",alpha=1),
                   "cfDNA methylation (promoter)"="#666666","cfDNA methylation (CpG island)"="#666666",
                   "cfDNA copy number"="#EEB4B4","cfDNA window protection score"="#EEB4B4","cfDNA end motif usage"="#EEB4B4","cfDNA fragment size"="#EEB4B4"
    )
    #scale_fill_manual(values=c(alpha("#5CACEE",alpha = 1),alpha("#EEB4B4",alpha = 0.45),alpha("#666666",alpha = 0.1),alpha("#CD0000",alpha = 0.8)))
    differential_events <- rbind(up_events,down_events)
    
    
    #differential_events$Alteration <- gsub("cf-miRNA abundance","cf-miRNA targeted gene",differential_events$Alteration)
    #up_events$Alteration <- gsub("cf-miRNA abundance","cf-miRNA targeted gene",up_events$Alteration)
    #differential_events$Alteration <- factor(differential_events$Alteration, levels = rev(c("miRNA targeted gene","DNA copy number","Allele specific expression",
    #                                                                                        "RNA expression","Alternative promoter","RNA SNP",
    #                                                                                        "RNA splicing","DNA nucleosome occupancy","RNA editing",
    #                                                                                        "DNA methylation","Chimeric RNA")))
    differential_events$Alteration <- factor(differential_events$Alteration, levels = rev(rev(up_events[order(up_events$Number-down_events$Number,decreasing = TRUE),]$Alteration)))
    p <- ggplot(differential_events,aes(x=Alteration,y=Number,fill=Alteration))+
      geom_bar(stat = "identity",colour = "black")+
      scale_fill_manual(values = RNA_color)+
      geom_hline(aes(yintercept = 0))+
      coord_flip()+
      xlab("")+
      ylab("")+
      scale_y_continuous(breaks = c(-30,-20,-10,-5,0,5,10,20,30),labels = c("30","20","10","5","0","5","10","20","30"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(15,5,10,45),units="pt"),
        legend.position='none',
        panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        #axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_blank(),
        #axis.line.y = element_line(color = "black"),
        axis.line = element_blank(),
        axis.text.x = element_text(color="black", size=16, angle = 45,hjust = 1,vjust = 1),
        #axis.text.y = element_blank(),
        axis.text.y = element_text( color="black", size=16, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))
    ggsave(paste0("./pathway_",comparison,".pdf"),p,width = 6.9, height = 4.14,device = "pdf")
  }
  
  #KEGG top3 pathway plot
  {
    KEGG_pathway_categroy <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/data/KEGG_pathway_categroy.csv")
    KEGG <- "KEGG"
    comparison <- "GIvsNC"
    {
      KEGG_expression <- read.csv(paste0("Expression/output/Expression_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE, row.names = 1)
      KEGG_expression_enriched <- KEGG_expression[KEGG_expression$p.adjust < 0.05,]
      KEGG_expression_enriched$Alteration <- "cfRNA abundance"
      
      KEGG_miRNA_target <- read.csv(paste0("miRNA/output/miRNA_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE, row.names = 1)
      KEGG_miRNA_target_enriched <- KEGG_miRNA_target[KEGG_miRNA_target$p.adjust < 0.05,]
      KEGG_miRNA_target_enriched$Alteration <- "cf-miRNA abundance"
      
      KEGG_Splicing <- read.csv(paste0("Splicing/output/Splicing_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE, row.names = 1)
      KEGG_Splicing_enriched <- KEGG_Splicing[KEGG_Splicing$p.adjust < 0.05,]
      KEGG_Splicing_enriched$Alteration <- "cfRNA splicing"
      
      #KEGG_Editing <- read.csv(paste0("Editing/output/Editing_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE, row.names = 1)
      #KEGG_Editing_enriched <- KEGG_Editing[KEGG_Editing$p.adjust < 0.05,]
      #KEGG_Editing_enriched$Alteration <- "RNA editing"
      
      KEGG_Mutation <- read.csv(paste0("Mutation/output/Mutation_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE, row.names = 1)
      KEGG_Mutation_enriched <- KEGG_Mutation[KEGG_Mutation$p.adjust < 0.05,]
      KEGG_Mutation_enriched$Alteration <- "cfRNA SNV"
      
      #KEGG_AlternativePromoter <- read.csv(paste0("Alternative promoter/output/AlternativePromoter_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE, row.names = 1)
      #KEGG_AlternativePromoter_enriched <- KEGG_AlternativePromoter[KEGG_AlternativePromoter$p.adjust < 0.05,]
      #KEGG_AlternativePromoter_enriched$Alteration <- "Alternative promoter"
      
      #KEGG_ASE <- read.csv(paste0("ASE/output/ASE_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE, row.names = 1)
      #KEGG_ASE_enriched <- KEGG_ASE[KEGG_ASE$p.adjust < 0.05,]
      #KEGG_ASE_enriched$Alteration <- "Allele specific expression"
      
      #KEGG_chimeric <- read.csv(paste0("chimeric/output/chimeric_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE, row.names = 1)
      #KEGG_chimeric_enriched <- KEGG_chimeric[KEGG_chimeric$p.adjust < 0.05,]
      #KEGG_chimeric_enriched$Alteration <- "Chimeric RNA"
      
      KEGG_CNV <- read.csv(paste0("CNV/output/CNV_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE)
      KEGG_CNV_enriched <- KEGG_CNV[KEGG_CNV$p.adjust < 0.05,]
      KEGG_CNV_enriched$Alteration <- "cfDNA copy number"
      
      KEGG_Methylation <- read.csv(paste0("Medip/output/Methylation_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE)
      KEGG_Methylation_enriched <- KEGG_Methylation[KEGG_Methylation$p.adjust < 0.05,]
      KEGG_Methylation_enriched$Alteration <- "cfDNA methylation (promoter)"
      
      KEGG_WPS <- read.csv(paste0("WPS/output/WPS_",comparison,"_",KEGG,".txt"),sep = "\t", header = TRUE)
      KEGG_WPS_enriched <- KEGG_WPS[KEGG_WPS$p.adjust < 0.05,]
      KEGG_WPS_enriched$Alteration <- "cfDNA window protection score"
      
      KEGG_enriched <- rbind(KEGG_expression_enriched,KEGG_miRNA_target_enriched,
                             KEGG_Splicing_enriched,#KEGG_Editing_enriched,
                             KEGG_Mutation_enriched,#KEGG_AlternativePromoter_enriched,
                             #KEGG_ASE_enriched,KEGG_chimeric_enriched,
                             KEGG_CNV_enriched,KEGG_Methylation_enriched,KEGG_WPS_enriched)
    }
    
    #Up pathways
    {
      top <- 10
      KEGG_res <- KEGG_expression
      filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_expression_up <- factor(up,levels=up)
      
      KEGG_res <- KEGG_miRNA_target
      filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_miRNA_target_up <- factor(up,levels=up)
      
      KEGG_res <- KEGG_Splicing
      filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_Splicing_up <- factor(up,levels=up)
      
      #KEGG_res <- KEGG_Editing
      #filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      #up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      #KEGG_Editing_up <- factor(up,levels=up)
      
      KEGG_res <- KEGG_Mutation
      filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_Mutation_up <- factor(up,levels=up)
      
      #KEGG_res <- KEGG_AlternativePromoter
      #filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      #up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      #KEGG_AlternativePromoter_up <- factor(up,levels=up)
      
      #KEGG_res <- KEGG_ASE
      #filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      #up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      #KEGG_ASE_up <- factor(up,levels=up)
      
      #KEGG_res <- KEGG_chimeric
      #filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      #up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      #KEGG_chimeric_up <- factor(up,levels=up)
      
      KEGG_res <- KEGG_CNV
      filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_CNV_up <- factor(up,levels=up)
      
      KEGG_res <- KEGG_Methylation
      filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_Methylation_up <- factor(up,levels=up)
      
      KEGG_res <- KEGG_WPS
      filtered_up <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Up regulated"), ]
      up <- head(filtered_up[order(filtered_up$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_WPS_up <- factor(up,levels=up)
      
      Up_pathways <- fct_c(KEGG_expression_up,KEGG_miRNA_target_up,KEGG_Splicing_up,
                           #KEGG_Editing_up,
                           KEGG_Mutation_up,
                           #KEGG_AlternativePromoter_up,KEGG_ASE_up,KEGG_chimeric_up,
                           KEGG_CNV_up,KEGG_Methylation_up,KEGG_WPS_up)
      
      KEGG_Up <- KEGG_enriched[which(KEGG_enriched$GeneEnrichedIn=="Up regulated"),]
      
      KEGG_Up <- filter(KEGG_Up, Description %in% Up_pathways)
      
      #KEGG_Up <- KEGG_Up[-which(KEGG_Up$Description=="Coronavirus disease - COVID-19"),]
      #KEGG_Up <- KEGG_Up[-which(KEGG_Up$Alteration=="RNA editing"),]
      KEGG_Up$Alteration <- factor(KEGG_Up$Alteration, levels = c("cfDNA copy number","cfDNA methylation (promoter)","cfDNA window protection score","cfRNA abundance","cf-miRNA abundance",
                                                                  "cfRNA splicing","cfRNA SNV"))
      up <- ggplot(KEGG_Up,aes(x=Alteration,y=Description))+
        geom_point(aes(size=-1*log10(p.adjust),color=as.numeric(Count)))+
        scale_colour_gradient2(low="#162252",high="#3B5E2B",mid="#008B00",midpoint = 0)+
        #geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
        #geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
        labs(color="GeneCount",
             size=expression(-log[10](p.adjust)),
             x="")+
        theme_bw()+
        theme(axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
              axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black",angle = 45,vjust=1,hjust=1),
              axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
              axis.title.y = element_blank())
    }
    
    #Down pathways
    {
      top <- 10
      KEGG_res <- KEGG_expression
      KEGG_res <- KEGG_res[-which(KEGG_res$Description=="Coronavirus disease - COVID-19"),]
      filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_expression_Down <- factor(Down,levels=Down)
      
      KEGG_res <- KEGG_miRNA_target
      filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_miRNA_target_Down <- factor(Down,levels=Down)
      
      KEGG_res <- KEGG_Splicing
      filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_Splicing_Down <- factor(Down,levels=Down)
      
      #KEGG_res <- KEGG_Editing
      #filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      #Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      #KEGG_Editing_Down <- factor(Down,levels=Down)
      
      KEGG_res <- KEGG_Mutation
      filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_Mutation_Down <- factor(Down,levels=Down)
      
      #KEGG_res <- KEGG_AlternativePromoter
      #KEGG_res <- KEGG_res[-which(KEGG_res$Description=="Coronavirus disease - COVID-19"),]
      #filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      #Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      #KEGG_AlternativePromoter_Down <- factor(Down,levels=Down)
      
      #KEGG_res <- KEGG_ASE
      #filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      #Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      #KEGG_ASE_Down <- factor(Down,levels=Down)
      
      #KEGG_res <- KEGG_chimeric
      #filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      #Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      #KEGG_chimeric_Down <- factor(Down,levels=Down)
      
      KEGG_res <- KEGG_CNV
      filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_CNV_Down <- factor(Down,levels=Down)
      
      KEGG_res <- KEGG_Methylation
      filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_Methylation_Down <- factor(Down,levels=Down)
      
      KEGG_res <- KEGG_WPS
      filtered_Down <- KEGG_res[(KEGG_res$p.adjust < 0.05) & (KEGG_res$GeneEnrichedIn=="Down regulated"), ]
      Down <- head(filtered_Down[order(filtered_Down$p.adjust,decreasing = FALSE),]$Description,top)
      KEGG_WPS_Down <- factor(Down,levels=Down)
      
      Down_pathways <- fct_c(KEGG_expression_Down,KEGG_miRNA_target_Down,KEGG_Splicing_Down,
                             #KEGG_Editing_Down,
                             KEGG_Mutation_Down, 
                             #KEGG_AlternativePromoter_Down,KEGG_ASE_Down,KEGG_chimeric_Down,
                             KEGG_CNV_Down,KEGG_Methylation_Down,KEGG_WPS_Down)
      
      KEGG_Down <- KEGG_enriched[which(KEGG_enriched$GeneEnrichedIn=="Down regulated"),]
      
      KEGG_Down <- filter(KEGG_Down, Description %in% Down_pathways)
      
      #KEGG_Down <- KEGG_Down[-which(KEGG_Down$Description=="Coronavirus disease - COVID-19"),]
      KEGG_Down$Alteration <- factor(KEGG_Down$Alteration, levels = c("cfDNA copy number","cfDNA methylation (promoter)","cfDNA window protection score","cfRNA abundance","cf-miRNA abundance",
                                                                      "cfRNA splicing","cfRNA SNV"))
      down <- ggplot(KEGG_Down,aes(x=Alteration,y=Description))+
        geom_point(aes(size=-1*log10(p.adjust),color=as.numeric(Count)))+
        scale_colour_gradient2(low="#162252",high="#3B5E2B",mid="#008B00",midpoint = 0)+
        #geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
        #geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
        labs(color="GeneCount",
             size=expression(-log[10](p.adjust)),
             x="")+
        theme_bw()+
        theme(axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
              axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black",angle = 45,vjust=1,hjust=1),
              axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
              axis.title.y = element_blank())
    }
    KEGG <- rbind(KEGG_Up,KEGG_Down)
    KEGG$Description <- factor(KEGG$Description,levels = rev(KEGG_pathway_categroy$Description))
    
    p <- ggplot(KEGG,aes(x=Alteration,y=Description))+
      facet_grid(.~KEGG$GeneEnrichedIn)+
      geom_point(aes(size=-1*log10(p.adjust),color=as.numeric(Count)))+
      scale_colour_gradient2(low="#162252",high="#3B5E2B",mid="#008B00",midpoint = 0)+
      #geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
      #geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
      labs(color="GeneCount",
           size=expression(-log[10](p.adjust)),
           x="")+
      theme_bw()+
      geom_hline(aes(yintercept = 45.5))+
      geom_hline(aes(yintercept = 37.5))+
      geom_hline(aes(yintercept = 34.5))+
      geom_hline(aes(yintercept = 30.5))+
      geom_hline(aes(yintercept = 28.5))+
      geom_hline(aes(yintercept = 25.5))+
      geom_hline(aes(yintercept = 20.5))+
      geom_hline(aes(yintercept = 13.5))+
      theme(strip.text = element_text(size = rel(1.3),colour = "black"),
            axis.text.y = element_text(size = rel(1.3),colour = "black"),
            axis.text.x = element_text(size=rel(1.3),colour = "black",angle = 45,vjust=1,hjust=1),
            axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
            axis.title.y = element_blank())
    p
    ggsave(plot = p, filename = "./KEGG_bubble.pdf", width = 9, height = 12,device = "pdf")
  }
  
  #ActivePathway
  {}
  
  #pathway boxplot
  {
    #make candidate list
    {
      gene_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/T cell receptor signaling pathway.csv",header = TRUE)
      values <- as.character(lapply(strsplit(as.character(gene_list$Gene.ID),".",fixed = TRUE),function(x) x[1]))
      mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
      candidate_list_raw <- getBM(attributes=c("hgnc_symbol","description","ensembl_gene_id","chromosome_name","start_position","end_position","strand"), 
                                  filters = "ensembl_gene_id", 
                                  values=unique(values), mart= mart,useCache = FALSE)
      candidate_list_raw$strand <- gsub("-1","-",fixed = TRUE,candidate_list_raw$strand)
      candidate_list_raw$strand <- gsub("1","+",fixed = TRUE,candidate_list_raw$strand)
      candidate_list_raw <- candidate_list_raw[-grep("_",fixed=TRUE,candidate_list_raw$chromosome_name),]
      candidate_list_raw$chromosome_name <- paste0("chr",candidate_list_raw$chromosome_name)
      colnames(candidate_list_raw) <- gsub("hgnc_symbol","Gene",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("description","Description",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("ensembl_gene_id","ensembl",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("chromosome_name","chromosome",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("start_position","start",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("end_position","end",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("strand","strand",fixed = TRUE,colnames(candidate_list_raw))
      write.csv(candidate_list_raw,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/T cell receptor signaling pathway.csv",row.names = FALSE)
    }
    {
      gene_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP.csv",header = TRUE)
      values <- as.character(lapply(strsplit(as.character(gene_list$GeneName),".",fixed = TRUE),function(x) x[1]))
      gene_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/23.TCGA and GTEx/STAD_colon_coupregulated_genes.csv",header = TRUE)
      values <- as.character(lapply(strsplit(as.character(gene_list$hgnc_symbol),".",fixed = TRUE),function(x) x[1])) 
      #gene_list <- gene_list[which(gene_list$Super.Category=="Antigen presentation"),]
      
      
      mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
      candidate_list_raw <- getBM(attributes=c("hgnc_symbol","description","ensembl_gene_id","chromosome_name","start_position","end_position","strand"), 
                                  filters = "hgnc_symbol", 
                                  values=unique(values), mart= mart,useCache = FALSE)
      candidate_list_raw$strand <- gsub("-1","-",fixed = TRUE,candidate_list_raw$strand)
      candidate_list_raw$strand <- gsub("1","+",fixed = TRUE,candidate_list_raw$strand)
      candidate_list_raw <- candidate_list_raw[-grep("_",fixed=TRUE,candidate_list_raw$chromosome_name),]
      candidate_list_raw$chromosome_name <- paste0("chr",candidate_list_raw$chromosome_name)
      colnames(candidate_list_raw) <- gsub("hgnc_symbol","Gene",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("description","Description",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("ensembl_gene_id","ensembl",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("chromosome_name","chromosome",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("start_position","start",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("end_position","end",fixed = TRUE,colnames(candidate_list_raw))
      colnames(candidate_list_raw) <- gsub("strand","strand",fixed = TRUE,colnames(candidate_list_raw))
      write.csv(candidate_list_raw,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/STAD_colon_coupregulated_genes.csv",row.names = FALSE)
    }
    #Plot siganture score of immune pathways
    {
      pathway <- "GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP"
      pathway <- "GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP"
      pathway <- "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP"
      pathway <- "T cell receptor signaling pathway"
      pathway <- "B cell receptor signaling pathway"
      pathway <- "Inhibotory immune genes"
      pathway <- "PELO"
      pathway <- "STAD_colon_coupregulated_genes"
      i=1
      
      prefix <- "Multiomics_Plasma"
      counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE174302_intron-spanning-available-samples_TPM.txt",sep = "\t", header = TRUE,check.names = FALSE)
      counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE133684_count_matrix_TPM.txt",sep = "\t", header = TRUE,check.names = FALSE)
      counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/multiomics_plasma_TPM_intron_spanning_20220427.txt",sep = "\t", header = TRUE,check.names = FALSE, row.names = 1)
      counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Plasma_TPM.txt",sep = "\t", header = TRUE,check.names = FALSE, row.names = 1)
      rownames(counts) <- as.character(lapply(strsplit(as.character(rownames(counts)),".",fixed = TRUE),function(x) x[1]))
      counts$Gene <- as.character(lapply(strsplit(as.character(rownames(counts)),".",fixed = TRUE),function(x) x[1]))
      counts <- aggregate(counts[,-which(colnames(counts)=="Gene")], list(Gene=counts$Gene), FUN = sum)
      rownames(counts) <- counts$Gene
      counts <- counts[,-which(colnames(counts)=="Gene")]
      #counts <- counts[,-which(colnames(counts)=="CRC-PKU-29-PBMC")]
      
      prefix <- "Multiomics_20211113"
      #counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/Expression/matrix/20211113_TPM.txt",sep = "\t", header = TRUE,check.names = FALSE)
      #counts$`Gene|Length` <- as.character(lapply(strsplit(as.character(counts$`Gene|Length`),".",fixed = TRUE),function(x) x[1]))
      counts <- aggregate(counts[,-1], list(Gene=counts[,1]), FUN = sum)
      rownames(counts) <- counts$Gene
      counts <- counts[,-1]
      
      pathway_gene <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/gene/",pathway,".csv"),header = T)
    }
    #get pathway genes
    {
      j=1
      pathway_gene_count={}
      while(j<=nrow(pathway_gene)){
        target <- pathway_gene[j,which(colnames(pathway_gene)=="ensembl")]
        gene_symbol <- pathway_gene[j,which(colnames(pathway_gene)=="Gene")]
        #target <- pathway_gene[j,which(colnames(pathway_gene)=="Gene.ID")]
        #gene_symbol <- pathway_gene[j,which(colnames(pathway_gene)=="Gene.Name")]
        if(length(grep(target,rownames(counts)))==0) {
          print(paste0("No ",target," in this dataset."))
          temp <- as.data.frame(array(,dim=c(1,ncol(counts))))
          temp[1,] <- 0
          rownames(temp) <- gene_symbol
          j=j+1
        } else {
          #temp <- counts[which(rownames(counts)==target),]
          temp <- counts[grep(target,rownames(counts),fixed=TRUE),]  #for ensg
          rownames(temp) <- gene_symbol
          pathway_gene_count <- rbind(pathway_gene_count,temp)
          j=j+1
        }
      }
      pathway_gene_count_log2 <- log2(pathway_gene_count+1)
      pathway_gene_count_log2_colmean <- colMeans(pathway_gene_count_log2)
      pathway_gene_count_log2_colmean <- as.data.frame(pathway_gene_count_log2_colmean)
      colnames(pathway_gene_count_log2_colmean) <- pathway
      
      pathway_gene_count_log2_colsum <- colSums(pathway_gene_count_log2)
      pathway_gene_count_log2_colsum <- as.data.frame(pathway_gene_count_log2_colsum)
      colnames(pathway_gene_count_log2_colsum) <- pathway
      
      pathway_gene_count_colmean <- colMeans(pathway_gene_count)
      pathway_gene_count_colmean <- as.data.frame(pathway_gene_count_colmean)
      colnames(pathway_gene_count_colmean) <- pathway
      
      dir.create(path = paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i]), recursive = TRUE)
      write.table(pathway_gene_count,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_matrix.txt"),quote = FALSE,sep = "\t")
      write.table(pathway_gene_count_log2_colmean,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_log2.txt"),quote = FALSE,sep = "\t")
      write.table(pathway_gene_count_log2_colsum,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_log2sum.txt"),quote = FALSE,sep = "\t")
      write.table(pathway_gene_count_colmean,paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,".txt"),quote = FALSE,sep = "\t")
    }
    #plot only compares group
    {
      prefix <- "GSE174302"
      des <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/group_info/",prefix[i],".csv"),header = TRUE, row.names = 1)
      log2_mean <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_log2sum.txt"),sep = "\t",header = TRUE, row.names = 1,check.names=FALSE)
      rownames(log2_mean) <- gsub(".","-",fixed= TRUE,rownames(log2_mean))
      des$sample_id <- rownames(des)
      des <- des[-which(des$group=="PBMC"),]
      #des <- des[-which(rownames(des)=="CRC-PKU-29-PBMC"),]
      log2_mean$sample_id <- rownames(log2_mean)
      boxplot <- left_join(des,log2_mean,by = c('sample_id'='sample_id'))
      
      #boxplot <- na.omit(boxplot)
      boxplot$group <- factor(boxplot$group,levels=c("negative","positive"))
      #boxplot$group <- factor(boxplot$group,levels=c("HD","CRC","STAD"))
      #my_comparisons <- list(c("STAD","HD"),c("CRC","HD"))
      #my_comparisons <- list(c("CRC","HD"))
      my_comparisons <- list(c("negative","positive"))
      p <- ggplot(boxplot,aes(x=group,y=boxplot[[pathway]],fill=group))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
        geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
        scale_fill_manual(values = c("CRC"="#FCB514","positive"="red","negative"="blue")) +
        #scale_fill_manual(values = c("CRC"="#FCB514","STAD"="red","HD"="blue")) +
        #scale_fill_manual(values = c("#EE7621","blue"))+
        ylab(pathway)+
        theme_bw()+
        theme(#legend.position="right",
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(size=1, colour = "black"),
          legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
          plot.title = element_text(hjust = 0.5,size=30,face="bold"),
          axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
          axis.text.y = element_text(face="bold",  color="black", size=24),
          axis.title.x = element_text(face="bold", color="black", size=24),
          axis.title.y = element_text(face="bold",color="black", size=24))
      
      dir.create(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"), recursive = TRUE)
      m=mean(boxplot[boxplot$group=="negative",ncol(boxplot)])-mean(boxplot[boxplot$group=="positive",ncol(boxplot)])
      if(m>0){
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  size = 8,
                                  vjust = 0.5,
                                  method.args = list(alternative = "greater"),
                                  label = "p.signif"
        )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.greater"), face="bold",fill="Type")
        ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_greater.pdf"),width = 6,height = 8, units = "in", dpi = 300)
      } else if(m==0){
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  size = 8,
                                  vjust = 0.5,
                                  method.args = list(alternative = "two.sided"),
                                  label = "p.signif"
        )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.twosided"), face="bold",fill="Type")
        ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_twosided.pdf"),width = 6,height = 8, units = "in", dpi = 300)
      } else {
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  size = 8,
                                  vjust = 0.5,
                                  method.args = list(alternative = "less"),
                                  label = "p.signif"
        )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.less"), face="bold",fill="Type")
        ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_less.pdf"),width = 6,height = 8, units = "in", dpi = 300)
      }
    }
    
    #plot compares stage in cancer for GSE174302
    {
      des <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/group_info/GSE174302_CRCSTAD_vs_NC_stage.csv"),header = TRUE, row.names = 1)
      log2_mean <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_matrix/",prefix[i],"/",prefix[i],"_",pathway,"_log2sum.txt"),sep = "\t",header = TRUE, row.names = 1,check.names=FALSE)
      rownames(log2_mean) <- gsub(".","-",fixed= TRUE,rownames(log2_mean))
      des$sample_id <- rownames(des)
      log2_mean$sample_id <- rownames(log2_mean)
      boxplot <- left_join(des,log2_mean,by = c('sample_id'='sample_id'))
      
      #boxplot <- na.omit(boxplot)
      boxplot$group <- factor(boxplot$group,levels=c("Cancer","HD"))
      #my_comparisons <- list(c("STAD","HD"),c("CRC","HD"))
      boxplot$Stage.simplified <- factor(boxplot$Stage.simplified,levels=c("No stage","Stage I","Stage II","Stage III","Stage IV"))
      my_comparisons <- list(c("Stage I","No stage"),c("Stage II","No stage"),c("Stage III","No stage"),c("Stage IV","No stage"))
      p <- ggplot(boxplot,aes(x=Stage.simplified,y=boxplot[[pathway]],fill=group))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
        geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
        scale_fill_manual(values = c("Cancer"="red","HD"="blue")) +
        #scale_fill_manual(values = c("CRC"="#FCB514","STAD"="red","HD"="blue")) +
        #scale_fill_manual(values = c("#EE7621","blue"))+
        ylab(pathway)+
        theme_bw()+
        theme(#legend.position="right",
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(size=1, colour = "black"),
          legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
          plot.title = element_text(hjust = 0.5,size=30,face="bold"),
          axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
          axis.text.y = element_text(face="bold",  color="black", size=24),
          axis.title.x = element_text(face="bold", color="black", size=24),
          axis.title.y = element_text(face="bold",color="black", size=24))
      
      dir.create(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"), recursive = TRUE)
      m=mean(boxplot[boxplot$group=="Cancer",ncol(boxplot)])-mean(boxplot[boxplot$group=="HD",ncol(boxplot)])
      if(m>0){
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  size = 8,
                                  vjust = 0.5,
                                  method.args = list(alternative = "greater"),
                                  label = "p.signif"
        )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.greater"), face="bold",fill="Type")
        ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_greater.pdf"),width = 6,height = 8, units = "in", dpi = 300)
      } else if(m==0){
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  size = 8,
                                  vjust = 0.5,
                                  method.args = list(alternative = "two.sided"),
                                  label = "p.signif"
        )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.twosided"), face="bold",fill="Type")
        ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_twosided.pdf"),width = 6,height = 8, units = "in", dpi = 300)
      } else {
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  size = 8,
                                  vjust = 0.5,
                                  method.args = list(alternative = "less"),
                                  label = "p.signif"
        )+labs(x="",y="log2(TPM+1)",title=paste0(prefix[i],"\nwilcox.test.less"), face="bold",fill="Type")
        ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_less.pdf"),width = 9,height = 8, units = "in", dpi = 300)
      }
    }
    
    #plot GSE27562 (PBMC RNA array)
    {
      pathway <- "Inhibitory immune genes"
      prefix <- "GSE27562"
      pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/GSE27562_T cell receptor_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
      pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/GSE27562_B cell receptor_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
      pathway_count <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/GSE27562_Inhibotory immune genes_TPM_pathway_level_log2.csv",header = TRUE , row.names = 1)
      pathway_count <- as.data.frame(t(pathway_count))
      pathway_count$ID <- rownames(pathway_count)
      
      group <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/01.cancer microbiome/05.selectbywilcox/Pathway_level/signature_mean_log2_TPM+1/GSE27562/GSE27562_des.csv",header = TRUE)
      boxplot <- left_join(pathway_count,group, by=c("ID"="ID"))
      
      my_comparisons <- list(c("Cancer","Healthy"))
      boxplot$group <- factor(boxplot$group,levels = c("Healthy","Cancer"))
      p <- ggplot(boxplot,aes(x=group,y=boxplot[[pathway]],fill=group))+geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
        geom_point(size = 1, position = position_jitterdodge(dodge.width=0.3,jitter.width = 0.3))+
        scale_fill_manual(values = c("Cancer"="red","Healthy"="blue")) +
        #scale_fill_manual(values = c("CRC"="#FCB514","STAD"="red","HD"="blue")) +
        #scale_fill_manual(values = c("#EE7621","blue"))+
        ylab(pathway)+
        theme_bw()+
        theme(#legend.position="right",
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(size=1, colour = "black"),
          legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
          plot.title = element_text(hjust = 0.5,size=30,face="bold"),
          axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
          axis.text.y = element_text(face="bold",  color="black", size=24),
          axis.title.x = element_text(face="bold", color="black", size=24),
          axis.title.y = element_text(face="bold",color="black", size=24))
      
      dir.create(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"), recursive = TRUE)
      m=mean(boxplot[boxplot$group=="Cancer",1])-mean(boxplot[boxplot$group=="Healthy",1])
      if(m>0){
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  size = 8,
                                  vjust = 0.5,
                                  method.args = list(alternative = "greater"),
                                  label = "p.signif"
        )+labs(x="",y="log2(Signal Intensity+1)",title=paste0(prefix[i],"\nwilcox.test.greater"), face="bold",fill="Type")
        ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_greater.pdf"),width = 6,height = 8, units = "in", dpi = 300)
      } else if(m==0){
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  size = 8,
                                  vjust = 0.5,
                                  method.args = list(alternative = "two.sided"),
                                  label = "p.signif"
        )+labs(x="",y="log2(Signal Intensity+1)",title=paste0(prefix[i],"\nwilcox.test.twosided"), face="bold",fill="Type")
        ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_twosided.pdf"),width = 6,height = 8, units = "in", dpi = 300)
      } else {
        p <- p+stat_compare_means(comparisons = my_comparisons,
                                  method = "wilcox.test",
                                  size = 8,
                                  vjust = 0.5,
                                  method.args = list(alternative = "less"),
                                  label = "p.signif"
        )+labs(x="",y="log2(Signal Intensity+1)",title=paste0(prefix[i],"\nwilcox.test.less"), face="bold",fill="Type")
        ggsave(p,device = "pdf",path=paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/",pathway,"_boxplot/"),filename = paste0(prefix[i],"_less.pdf"),width = 9,height = 8, units = "in", dpi = 300)
      }
    }
  }
  
  #immune gene boxplot
  {
    #test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/14.differential analysis for each alterations/Expression/matrix/20211113_TPM.txt",sep = "\t", header = TRUE, row.names = 1)
    #plot <- as.data.frame(t(test[which(rownames(test)=="lncRNA-GC1|2145"),]))
    #test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20210722_multiomics/Altpromoter/Altpromoter_ML.txt",sep = "\t", header = TRUE, row.names = 1)
    
    #test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/pico_PBMC_featurecount_intron_spanning_noMTRNA_TPM.txt",sep = "\t", header = TRUE, row.names = 1)  
    #test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/pico_tissue_featurecount_intron_spanning_noMTRNA_TPM.txt",sep = "\t", header = TRUE, row.names = 1)  
    test <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Plasma_TPM.txt",sep = "\t", header = TRUE, row.names = 1)
    
    #wangtantan HPK1, MAP4K1, ENSG00000104814
    #mTOR: ENSG00000198793
    
    #Upregulated inhibitory receptor
    #ENSG00000120217: CD274(PD-L1)
    #ENSG00000103855: CD276(B7-H3)
    #ENSG00000134258: VTCN1(B7-H4)
    #ENSG00000114455: HHLA2(B7-H5)
    #ENSG00000188389: PDCD1(PD1)
    #ENSG00000163599: CTLA4
    #ENSG00000089692: LAG3
    #ENSG00000135077: TIM3
    #ENSG00000079385: CEACAM1(TIM3 ligand)
    #ENSG00000181847: TIGIT
    
    #Downregulated immune activation genes
    #ENSG00000153563: CD8A
    #ENSG00000172116: CD8B
    #ENSG00000010610: CD4
    #ENSG00000115085: ZAP70
    #ENSG00000113263: ITK
    #ENSG00000111537: IFNG
    #ENSG00000074966: RLK
    #ENSG00000170345: FOS
    
    #CD8 T cell marker
    #ENSG00000211789: TRAV12-2
    #ENSG00000153563: CD8A
    #ENSG00000172116: CD8B
    
    
    
    #CD4 memory naive T cell marker
    #"ENSG00000231160"="KLF3-AS1"
    #"ENSG00000176293"="ZNF135"
    #"ENSG00000125864"="BFSP1"
    #"ENSG00000166573"="GALR1"
    #"ENSG00000168067"="MAP4K2"
    #"ENSG00000134765"="DSC1"
    #"ENSG00000204789"="ZNF204P"
    #"ENSG00000146904"="EPHA1"
    #"ENSG00000090554"="FLT3LG"
    #"ENSG00000138795"="LEF1"
    #"ENSG00000197093"="GAL3ST4"
    #"ENSG00000154764"="WNT7A"
    #"ENSG00000129158"="SERGEF"
    #"ENSG00000164512"="ANKRD55"
    #"ENSG00000142102"="PGGHG"
    #"ENSG00000083812"="ZNF324"
    
    #CD 4 memory resting cell marker
    #"ENSG00000163346"="RCAN3"
    #"ENSG00000152518"="ZFP36L2"
    #"ENSG00000159023"="EPB41"
    #"ENSG00000134954"="ETS1"
    #"ENSG00000135722"="FBXL8"
    #"ENSG00000165496"="RPL10L"
    #"ENSG00000197635"="DPP4"
    #"ENSG00000117602"="RCAN3"
    #"ENSG00000139193"="CD27"
    #"ENSG00000170128"="GPR25"
    #"ENSG00000177272"="KCNA3"
    #"ENSG00000100100"="PIK3IP1"
    #"ENSG00000198851"="CD3E"
    #"ENSG00000107742"="SPOCK2"
    
    #gama delta
    "ENSG00000109943"="CRTAM"
    "ENSG00000158050"="DUSP2"
    "ENSG00000113088"="GZMK"
    "ENSG00000139187"="KLRG1"
    "ENSG00000147231"="RADX"
    "ENSG00000233402"="TARDBPP1"
    "ENSG00000160791"="CCR5"
    "ENSG00000198342"="ZNF442"
    "ENSG00000211829"="TRDC"
    "ENSG00000211804"="TRDV1"
    "ENSG00000211804"="TRDV2"
    
    List <-  c("ENSG00000120217"="CD274 (PD-L1)",
               "ENSG00000089692"="LAG3",
               "ENSG00000109943"="CRTAM",
               "ENSG00000158050"="DUSP2",
               "ENSG00000113088"="GZMK",
               "ENSG00000139187"="KLRG1",
               "ENSG00000147231"="RADX",
               "ENSG00000233402"="TARDBPP1",
               "ENSG00000160791"="CCR5",
               "ENSG00000198342"="ZNF442",
               "ENSG00000211829"="TRDC",
               "ENSG00000211804"="TRDV1",
               "ENSG00000211804"="TRDV2",
               "ENSG00000211789"="TRAV12-2",
               "ENSG00000153563"="CD8A",
               "ENSG00000172116"="CD8B",
               "ENSG00000010610"="CD4",
               "ENSG00000178562"="CD28",
               "ENSG00000177455"="CD19",
               "ENSG00000115085"="ZAP70",
               "ENSG00000113263"="ITK",
               "ENSG00000111537"="IFNG",
               "ENSG00000074966"="RLK",
               "ENSG00000170345"="FOS",
               "ENSG00000103855"="CD276(B7-H3)",
               "ENSG00000134258"="VTCN1(B7-H4)",
               "ENSG00000114455"="HHLA2(B7-H5)",
               "ENSG00000188389"="PDCD1(PD1)",
               "ENSG00000163599"="CTLA4",
               "ENSG00000089692"="LAG3",
               "ENSG00000135077"="TIM3",
               "ENSG00000079385"="CEACAM1(TIM3 ligand)",
               "ENSG00000181847"="TIGIT",
               "ENSG00000198793"="mTOR",
               "ENSG00000213809"="NKG2D",
               "ENSG00000134545"="NKG2A")
    
    #List <- as.list(rownames(LM22))
    #names(List) <- rownames(LM22)
    
    i=1
    while(i <= length(List)){
      Gene_ID <- names(List)[i]
      Gene_name <- as.character(List[i])
      {
        plot <- as.data.frame(t(test[grep(Gene_ID,rownames(test)),]))
        plot$sample <- rownames(plot)
        plot <- plot[grep("pico",rownames(plot)),]
        plot <- plot[-grep("NC.PKU.mix..pico|CRC.PKU.mix1.pico|CRC.PKU.5.pico|NC.PKU.mix17.pico|STAD.PKU.4.pico",rownames(plot)),]
        
        plot$group <- as.character(lapply(strsplit(rownames(plot),".",fixed=TRUE),function(x) x[1]))
        #plot$group <- as.character(lapply(strsplit(rownames(plot),".",fixed=TRUE),function(x) tail(x,n=1)))
        
        plot$group <- gsub("NC","HD",plot$group)
        
        my_comparisons <- list(c("HD","STAD"),c("HD","CRC"))
        #my_comparisons <- list(c("CRC","HD"))
        #my_comparisons <- list(c("T","N"))
        plot$group <- factor(plot$group,levels=c("HD","STAD","CRC"))
        
        #plot$group <- gsub("CRC","GIC",plot$group)
        #plot$group <- gsub("STAD","GIC",plot$group)
        #my_comparisons <- list(c("HD","GIC"))
        #plot$group <- factor(plot$group,levels=c("HD","GIC"))
        forplot <- plot[,c(1,ncol(plot))]
        colnames(forplot) <- c("value","group")
        #p <- ggplot(forplot[-which(rownames(forplot)=="CRC.PKU.29.pico"),],aes(x=group,y=value,fill = group))+
        p <- ggplot(forplot,aes(x=group,y=value,fill = group))+
          geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
          geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1))+
          scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue"))+
          #scale_fill_manual(values=c("#87CEEB","#FCB514"))+
          #scale_fill_brewer(palette="Blues") +
          #ylim(0,25)+
          theme_bw()+
          xlab("")+
          ylab(Gene_name)+
          #ylab(colnames(plot)[1])+
          theme(#legend.position="right",
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            axis.line = element_line(size=1, colour = "black"),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
            plot.title = element_text(hjust = 0.5,size=36,face="bold"),
            axis.text.x = element_blank(),
            #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
            axis.text.y = element_text(face="bold",  color="black", size=24),
            axis.title.x = element_text(face="bold", color="black", size=24),
            axis.title.y = element_text(face="bold",color="black", size=24))+
          stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test",
                             method.args = list(alternative = "two.sided",paired = TRUE),
                             label = "p.signif",
                             size = 10,
                             vjust = 0.5)
        ggsave(plot = p,path = "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 3/Marker", filename = paste0(List[i],".pdf"),device = "pdf",height = 5.0,width = 3.8)
        i=i+1
      }
    }
  }
}

#Figure 4
{
#Machine learning result comparison
{
  #ROC
  RNA_color <- c("cfDNA methylation (Promoter)"="#666666",
                 "cfDNA methylation (CpG island)"="#666666",
                 
                 "cfRNA alternative promoter usage"="#5CACEE",
                 "cfRNA abundance"="#5CACEE",
                 "cfRNA allele specific expression"="#5CACEE",
                 "cfRNA editing"="#5CACEE",
                 "Chimeric cfRNA" ="#5CACEE",
                 "Microbe genus abundance\n (microbial cfRNA sequence)"="#5CACEE",
                 "cf-miRNA"="#CD0000",
                 "cfRNA SNV"="#5CACEE",
                 "Transposable elements abudance\n(cfRNA expression)"="#5CACEE",
                 "cfRNA splicing"="#5CACEE",
                 
                 "cfDNA copy number"="#EEB4B4",
                 "cfDNA nucleosome occupancy"="#EEB4B4",
                 "cfDNA end motif usage"="#EEB4B4",
                 "cfDNA fragment size"="#EEB4B4",
                 "WPS-divide-bgCOV"="#EEB4B4",
                 "WPS-divide-COV"="#EEB4B4",
                 "WPS-divide-bgWPS"="#EEB4B4",
                 "WPS-substract-bgWPS"="#EEB4B4",
                 "WPS-subtractWPS-divideCOV"="#EEB4B4",
                 "cfDNA window protection score"="#EEB4B4",
                 
                 "Integrated"="#FFFFFF",
                 "DNA+RNA alterations combined"="#000000",
                 "DNA alterations combined"=alpha("#EEB4B4",alpha=1),
                 "RNA alterations combined"=alpha("#5CACEE",alpha=1)
                 #"DNA+RNA alterations combined"="#FFFFFF",
                 #"DNA alterations combined"="#FFFFFF",
                 #"RNA alterations combined"="#FFFFFF"
  )

  timepoint <- date()
  timepoint <- gsub(":","_",timepoint)
  timepoint <- gsub(" ","_",timepoint)
  
  #other metric plot (include AUC, just change metric)
  {
    metric <- "AUC"
    library(dplyr)
    library(ggrepel)
    #,"STAD_vs_CRC"
    dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220718_yuhuan/All/late_integration_result_300f_50b_30bsn"
    for(compare in c("CRC_vs_HD","STAD_vs_HD")){
      Early_AUC_integrated <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220718_yuhuan/early/late_integration_result_300f_50b_30bsn/",compare,"/1_late_integrated_AUC.txt"),sep = "\t",header = TRUE,check.names = FALSE)
      Early_AUC_integrated <- Early_AUC_integrated[,grep(metric,colnames(Early_AUC_integrated)),drop = FALSE]
      colnames(Early_AUC_integrated) <- c("Integrated")
      Early_AUC <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220718_yuhuan/early/late_integration_result_300f_50b_30bsn/",compare,"/1_test_",metric,".txt"),sep = "\t",header = TRUE,check.names = FALSE)
      Early_AUC <- Early_AUC[,-grep("seed",colnames(Early_AUC))]
      #Early_AUC <- Early_AUC[,-which(colnames(Early_AUC)=="comparison")]
      Early_AUC <- cbind(Early_AUC,Early_AUC_integrated)
      plot_Early_AUC <- data.frame("Alteration"=colnames(Early_AUC),"mean"=colMeans(Early_AUC),"sd"=colStdevs(Early_AUC))
      
      AUC_integrated <- read.csv(paste0(dir,"/",compare,"/1_late_integrated_AUC.txt"),sep = "\t",header = TRUE,check.names = FALSE)
      AUC_integrated <- AUC_integrated[,grep(metric,colnames(AUC_integrated)),drop = FALSE]
      colnames(AUC_integrated) <- c("Integrated")
      AUC_plot <- read.csv(paste0(dir,"/",compare,"/1_test_",metric,".txt"),sep = "\t",header = TRUE,check.names = FALSE)
      Comparison <- compare
      AUC_plot <- AUC_plot[,-grep("seed",colnames(AUC_plot))]
      #AUC_plot <- AUC_plot[,-which(colnames(AUC_plot)=="comparison")]
      AUC_plot <- cbind(AUC_plot,AUC_integrated)
      plot_AUC <- data.frame("Alteration"=colnames(AUC_plot),"mean"=colMeans(AUC_plot),"sd"=colStdevs(AUC_plot))
      
      k=1
      wilcox_early_all <- data.frame("Alteration"=colnames(AUC_plot),"pvalue"=NA)
      while(k<=ncol(AUC_plot)){
        test <- wilcox.test(AUC_plot[,k],Early_AUC[,k])
        wilcox_early_all[k,"pvalue"]<- test$p.value
        k=k+1
      }
      write.csv(wilcox_early_all,paste0(dir,"/",timepoint,"/",metric,"_",Comparison,".csv"))
      
      plot_AUC <- left_join(plot_AUC, plot_Early_AUC,by = c("Alteration"="Alteration"))
      #plot_AUC$group <- "All samples"
      #plot_Early_AUC$group <- "Early stage samples"
      #plot_AUC <- rbind(plot_AUC,plot_Early_AUC) 
      
      {
        plot_AUC$Alteration <- gsub(paste0(metric,"-CNV_ML"),"cfDNA copy number",plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-Expression_ML"),"cfRNA abundance",plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-AlternativePromoter_ML"),"cfRNA alternative promoter usage",plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-ASE_ML"),"cfRNA allele specific expression",plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-chimeric_ML"),"Chimeric cfRNA",plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-Editing_ML"),"cfRNA editing",plot_AUC$Alteration)
        #plot_AUC$Alteration <- gsub(paste0(metric,"-Nucleosome_ML"),"cfDNA nucleosome occupancy",plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-Mutation_ML"),"cfRNA SNV",plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-Splicing_ML"),"cfRNA splicing",plot_AUC$Alteration)
        #plot_AUC$Alteration <- gsub(paste0(metric,"-Microbe_ML"),"Microbe genus abundance\n (microbial cfRNA sequence)",plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-EndMotif_ML"),"cfDNA end motif usage",plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-FragSize_ML"),"cfDNA fragment size",plot_AUC$Alteration)
        #plot_AUC$Alteration <- gsub(paste0(metric,"-TE_ML"),"Transposable elements abudance\n(cfRNA expression)",plot_AUC$Alteration)
        #AUC_plot$data_type <- gsub("microcf","cf-miRNA",AUC_plot$data_type)
        plot_AUC$Alteration <- gsub(paste0(metric,"-Medip_ML"),"cfDNA methylation (Promoter)",plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-Medip_bin_ML"),"cfDNA methylation (CpG island)",fixed = TRUE,plot_AUC$Alteration)
        plot_AUC$Alteration <- gsub(paste0(metric,"-WPS-divide-bgCOV_matrix_150TSS50"),"cfDNA window protection score",fixed = TRUE,plot_AUC$Alteration)
        #plot_AUC$Alteration <- gsub("spec_youden","Integrated",fixed = TRUE,plot_AUC$Alteration)
        
      }
      
      #ROC_barplot
      {
        rank <- plot_AUC[order(plot_AUC[1:length(unique(plot_AUC$Alteration)),]$mean.x,decreasing = TRUE),]$Alteration
        rank <- grep("Integrated",rank,invert = TRUE, value = TRUE)
        plot_AUC$Alteration <- factor(plot_AUC$Alteration,levels = rev(c("Integrated",rank)))
        p<-ggplot(plot_AUC,aes(x=Alteration,y=mean.x,fill=Alteration))+
          #facet_wrap(~metric,nrow = 4)+
          geom_bar(stat = "identity",colour = "black",width = 0.4,position = position_nudge(x = 0.2))+
          geom_bar_pattern(aes(x=Alteration,y=mean.y,fill = Alteration),stat = "identity",colour = "black",width = 0.4,position = position_nudge(x = -0.2))+
          geom_text(aes(x=Alteration,y=0.1,label=round(mean.x,digits = 3)),size = 4,angle = 90,position = position_nudge(x = 0.2))+
          geom_text(aes(x=Alteration,y=0.1,label=round(mean.y,digits = 3)),size = 4,angle = 90,position = position_nudge(x = -0.2))+
          geom_errorbar(aes(x=Alteration, ymin=mean.x-sd.x, ymax=mean.x+sd.x), width=0.4, colour="black", alpha=0.5, size=0.5,position = position_nudge(x = 0.2))+
          geom_errorbar(aes(x=Alteration, ymin=mean.y-sd.y, ymax=mean.y+sd.y), width=0.4, colour="black", alpha=0.5, size=0.5,position = position_nudge(x = -0.2))+
          xlab("")+
          labs(title = Comparison)+
          scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("0","0.2","0.4","0.6","0.8","1.0"),expand = c(0,0),limits = c(0,1))+
          scale_fill_manual(values = RNA_color)+
          #scale_fill_manual(values = rev(blues20))+
          #geom_vline(xintercept=1.5,linetype = 2)+
          #geom_vline(xintercept=3.5,linetype = 2)+
          xlab("")+
          ylab(metric)+
          theme_bw()+
          theme(
            plot.margin = unit(x=c(15,5,-10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(face="bold", color="black", size=18),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        dir.create(paste0(dir,"/",timepoint),recursive = TRUE)
        ggsave(paste0(dir,"/",timepoint,"/",metric,"_",Comparison,".pdf"),p,width = 12.26,height = 6.57,device = "pdf")
      }
    }
  }
  
  #Alteration contribution plot
  {
    dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220706_yuhuan/DNA+RNA/late_integration_result_300f_50b_50bsn"
    for(compare in c("CRC_vs_HD","STAD_vs_HD","STAD_vs_CRC")){
      Coefficient_plot <- read.csv(paste0(dir,"/",compare,"/1_coefficients_for_different_alterations.txt"),sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
      Comparison <- compare
      
      {
        rownames(Coefficient_plot) <- gsub("CNV_ML.txt","cfDNA copy number",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("Expression_ML.txt","Human gene abundance\n(cfRNA expression)",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("AlternativePromoter_ML.txt","cfRNA alternative promoter usage",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("ASE_ML.txt","cfRNA allele specific expression",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("chimeric_ML.txt","Chimeric cfRNA",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("Editing_ML.txt","cfRNA editing",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("Nucleosome_ML.txt","cfDNA nucleosome occupancy",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("Mutation_ML.txt","cfRNA SNV",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("Splicing_ML.txt","cfRNA splicing",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("Microbe_ML.txt","Microbe genus abundance\n (microbial cfRNA sequence)",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("EndMotif_ML.txt","cfDNA end motif usage",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("FragSize_ML.txt","cfDNA fragment size",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("TE_ML.txt","Transposable elements abudance\n(cfRNA expression)",rownames(Coefficient_plot))
        #AUC_plot$data_type <- gsub("microcf","cf-miRNA",AUC_plot$data_type)
        rownames(Coefficient_plot) <- gsub("Medip_ML.txt","cfDNA methylation (Promoter)",rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("Medip_bin_ML.txt","cfDNA methylation (CpG island)",fixed = TRUE,rownames(Coefficient_plot))
        rownames(Coefficient_plot) <- gsub("WPS-divide-bgCOV_matrix_150TSS50.txt","cfDNA window protection score",fixed = TRUE,rownames(Coefficient_plot))
        Coefficient_plot <- Coefficient_plot/colSums(Coefficient_plot)
        plot_Coefficient <- data.frame("Alteration"=rownames(Coefficient_plot),"mean"=rowMeans(Coefficient_plot),"sd"=apply(Coefficient_plot,1, sd))
      }
      
      #plot alterations contribution
      {
        rank <- plot_Coefficient[order(plot_Coefficient$mean,decreasing = TRUE),]$Alteration
        plot_Coefficient$Alteration <- factor(plot_Coefficient$Alteration,levels = rev(rank))
        p<-ggplot(plot_Coefficient,aes(x=Alteration,y=mean,fill=Alteration))+
          #facet_wrap(~metric,nrow = 4)+
          geom_bar(stat = "identity",colour = "black")+
          geom_text(aes(x=Alteration,y=mean+0.05,label=round(mean,digits = 3)),size = 4,angle = 0)+
          geom_errorbar(aes(x=Alteration, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
          xlab("")+
          labs(title = Comparison)+
          coord_flip()+
          #scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("0","0.2","0.4","0.6","0.8","1.0"),expand = c(0,0),limits = c(0,1))+
          scale_fill_manual(values = RNA_color)+
          #scale_fill_manual(values = rev(blues20))+
          #geom_vline(xintercept=1.5,linetype = 2)+
          #geom_vline(xintercept=3.5,linetype = 2)+
          xlab("")+
          ylab("Contribution")+
          theme_bw()+
          theme(
            plot.margin = unit(x=c(15,5,10,70),units="pt"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(face="bold", color="black", size=18),
            #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
            axis.ticks.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.x = element_text(face="bold", color="black", size=12, angle = 45,hjust = 1,vjust = 1),
            axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))
        #dir.create(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220706_yuhuan/",Comparison),recursive = TRUE)
        ggsave(paste0(dir,"/Coefficient_",Comparison,".pdf"),p,width = 10,height = 8,device = "pdf")
      }
    }
  }
  
  #scatter plot
  {
    Group_color <- c(CRC=alpha("#FCB514",alpha = 1),STAD=alpha("red",alpha = 1),HD=alpha("blue",alpha = 1),GIC="#FCB514")
    for_scatter <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220706_yuhuan/DNA+RNA/late_integration_result_300f_50b_50bsn/STAD_vs_CRC/1_test_prob.txt",sep = "\t", header = TRUE, check.names = FALSE)
    for_scatter_plot <- aggregate(for_scatter[,-grep("seed|label|samples",colnames(for_scatter))], list(samples = for_scatter$samples),FUN = mean)
    for_scatter_plot$label <- as.character(lapply(strsplit(as.character(for_scatter_plot$samples),"-",fixed = TRUE), function(x) x[1]))
    for_scatter_plot$label <- gsub("NC","HD",for_scatter_plot$label)
    p <- ggplot(for_scatter_plot,aes(x=`STAD-Microbe_ML`,y=`STAD-CNV_ML`,color=label))+geom_point()+
      scale_color_manual(values = Group_color)+
      theme_bw()
    ggsave(filename="/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20220706_yuhuan/DNA+RNA/late_integration_result_300f_50b_50bsn/STAD_vs_CRC_Microbe_EndMotif.pdf",p,width = 5.99, height = 4.81,device = "pdf")
  }
}

#complementary of different alterations
{
  samples <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/ML_samples.csv",row.names = NULL)
  positive <- "STAD|CRC"
  negative <- "HD"
  positive_samples <- as.character(samples[grep(positive,samples$Group),]$sample_id)
  n_positive <- length(positive_samples)
  negative_samples <- as.character(samples[grep(negative,samples$Group),]$sample_id)
  n_negative <- length(negative_samples)
  
  library(scatterplot3d)
  
  data1 <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/RNA/Plasma_TPM.txt",sep = "\t",row.names = 1,check.names = FALSE,header = TRUE)
  data2 <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Alterations/matrix_raw/DNA/CPM-TMM_matrix_gene.correctGC.txt",sep = "\t",row.names = 1,check.names = FALSE,header = TRUE)
  colnames(data1) <- gsub("-pico","",colnames(data1))
  colnames(data2) <- gsub("-wgs","",colnames(data2))
  
  grep("ENSG00000141510",rownames(data1),value = TRUE)
  grep("ENSG00000141510",rownames(data2),value = TRUE)
  
  
  
  Group_color <- c(CRC=alpha("#FCB514",alpha = 1),STAD=alpha("red",alpha = 1),HD=alpha("blue",alpha = 1),GIC="#FCB514")
  event_plot <- t(rbind(data1[grep("ENSG00000141510",rownames(data1)),c(positive_samples,negative_samples)],data2[grep("ENSG00000141510",rownames(data2)),c(positive_samples,negative_samples)]))
  event_plot <- as.data.frame(event_plot)
  event_plot$group <- as.character(lapply(strsplit(rownames(event_plot),"-",fixed = TRUE), function(x) x[1]))
  event_plot$group <- gsub("NC","HD",event_plot$group)
  
  xvalues <- event_plot[event_plot$group=="HD",]$`ENSG00000141510.16|4576`
  yvalues <- event_plot[event_plot$group=="HD",]$`ENSG00000141510.16|25772`
  p <- ggplot(event_plot,aes(x=log10(event_plot$`ENSG00000141510.16|4576`),y=log10(event_plot$`ENSG00000141510.16|25772`)))+
    scale_color_manual(values = Group_color)+
    xlab("TP53 cfRNA abundance")+
    ylab("TP53 cfDNA copy number")+
    scale_x_reverse()+
    scale_y_reverse()+
    geom_vline(xintercept = log10(min(xvalues[xvalues!=min(xvalues)])), linetype='dashed')+
    geom_hline(yintercept = log10(min(yvalues[yvalues!=min(yvalues)])), linetype='dashed')+
    geom_point(aes(color = event_plot$group))+
    theme_bw()+
    theme(
      plot.margin = unit(x=c(5,5,5,5),units="pt"),
      legend.position="none",
      panel.grid=element_blank(),
      panel.border=element_rect(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face="bold", color="black", size=18),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
      axis.text.x = element_text(face="bold",color="black", size=12, angle = 0,hjust = 0.5,vjust = 1),
      axis.text.y = element_text(face="bold",color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_text(color="black", size=16),
      axis.title.y = element_text(color="black", size=16))
  ggsave(p,filename = "./TP53_example.pdf",device = "pdf",width = 3.74,height = 3.41)
 
  combined_3d <- as.data.frame(t(rbind(
    data1[event1,as.character(samples$sample_id)],
    data2[event2,as.character(samples$sample_id)],
    data3[event3,as.character(samples$sample_id)]
  )))
  combined_3d$group <- as.character(samples$Group)
  colors <- c("#FCB514", "blue", "red")
  colors <- colors[as.numeric(as.factor(combined_3d$group))]
  s3d <- scatterplot3d(combined_3d[,1:3], pch=16, color = colors,angle = 90,#type="h",
                       grid=TRUE, box=TRUE)
  legend("topright", legend = levels(as.factor(combined_3d$group)),
         col =  c("#FCB514", "blue", "red"), pch = 16, bg = "transparent")
  text(s3d$xyz.convert(combined_3d[,1:3]), labels = rownames(combined_3d),
       cex= 0.5, col = "black")
}

#PCAWG plot
{
  #Functional annotation
  {
    #function for enrichment analysis for RNA expression: clusterprofiler
    KEGG_GO<- function(ensembel_genes,output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC){
      
      forenrich <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                              filters = "ensembl_gene_id",
                              values=ensembel_genes, mart= mart,useCache = FALSE)
      if(length(which(forenrich$entrezgene_id==""))==0){
        forenrich <- forenrich
      } else {
        forenrich <- forenrich[-which(forenrich$entrezgene_id==""),]
      }
      forenrich <- forenrich$entrezgene_id
      
      #KEGG
      {
        KEGG_res <- enrichKEGG(
          forenrich,
          organism = "hsa",
          keyType = "kegg",
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          #universe=as.character(background),
          minGSSize = 0,
          maxGSSize = 500,
          qvalueCutoff = 1,
          use_internal_data = FALSE)
        KEGG_output <- KEGG_res@result
        write.table(KEGG_output,output_KEGG,quote = FALSE,sep = "\t")
      }
      #GO_BP
      {
        GO_BP_res <- enrichGO(
          forenrich,
          'org.Hs.eg.db',
          ont = "BP",
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          #universe=as.character(background),
          minGSSize = 0,
          maxGSSize = 500,
          qvalueCutoff = 1)
        GO_BP_output <- GO_BP_res@result
        write.table(GO_BP_output,output_GO_BP,quote = FALSE,sep = "\t")
      }
      #GO_CC
      {
        GO_CC_res <- enrichGO(
          forenrich,
          'org.Hs.eg.db',
          ont = "CC",
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          #universe=as.character(background),
          minGSSize = 0,
          maxGSSize = 500,
          qvalueCutoff = 1)
        GO_CC_output <- GO_CC_res@result
        write.table(GO_CC_output,output_GO_CC,quote = FALSE,sep = "\t")
      }
      #GO_MF
      {
        GO_MF_res <- enrichGO(
          forenrich,
          'org.Hs.eg.db',
          ont = "MF",
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          #universe=as.character(background),
          minGSSize = 0,
          maxGSSize = 500,
          qvalueCutoff = 1)
        GO_MF_output <- GO_MF_res@result
        write.table(GO_MF_output,output_GO_MF,quote = FALSE,sep = "\t")
      }
    }
  }
  #Outlier Finder
  {
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = length(gene_list), clear = FALSE, width= 60)
    
    i=1
    CRC_final <- as.data.frame(matrix(numeric(0),ncol=17))
    colnames(CRC_final) <- c("gene","FDR_up0.05","FDR_up0.01","FDR_down0.05","FDR_down0.01",
                             "CRC_up0.05","CRC_up0.05_ratio","CRC_up0.05_sample",
                             "CRC_up0.01","CRC_up0.01_ratio","CRC_up0.01_sample",
                             "CRC_down0.05","CRC_down0.05_ratio","CRC_down0.05_sample",
                             "CRC_down0.01","CRC_down0.01_ratio","CRC_down0.01_sample")
    STAD_final <- as.data.frame(matrix(numeric(0),ncol=13))
    colnames(STAD_final) <- c("gene","STAD_up0.05","STAD_up0.05_ratio","STAD_up0.05_sample",
                              "STAD_up0.01","STAD_up0.01_ratio","STAD_up0.01_sample",
                              "STAD_down0.05","STAD_down0.05_ratio","STAD_down0.05_sample",
                              "STAD_down0.01","STAD_down0.01_ratio","STAD_down0.01_sample")
    panCancer_final <- as.data.frame(matrix(numeric(0),ncol=13))
    colnames(panCancer_final) <- c("gene","panCancer_up0.05","panCancer_up0.05_ratio","panCancer_up0.05_sample",
                                   "panCancer_up0.01","panCancer_up0.01_ratio","panCancer_up0.01_sample",
                                   "panCancer_down0.05","panCancer_down0.05_ratio","panCancer_down0.05_sample",
                                   "panCancer_down0.01","panCancer_down0.01_ratio","panCancer_down0.01_sample")
    
    all_final <- as.data.frame(matrix(numeric(0),ncol=43))
    colnames(all_final) <- c("gene","FDR_up0.05","FDR_up0.01","FDR_down0.05","FDR_down0.01",
                             "CRC_up0.05","CRC_up0.05_ratio","CRC_up0.05_sample",
                             "CRC_up0.01","CRC_up0.01_ratio","CRC_up0.01_sample",
                             "CRC_down0.05","CRC_down0.05_ratio","CRC_down0.05_sample",
                             "CRC_down0.01","CRC_down0.01_ratio","CRC_down0.01_sample",
                             "gene","STAD_up0.05","STAD_up0.05_ratio","STAD_up0.05_sample",
                             "STAD_up0.01","STAD_up0.01_ratio","STAD_up0.01_sample",
                             "STAD_down0.05","STAD_down0.05_ratio","STAD_down0.05_sample",
                             "STAD_down0.01","STAD_down0.01_ratio","STAD_down0.01_sample",
                             "gene","panCancer_up0.05","panCancer_up0.05_ratio","panCancer_up0.05_sample",
                             "panCancer_up0.01","panCancer_up0.01_ratio","panCancer_up0.01_sample",
                             "panCancer_down0.05","panCancer_down0.05_ratio","panCancer_down0.05_sample",
                             "panCancer_down0.01","panCancer_down0.01_ratio","panCancer_down0.01_sample")
    
    while(i<=length(gene_list)){
      gene <- gene_list[i]
      
      transfer <- transfer_all[which(rownames(transfer_all)==gene),]
      
      NC <- as.data.frame(t(transfer[,grep("NC",colnames(transfer))]))
      NC$sample <- rownames(NC)
      CRC <- as.data.frame(t(transfer[,grep("CRC",colnames(transfer))]))
      CRC$sample <- rownames(CRC)
      STAD <- as.data.frame(t(transfer[,grep("STAD",colnames(transfer))]))
      STAD$sample <- rownames(STAD)
      panCancer <- as.data.frame(t(transfer[,-grep("NC",colnames(transfer))]))
      panCancer$sample <- rownames(panCancer)
      
      NC$sample <- rownames(NC)
      NC_sort <- NC[order(NC[1],decreasing = TRUE),]
      FDR_up0.05 <- NC_sort[ceiling(nrow(NC)*0.05),1]
      FDR_up0.01 <- NC_sort[ceiling(nrow(NC)*0.01),1]
      FDR_down0.05 <- NC_sort[ceiling(nrow(NC)*0.95),1]
      FDR_down0.01 <- NC_sort[ceiling(nrow(NC)*0.99),1] #ceiling is more strict for down
      
      #CRC
      {
        CRC_up0.05 <- nrow(CRC[CRC[1]>FDR_up0.05,])
        CRC_up0.05_sample <- paste(CRC[CRC[1]>FDR_up0.05,]$sample,collapse = "\n")
        CRC_up0.05_ratio <- CRC_up0.05/nrow(CRC)
        
        CRC_up0.01 <- nrow(CRC[CRC[1]>FDR_up0.01,])
        CRC_up0.01_sample <- paste(CRC[CRC[1]>FDR_up0.01,]$sample,collapse = "\n")
        CRC_up0.01_ratio <- CRC_up0.01/nrow(CRC)
        
        CRC_down0.05 <- nrow(CRC[CRC[1]<FDR_down0.05,])
        CRC_down0.05_sample <- paste(CRC[CRC[1]<FDR_down0.05,]$sample,collapse = "\n")
        CRC_down0.05_ratio <- CRC_down0.05/nrow(CRC)
        
        CRC_down0.01 <- nrow(CRC[CRC[1]<FDR_down0.01,])
        CRC_down0.01_sample <- paste(CRC[CRC[1]<FDR_down0.01,]$sample,collapse = "\n")
        CRC_down0.01_ratio <- CRC_down0.01/nrow(CRC)
        
        CRC_tmp <- data.frame(gene,FDR_up0.05,FDR_up0.01,FDR_down0.05,FDR_down0.01,
                              CRC_up0.05,CRC_up0.05_ratio,CRC_up0.05_sample,
                              CRC_up0.01,CRC_up0.01_ratio,CRC_up0.01_sample,
                              CRC_down0.05,CRC_down0.05_ratio,CRC_down0.05_sample,
                              CRC_down0.01,CRC_down0.01_ratio,CRC_down0.01_sample)
      }
      
      #STAD
      {
        STAD_up0.05 <- nrow(STAD[STAD[1]>FDR_up0.05,])
        STAD_up0.05_sample <- paste(STAD[STAD[1]>FDR_up0.05,]$sample,collapse = "\n")
        STAD_up0.05_ratio <- STAD_up0.05/nrow(STAD)
        
        STAD_up0.01 <- nrow(STAD[STAD[1]>FDR_up0.01,])
        STAD_up0.01_sample <- paste(STAD[STAD[1]>FDR_up0.01,]$sample,collapse = "\n")
        STAD_up0.01_ratio <- STAD_up0.01/nrow(STAD)
        
        STAD_down0.05 <- nrow(STAD[STAD[1]<FDR_down0.05,])
        STAD_down0.05_sample <- paste(STAD[STAD[1]<FDR_down0.05,]$sample,collapse = "\n")
        STAD_down0.05_ratio <- STAD_down0.05/nrow(STAD)
        
        STAD_down0.01 <- nrow(STAD[STAD[1]<FDR_down0.01,])
        STAD_down0.01_sample <- paste(STAD[STAD[1]<FDR_down0.01,]$sample,collapse = "\n")
        STAD_down0.01_ratio <- STAD_down0.01/nrow(STAD)
        
        STAD_tmp <- data.frame(gene,STAD_up0.05,STAD_up0.05_ratio,STAD_up0.05_sample,
                               STAD_up0.01,STAD_up0.01_ratio,STAD_up0.01_sample,
                               STAD_down0.05,STAD_down0.05_ratio,STAD_down0.05_sample,
                               STAD_down0.01,STAD_down0.01_ratio,STAD_down0.01_sample)
      }
      
      #panCancer
      {
        panCancer_up0.05 <- nrow(panCancer[panCancer[1]>FDR_up0.05,])
        panCancer_up0.05_sample <- paste(panCancer[panCancer[1]>FDR_up0.05,]$sample,collapse = "\n")
        panCancer_up0.05_ratio <- panCancer_up0.05/nrow(panCancer)
        
        panCancer_up0.01 <- nrow(panCancer[panCancer[1]>FDR_up0.01,])
        panCancer_up0.01_sample <- paste(panCancer[panCancer[1]>FDR_up0.01,]$sample,collapse = "\n")
        panCancer_up0.01_ratio <- panCancer_up0.01/nrow(panCancer)
        
        panCancer_down0.05 <- nrow(panCancer[panCancer[1]<FDR_down0.05,])
        panCancer_down0.05_sample <- paste(panCancer[panCancer[1]<FDR_down0.05,]$sample,collapse = "\n")
        panCancer_down0.05_ratio <- panCancer_down0.05/nrow(panCancer)
        
        panCancer_down0.01 <- nrow(panCancer[panCancer[1]<FDR_down0.01,])
        panCancer_down0.01_sample <- paste(panCancer[panCancer[1]<FDR_down0.01,]$sample,collapse = "\n")
        panCancer_down0.01_ratio <- panCancer_down0.01/nrow(panCancer)
        
        panCancer_tmp <- data.frame(gene,panCancer_up0.05,panCancer_up0.05_ratio,panCancer_up0.05_sample,
                                    panCancer_up0.01,panCancer_up0.01_ratio,panCancer_up0.01_sample,
                                    panCancer_down0.05,panCancer_down0.05_ratio,panCancer_down0.05_sample,
                                    panCancer_down0.01,panCancer_down0.01_ratio,panCancer_down0.01_sample)
      }
      
      if(i%%3000==0){
        CRC_final <- rbind(CRC_final,CRC_tmp)
        STAD_final <- rbind(STAD_final,STAD_tmp)
        panCancer_final <- rbind(panCancer_final,panCancer_tmp)
        final <- cbind(CRC_final,STAD_final,panCancer_final)
        all_final <- rbind(all_final,final)
        write.csv(CRC_final,paste0(outdir,i/3000,"_CRC3000.csv"))
        write.csv(STAD_final,paste0(outdir,i/3000,"_STAD3000.csv"))
        write.csv(panCancer_final,paste0(outdir,i/3000,"_panCancer3000.csv"))
        write.csv(final,paste0(outdir,i/3000,"_3000.csv"))
        
        CRC_final <- as.data.frame(matrix(numeric(0),ncol=17))
        colnames(CRC_final) <- c("gene","FDR_up0.05","FDR_up0.01","FDR_down0.05","FDR_down0.01",
                                 "CRC_up0.05","CRC_up0.05_ratio","CRC_up0.05_sample",
                                 "CRC_up0.01","CRC_up0.01_ratio","CRC_up0.01_sample",
                                 "CRC_down0.05","CRC_down0.05_ratio","CRC_down0.05_sample",
                                 "CRC_down0.01","CRC_down0.01_ratio","CRC_down0.01_sample")
        STAD_final <- as.data.frame(matrix(numeric(0),ncol=13))
        colnames(STAD_final) <- c("gene","STAD_up0.05","STAD_up0.05_ratio","STAD_up0.05_sample",
                                  "STAD_up0.01","STAD_up0.01_ratio","STAD_up0.01_sample",
                                  "STAD_down0.05","STAD_down0.05_ratio","STAD_down0.05_sample",
                                  "STAD_down0.01","STAD_down0.01_ratio","STAD_down0.01_sample")
        panCancer_final <- as.data.frame(matrix(numeric(0),ncol=13))
        colnames(panCancer_final) <- c("gene","panCancer_up0.05","panCancer_up0.05_ratio","panCancer_up0.05_sample",
                                       "panCancer_up0.01","panCancer_up0.01_ratio","panCancer_up0.01_sample",
                                       "panCancer_down0.05","panCancer_down0.05_ratio","panCancer_down0.05_sample",
                                       "panCancer_down0.01","panCancer_down0.01_ratio","panCancer_down0.01_sample")
        i=i+1
      } else if(i==length(gene_list)){
        CRC_final <- rbind(CRC_final,CRC_tmp)
        STAD_final <- rbind(STAD_final,STAD_tmp)
        panCancer_final <- rbind(panCancer_final,panCancer_tmp)
        final <- cbind(CRC_final,STAD_final,panCancer_final)
        all_final <- rbind(all_final,final)
        
        write.csv(CRC_final,paste0(outdir,i%/%3000+1,"_CRC",i%%3000,".csv"))
        write.csv(STAD_final,paste0(outdir,i%/%3000+1,"_STAD",i%%3000,".csv"))
        write.csv(panCancer_final,paste0(outdir,i%/%3000+1,"_panCancer",i%%3000,".csv"))
        write.csv(final,paste0(outdir,i%/%3000+1,"_",i%%3000,".csv"))
        
        i=i+1
      } else {
        CRC_final <- rbind(CRC_final,CRC_tmp)
        STAD_final <- rbind(STAD_final,STAD_tmp)
        panCancer_final <- rbind(panCancer_final,panCancer_tmp)
        i=i+1
      }
      pb$tick()
      Sys.sleep(1 / 100)
    }
    write.csv(all_final,paste0(outdir,"all.csv"))
  }
  
  #make candidate list
  {
    gene_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 4/All_genes.csv",header = TRUE)
    values <- as.character(lapply(strsplit(as.character(gene_list$Gene.ID),".",fixed = TRUE),function(x) x[1]))
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    candidate_list_raw <- getBM(attributes=c("hgnc_symbol","description","ensembl_gene_id","chromosome_name","start_position","end_position","strand"), 
                                filters = "ensembl_gene_id", 
                                values=unique(values), mart= mart,useCache = FALSE)
    candidate_list_raw$strand <- gsub("-1","-",fixed = TRUE,candidate_list_raw$strand)
    candidate_list_raw$strand <- gsub("1","+",fixed = TRUE,candidate_list_raw$strand)
    #candidate_list_raw <- candidate_list_raw[-grep("_",fixed=TRUE,candidate_list_raw$chromosome_name),]
    candidate_list_raw$chromosome_name <- paste0("chr",candidate_list_raw$chromosome_name)
    colnames(candidate_list_raw) <- gsub("hgnc_symbol","Gene",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("description","Description",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("ensembl_gene_id","ensembl",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("chromosome_name","chromosome",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("start_position","start",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("end_position","end",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("strand","strand",fixed = TRUE,colnames(candidate_list_raw))
    write.csv(candidate_list_raw,"/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 4/All_genes_full_annotated.csv",row.names = FALSE)
  }
  {
    gene_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/78_ImmuneModulatory_Genes.csv",header = TRUE)
    #gene_list <- gene_list[which(gene_list$Super.Category=="Antigen presentation"),]
    values <- as.character(lapply(strsplit(as.character(gene_list$HGNC.Symbol),".",fixed = TRUE),function(x) x[1]))
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    candidate_list_raw <- getBM(attributes=c("hgnc_symbol","description","ensembl_gene_id","chromosome_name","start_position","end_position","strand"), 
                                filters = "hgnc_symbol", 
                                values=unique(values), mart= mart,useCache = FALSE)
    candidate_list_raw$strand <- gsub("-1","-",fixed = TRUE,candidate_list_raw$strand)
    candidate_list_raw$strand <- gsub("1","+",fixed = TRUE,candidate_list_raw$strand)
    candidate_list_raw <- candidate_list_raw[-grep("_",fixed=TRUE,candidate_list_raw$chromosome_name),]
    candidate_list_raw$chromosome_name <- paste0("chr",candidate_list_raw$chromosome_name)
    colnames(candidate_list_raw) <- gsub("hgnc_symbol","Gene",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("description","Description",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("ensembl_gene_id","ensembl",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("chromosome_name","chromosome",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("start_position","start",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("end_position","end",fixed = TRUE,colnames(candidate_list_raw))
    colnames(candidate_list_raw) <- gsub("strand","strand",fixed = TRUE,colnames(candidate_list_raw))
    write.csv(candidate_list_raw,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/78_ImmuneModulatory_Genes_full.csv",row.names = FALSE)
  }
  
  #PCAWG plot
  {
    #load functions
    {
      #select best performance
      select_bestperformance <- function(data,level,group){
        print("data: rows are chr_position, or chr|position, chromsome and position are seperated by '_' or '|'. Or gene ensembl ids.")
        print("level: confident level 0.05 or 0.01.")
        print("group: CRC or STAD or panCancer")
        
        data$BP_level <- NA
        data$BP_level_ratio <- NA
        data$BP_level_sample <- NA
        
        i=1
        while(i<=nrow(data)){
          up_ratio <- data[i,which(colnames(data)==paste0(group,"_up",level,"_ratio"))]
          down_ratio <- data[i,which(colnames(data)==paste0(group,"_down",level,"_ratio"))]
          up_n <- which(colnames(data)==paste0(group,"_up",level,"_sample"))
          down_n <- which(colnames(data)==paste0(group,"_down",level,"_sample"))
          up_sample <- as.character(data[i,as.numeric(up_n)])
          down_sample  <- as.character(data[i,as.numeric(down_n)])
          up <- data[i,which(colnames(data)==paste0(group,"_up",level))]
          down <- data[i,which(colnames(data)==paste0(group,"_down",level))]
          
          if(up_ratio>=down_ratio){
            data[i,which(colnames(data)=="BP_level")] <- up
            data[i,which(colnames(data)=="BP_level_ratio")] <- up_ratio
            data[i,which(colnames(data)=="BP_level_sample")] <- up_sample
          } else {
            data[i,which(colnames(data)=="BP_level")] <- down
            data[i,which(colnames(data)=="BP_level_ratio")] <- down_ratio
            data[i,which(colnames(data)=="BP_level_sample")] <- down_sample
          }
          
          #show remain jobs
          if(i%%10000==0){
            print(paste0("Finished: ",round(i/nrow(data)*100,2),"%"))
          } else if(i==nrow(data)) {
            print("Finished!")
          } 
          i=i+1
          
        }
        
        cutoff_tmp <- paste0("BP_",group,"_",level)
        colnames(data) <- c(colnames(data[,1:(ncol(data)-3)]),cutoff_tmp,paste0(cutoff_tmp,"_ratio"),paste0(cutoff_tmp,"_sample"))
        return(data)
      }
      
      #chromosome location2ensembl
      location2ensembl <- function(data,candidate_list){
        print("data: rows are chr_position, or chr|position, chromsome and position are seperated by '_' or '|'.")
        print("list must contain 5 columns: Gene, ensembl, chromosome, start, end.")
        
        #prepare input matrix rownames
        rownames(data) <- gsub("|","_",fixed=TRUE,rownames(data))
        data$gene <- as.character(data$gene)
        data$chromosome <- as.character(lapply(strsplit(as.character(rownames(data)),"_",fixed = TRUE),function(x) x[1]))
        data$position <- as.integer(lapply(strsplit(as.character(rownames(data)),"_",fixed = TRUE),function(x) x[2]))
        
        i=1
        while(i<=nrow(candidate_list)){
          candidate <- candidate_list$ensembl[i]
          name <- candidate_list$Gene[i]
          chromosome <- as.character(candidate_list$chromosome[i])
          start <- candidate_list$start[i]
          end <- candidate_list$end[i]
          
          if(nrow(data[data$chromosome %in% chromosome & data$position>=start & data$position<=end & !is.na(data$position),])==0) {
            print(paste0(i,": No ",candidate,":",name," in this dataset."))
            #data[data$chromosome %in% chromosome & data$position>=start & data$position<=end & !is.na(data$position), which(colnames(data)=="gene")] <- paste(name,"|",candidate)
            i=i+1
          } else {
            data[data$chromosome %in% chromosome & data$position>=start & data$position<=end & !is.na(data$position), which(colnames(data)=="gene")] <- paste(name,"|",candidate)
            i=i+1
          } 
        }
        data <- data[grep("ENSG",data$gene),]
        return(data)
      }
      
      #subset by ensembl id
      subset_row <- function(data,candidate_list,cutoff){
        print("data: rows are gene ensembl ids.")
        print("list must contain 2 columns: Gene, ensembl.")
        i=1
        output_matrix <- {}
        while(i<=nrow(candidate_list)){
          candidate <- candidate_list$ensembl[i]
          name <- candidate_list$Gene[i]
          if(length(grep(candidate,data$gene))==0) {
            #print(paste0("No ",candidate," in this dataset."))
            temp <- as.data.frame(array(,dim=c(1,ncol(data))))
            temp[1,] <- NA
            colnames(temp) <- colnames(data)
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          } else if(length(grep(candidate,data$gene))==1) {
            temp <- data[grep(candidate,data$gene,fixed=TRUE),]
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          } else {
            ##max
            minidata <- data[grep(candidate,data$gene),]
            minidata <- minidata[order(minidata[,cutoff],decreasing = TRUE),]
            temp <- minidata[1,]
            #average
            #
            #union
            #
            #take 1st as represent
            #temp <- data[grep(candidate,data$gene,fixed=TRUE)[1],]
            
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          }
        }
        return(output_matrix)
      }
      
      #chimeric subset by gene_symbol
      chimeric_subset_row <- function(data,candidate_list,cutoff){
        print("data: rows are gene enembl ids.")
        print("list must contain 2 columns: Gene, ensembl.")
        i=1
        output_matrix <- {}
        while(i<=nrow(candidate_list)){
          candidate <- candidate_list$ensembl[i]
          name <- candidate_list$Gene[i]
          
          toMatch <- c(paste0("^",name,"\\|"),paste0("--",name,"\\|"))
          pattern <- paste(toMatch,collapse="|")
          
          #pattern <- paste0("^",candidate,"\\|")
          
          if(length(grep(pattern,data$gene))==0) {
            #print(paste0("No ",candidate," in this dataset."))
            temp <- as.data.frame(array(,dim=c(1,ncol(data))))
            temp[1,] <- NA
            colnames(temp) <- colnames(data)
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          } else if(length(grep(pattern,data$gene))==1) {
            temp <- data[grep(pattern,data$gene),]
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          } else {
            ##max
            minidata <- data[grep(pattern,data$gene),]
            minidata <- minidata[order(minidata[,cutoff],decreasing = TRUE),]
            temp <- minidata[1,]
            ##average
            #
            ##union
            #
            ##take 1st as represent
            #temp <- data[grep(pattern,data$gene)[1],]
            
            rownames(temp) <- paste(name,"|",candidate)
            output_matrix <- rbind(output_matrix,temp)
            i=i+1
          }
        }
        return(output_matrix)
      }
    }
    
    #settings
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 4/Detection ratio result/")
    #candidate_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/06.All features/Colorectal Cancer Genes.csv",header = TRUE)
    candidate_list <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 4/All_genes_full_annotated.csv",header = TRUE)
    genes <- "All_Genes"
    #choose a confidence level
    level <- 0.05 #c(0.05,0.01)
    #choose an interested cancer type
    group <- "panCancer" #c("CRC","STAD","panCancer")
    #choose a cutoff, then subset
    cutoff <- "BP_panCancer_0.05" #c("BP_CRC_0.05","BP_CRC_0.01","BP_STAD_0.05","BP_STAD_0.01","BP_panCancer_0.05","BP_panCancer_0.01","CRC_up0.05","CRC_up0.01","CRC_down0.05","CRC_down0.01","STAD_up0.05","STAD_up0.01","STAD_down0.05","STAD_down0.01","panCancer_up0.05","panCancer_up0.01","panCancer_down0.05","panCancer_down0.01")
    
    #give the total sample number
    
    #CRC sample number
    #{
    #sampleN_RNA <- 23
    #sampleN_DNA <- 23
    #}
    #STAD sample number
    #{
    #  sampleN_RNA <- 30
    #  sampleN_DNA <- 30
    #}
    
    #panCancer sample number
    {
      sampleN_RNA <- 53
      sampleN_DNA <- 53
    }
    
    #read in feature outliers
    {
      #Alt.promoter
      Alt <- read.csv("AlternativePromoter/all.csv",header = TRUE)
      Alt <- Alt[,-1]
      rownames(Alt) <- Alt$gene
      Alt <- select_bestperformance(Alt,level,group)
      Alt_targetgenes <- subset_row(Alt,candidate_list,cutoff)
      
      #Expression
      Expression <- read.csv("Expression/all.csv",header = TRUE)
      Expression <- Expression[,-1]
      rownames(Expression) <- Expression$gene
      Expression <- select_bestperformance(Expression,level,group)
      Expression_targetgenes <- subset_row(Expression,candidate_list,cutoff)
      
      #Splicing
      Splicing <- read.csv("Splicing/all.csv",header = TRUE)
      Splicing <- Splicing[,-1]
      rownames(Splicing) <- Splicing$gene
      Splicing <- select_bestperformance(Splicing,level,group)
      Splicing_targetgenes <- subset_row(Splicing,candidate_list,cutoff)
      
      #APA
      #APA <- read.csv("APA/all.csv",header = TRUE)
      #APA <- APA[,-1]
      #rownames(APA) <- APA$gene
      #APA <- select_bestperformance(APA,level,group)
      #APA_targetgenes <- subset_row(APA,candidate_list,cutoff)
      
      #chimeric
      chimeric <- read.csv("chimeric/all.csv",header = TRUE)
      chimeric <- chimeric[,-1]
      rownames(chimeric) <- chimeric$gene
      chimeric <- select_bestperformance(chimeric,level,group)
      chimeric_targetgenes <- chimeric_subset_row(chimeric,candidate_list,cutoff)
      
      #Editing
      Editing <- read.csv("Editing/all.csv",header = TRUE)
      Editing <- Editing[,-1]
      rownames(Editing) <- Editing$gene
      Editing <- location2ensembl(Editing,candidate_list)
      Editing <- select_bestperformance(Editing,level,group)
      Editing_targetgenes <- subset_row(Editing,candidate_list,cutoff)
      
      #ASE
      ASE <- read.csv("ASE/all.csv",header = TRUE)
      ASE <- ASE[,-1]
      rownames(ASE) <- ASE$gene
      ASE <- location2ensembl(ASE,candidate_list)
      ASE <- select_bestperformance(ASE,level,group)
      ASE_targetgenes <- subset_row(ASE,candidate_list,cutoff)
      
      #SNP
      SNP <- read.csv("SNP/all.csv",header = TRUE)
      SNP <- SNP[,-1]
      rownames(SNP) <- SNP$gene
      SNP <- location2ensembl(SNP,candidate_list)
      SNP <- select_bestperformance(SNP,level,group)
      SNP_targetgenes <- subset_row(SNP,candidate_list,cutoff)
      
      #DNA-CNV
      DNA_CNV <- read.csv("CNV/all.csv",header = TRUE)
      DNA_CNV <- DNA_CNV[,-1]
      rownames(DNA_CNV) <- DNA_CNV$gene
      DNA_CNV <- select_bestperformance(DNA_CNV,level,group)
      DNA_CNV_targetgenes <- subset_row(DNA_CNV,candidate_list,cutoff)
      
      #DNA-WPS
      DNA_WPS <- read.csv("WPS/all.csv",header = TRUE)
      DNA_WPS <- DNA_WPS[,-1]
      rownames(DNA_WPS) <- DNA_WPS$gene
      DNA_WPS <- select_bestperformance(DNA_WPS,level,group)
      DNA_WPS_targetgenes <- subset_row(DNA_WPS,candidate_list,cutoff)
      
      #MeDIP
      DNA_MeDIP <- read.csv("Medip_promoter/all.csv",header = TRUE)
      DNA_MeDIP <- DNA_MeDIP[,-1]
      rownames(DNA_MeDIP) <- DNA_MeDIP$gene
      DNA_MeDIP <- select_bestperformance(DNA_MeDIP,level,group)
      DNA_MeDIP_targetgenes <- subset_row(DNA_MeDIP,candidate_list,cutoff)
      
    }
    
    
    #read in alterations outliers (only for whole genome candidate list)
    {
      Alt_targetgenes <- read.csv("AlternativePromoter/Alt_targetgenes.csv",row.names = 1)
      Editing_targetgenes <- read.csv("Editing/Editing_targetgenes.csv",row.names = 1)
      Expression_targetgenes <- read.csv("Expression/Expression_targetgenes.csv",row.names = 1)
      Splicing_targetgenes <- read.csv("Splicing/Splicing_targetgenes.csv",row.names = 1)
      chimeric_targetgenes <- read.csv("chimeric/chimeric_targetgenes.csv",row.names = 1)
      ASE_targetgenes <- read.csv("ASE/ASE_targetgenes.csv",row.names = 1)
      SNP_targetgenes <- read.csv("SNP/SNP_targetgenes.csv",row.names = 1)
      DNA_CNV_targetgenes <- read.csv("CNV/CNV_targetgenes.csv",row.names = 1)
      DNA_WPS_targetgenes <- read.csv("WPS/WPS_targetgenes.csv",row.names = 1)
      DNA_MeDIP_targetgenes <- read.csv("Medip_promoter/MeDIP_targetgenes.csv",row.names = 1)
    }
    
    #plot
    {
      #summarize outliers (RNA+DNA merged barplot for complementary)
      {
        #input features
        Alt_cutoff <- Alt_targetgenes[grep(cutoff,colnames(Alt_targetgenes))]
        Alt_cutoff[,3] <- gsub(".pico", "", Alt_cutoff[,3])
        
        Expression_cutoff <- Expression_targetgenes[grep(cutoff,colnames(Expression_targetgenes))]
        Expression_cutoff[,3] <- gsub(".pico", "", Expression_cutoff[,3])
        
        Splicing_cutoff <- Splicing_targetgenes[grep(cutoff,colnames(Splicing_targetgenes))]
        Splicing_cutoff[,3] <- gsub(".pico", "", Splicing_cutoff[,3])
        
        #APA_cutoff <- APA_targetgenes[grep(cutoff,colnames(APA_targetgenes))]
        #APA_cutoff[,3] <- gsub(".pico", "", APA_cutoff[,3])
        
        chimeric_cutoff <- chimeric_targetgenes[grep(cutoff,colnames(chimeric_targetgenes))]
        chimeric_cutoff[,3] <- gsub(".pico", "", chimeric_cutoff[,3])
        
        Editing_cutoff <- Editing_targetgenes[grep(cutoff,colnames(Editing_targetgenes))]
        Editing_cutoff[,3] <- gsub(".pico", "", Editing_cutoff[,3])
        
        ASE_cutoff <- ASE_targetgenes[grep(cutoff,colnames(ASE_targetgenes))]
        ASE_cutoff[,3] <- gsub(".pico", "", ASE_cutoff[,3])
        
        SNP_cutoff <- SNP_targetgenes[grep(cutoff,colnames(SNP_targetgenes))]
        SNP_cutoff[,3] <- gsub(".pico", "", SNP_cutoff[,3])
        
        DNA_CNV_cutoff <- DNA_CNV_targetgenes[grep(cutoff,colnames(DNA_CNV_targetgenes))]
        DNA_CNV_cutoff[,3] <- gsub(".wgs", "", DNA_CNV_cutoff[,3])
        
        DNA_WPS_cutoff <- DNA_WPS_targetgenes[grep(cutoff,colnames(DNA_WPS_targetgenes))]
        DNA_WPS_cutoff[,3] <- gsub(".wgs", "", DNA_WPS_cutoff[,3])
        
        DNA_MeDIP_cutoff <- DNA_MeDIP_targetgenes[grep(cutoff,colnames(DNA_MeDIP_targetgenes))]
        DNA_MeDIP_cutoff[,3] <- gsub(".me", "", DNA_MeDIP_cutoff[,3])
        
        features <- cbind(Alt_cutoff,Expression_cutoff,Splicing_cutoff,
                          chimeric_cutoff,Editing_cutoff,ASE_cutoff,SNP_cutoff,
                          DNA_CNV_cutoff,DNA_WPS_cutoff,DNA_MeDIP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
        
        features$UnionN <- NA 
        features$UnionN_sample <- NA
        features$UnionN_early <- NA 
        features$UnionN_early_sample <- NA 
        
        #Union number of input features
        i=1
        while(i<=nrow(features)){
          x <- features[i,grep(paste0(cutoff,"_sample"),colnames(features))]
          print(x)
          y <- lapply(x,as.character)
          z <- paste(y,sep = "\n",collapse = "\n")
          u <- unique(unlist(strsplit(z,"\n")))
          u <- u[u!=""]
          u <- u[u!="NA"]
          features[i,which(colnames(features)=="UnionN")] <- length(u)
          features[i,which(colnames(features)=="UnionN_sample")] <- paste(u,collapse = "\n")
          
          #remove late stage samples
          u <- gsub("STAD-PKU-11","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-12","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-13","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-16","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-17","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-18","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-19","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-25","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-27","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-28","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-29","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-34","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-35","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-36","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-37","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-38","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-10","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-27","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-29","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-30","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-6","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-9","",u,fixed = TRUE)
          u <- unique(unlist(strsplit(u,"\n")))
          u <- u[u!=""]
          u <- u[u!="NA"]
          features[i,which(colnames(features)=="UnionN_early")] <- length(u)
          features[i,which(colnames(features)=="UnionN_early_sample")] <- paste(u,collapse = "\n")
          
          i=i+1
        }
        
        #melt
        library(reshape)
        Outlier_tmp <- features[,-grep("_sample",colnames(features))]
        Outlier <- Outlier_tmp[,-grep("_ratio",colnames(Outlier_tmp))]
        colnames(Outlier) <- c("Alternative promoter outliers","Expression outliers","Splice outliers",
                               "Chimeric RNA","RNA editing outliers","ASE outliers","SNP",
                               "DNA: Copy-number outlier","DNA: Window protection score outliers","DNA: Methylation outliers","UnionN","UnionN_early")
        Outlier[is.na(Outlier)] <- 0
        Outlier$Gene <- rownames(Outlier)
        Outlier_plot <- melt(Outlier)
        Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
        colnames(Outlier_plot) <- c("Gene","Feature","OutlierN","Gene_symbol")
        
        write.csv(Outlier_plot,paste0("./",genes,"_",group,"_Outlier_plot_all_with_early_stage.csv"),quote = FALSE)
        
        #plot
        #Union
        Union <- Outlier_plot[grep("UnionN",Outlier_plot$Feature),]
        Union$Gene <- factor(Union$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        Union$Gene_symbol <- factor(Union$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        p_Union <- ggplot(Union,aes(x=Gene_symbol,y=OutlierN/sampleN_RNA,fill = Feature))+
          scale_fill_manual(values = c("#666666"))+
          geom_bar(stat = "identity")+
          xlab("")+
          ylab(paste0("RNA+DNA Outliers: ",cutoff,"\nTotal Number: ",sampleN_RNA))+
          guides(fill=FALSE)+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(10,5,-10,5),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            axis.line.y = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            #axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            axis.text.y = element_text(face="bold",  color="black", size=18, angle = 0,hjust=1,vjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0.0","","0.5","","1.0"),expand = c(0,0),limits = c(0,1))
        
        
        #different features contribution
        Features <- Outlier_plot[-grep("UnionN",Outlier_plot$Feature),]
        Features$Gene <- factor(Features$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        Features$Gene_symbol <- factor(Features$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        p_Features <- ggplot(Features,aes(x=Gene_symbol,y=OutlierN/sampleN_RNA,fill=Feature))+
          geom_bar(stat = "identity",position = "fill")+
          xlab("")+
          ylab(paste0("Features: ",cutoff,"\nTotal Number: ",sampleN_RNA))+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(5,5,5,5),units="pt"),
            legend.position="bottom",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            #axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            #axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0.01),limits = c(0,1))+
          scale_fill_simpsons(alpha = 0.8)
        
        library(ggpubr)
        ggarrange(p_Union, p_Features,
                  ncol = 1, nrow = 2,heights = c(5,5),align = c("v"))
      }  
      
      #summarize outliers (RNA)
      {
        #input features
        Alt_cutoff <- Alt_targetgenes[grep(cutoff,colnames(Alt_targetgenes))]
        Expression_cutoff <- Expression_targetgenes[grep(cutoff,colnames(Expression_targetgenes))]
        Splicing_cutoff <- Splicing_targetgenes[grep(cutoff,colnames(Splicing_targetgenes))]
        
        #APA_cutoff <- APA_targetgenes[grep(cutoff,colnames(APA_targetgenes))]
        chimeric_cutoff <- chimeric_targetgenes[grep(cutoff,colnames(chimeric_targetgenes))]
        
        Editing_cutoff <- Editing_targetgenes[grep(cutoff,colnames(Editing_targetgenes))]
        ASE_cutoff <- ASE_targetgenes[grep(cutoff,colnames(ASE_targetgenes))]
        SNP_cutoff <- SNP_targetgenes[grep(cutoff,colnames(SNP_targetgenes))]
        
        features <- cbind(Alt_cutoff,Expression_cutoff,Splicing_cutoff,chimeric_cutoff,Editing_cutoff,ASE_cutoff,SNP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
        
        features$UnionN <- NA 
        features$UnionN_sample <- NA
        features$UnionN_early <- NA 
        features$UnionN_early_sample <- NA 
        
        #Union number of input features
        i=1
        while(i<=nrow(features)){
          x <- features[i,grep(paste0(cutoff,"_sample"),colnames(features))]
          print(x)
          y <- lapply(x,as.character)
          z <- paste(y,sep = "\n",collapse = "\n")
          u <- unique(unlist(strsplit(z,"\n")))
          u <- u[u!=""]
          u <- u[u!="NA"]
          features[i,which(colnames(features)=="UnionN")] <- length(u)
          features[i,which(colnames(features)=="UnionN_sample")] <- paste(u,collapse = "\n")
          
          #remove late stage samples
          u <- gsub("STAD-PKU-11","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-12","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-13","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-16","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-17","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-18","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-19","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-25","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-27","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-28","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-29","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-34","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-35","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-36","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-37","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-38","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-10","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-27","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-29","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-30","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-6","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-9","",u,fixed = TRUE)
          u <- unique(unlist(strsplit(u,"\n")))
          u <- u[u!=""]
          u <- u[u!="NA"]
          features[i,which(colnames(features)=="UnionN_early")] <- length(u)
          features[i,which(colnames(features)=="UnionN_early_sample")] <- paste(u,collapse = "\n")
          
          i=i+1
        }
        
        #melt
        library(reshape)
        Outlier_tmp <- features[,-grep("_sample",colnames(features))]
        Outlier <- Outlier_tmp[,-grep("_ratio",colnames(Outlier_tmp))]
        colnames(Outlier) <- c("Alternative promoter outliers","Expression outliers","Splice outliers","Chimeric RNA","RNA editing outliers","ASE outliers","SNP","UnionN","UnionN_early")
        Outlier[is.na(Outlier)] <- 0
        Outlier$Gene <- rownames(Outlier)
        Outlier_plot <- melt(Outlier)
        Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
        colnames(Outlier_plot) <- c("Gene","Feature","OutlierN","Gene_symbol")
        
        write.csv(Outlier_plot,paste0("./",genes,"_",group,"_Outlier_plot_all_RNA_with_early.csv"),quote = FALSE)
        
        #plot
        #Union
        RNA_Union <- Outlier_plot[grep("UnionN",Outlier_plot$Feature),]
        RNA_Union$Gene <- factor(Union$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        RNA_Union$Gene_symbol <- factor(Union$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        RNA_Union$background <- NA
        RNA_Union$RNA_unique <- NA
          
        background <- Union
        background$background <- Union$OutlierN
        background$Feature <- "Background"
        background$OutlierN <- NA
        background$RNA_unique <- NA
        
        RNA_unique <- Union
        RNA_unique$background <- NA
        RNA_unique$Feature <- "RNA_unique"
        RNA_unique$OutlierN <- NA
        
        {
          #input features
          DNA_CNV_cutoff <- DNA_CNV_targetgenes[grep(cutoff,colnames(DNA_CNV_targetgenes))]
          DNA_CNV_cutoff[,3] <- gsub(".wgs", "", DNA_CNV_cutoff[,3])
          DNA_WPS_cutoff <- DNA_WPS_targetgenes[grep(cutoff,colnames(DNA_WPS_targetgenes))]
          DNA_WPS_cutoff[,3] <- gsub(".wgs", "", DNA_WPS_cutoff[,3])
          DNA_MeDIP_cutoff <- DNA_MeDIP_targetgenes[grep(cutoff,colnames(DNA_MeDIP_targetgenes))]
          DNA_MeDIP_cutoff[,3] <- gsub(".me", "", DNA_MeDIP_cutoff[,3])
          
          features <- cbind(DNA_CNV_cutoff,DNA_WPS_cutoff,DNA_MeDIP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
          features$UnionN <- NA 
          
          #Union number of input features
          i=1
          while(i<=nrow(features)){
            x <- features[i,grep(paste0(cutoff,"_sample"),colnames(features))]
            y <- lapply(x,as.character)
            z <- paste(y,sep = "\n",collapse = "\n")
            u <- unique(unlist(strsplit(z,"\n")))
            u <- u[u!=""]
            u <- u[u!="NA"]
            features[i,which(colnames(features)=="UnionN")] <- length(u)
            i=i+1
          }
          
          #melt
          library(reshape)
          Outlier_tmp <- features[,-grep("_sample",colnames(features))]
          Outlier <- Outlier_tmp[,-grep("_ratio",colnames(Outlier_tmp))]
          colnames(Outlier) <- c("DNA: Copy-number outliers","DNA: Window protection score outliers","DNA: Methylation outliers","UnionN")
          Outlier[is.na(Outlier)] <- 0
          Outlier$Gene <- rownames(Outlier)
          Outlier_plot <- melt(Outlier)
          Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
          colnames(Outlier_plot) <- c("Gene","Feature","OutlierN","Gene_symbol")
          
          
          #plot
          #Union
          DNA_Union <- Outlier_plot[grep("UnionN",Outlier_plot$Feature),]
          DNA_Union$Gene <- factor(DNA_Union$Gene,levels=DNA_Union[order(DNA_Union$OutlierN,decreasing = FALSE),]$Gene)
          DNA_Union$Gene_symbol <- factor(DNA_Union$Gene_symbol,levels=DNA_Union[order(DNA_Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        }
        RNA_unique$RNA_unique <- background$background-DNA_Union$OutlierN
        
        RNA_Union <- rbind(RNA_Union,background,RNA_unique)
        p_Union_RNA <- ggplot(RNA_Union)+
          geom_bar(aes(x=Gene_symbol,y=background/sampleN_RNA,fill = Feature),stat = "identity")+
          geom_bar(aes(x=Gene_symbol,y=OutlierN/sampleN_RNA,fill = Feature),stat = "identity")+
          geom_bar(aes(x=Gene_symbol,y=RNA_unique/sampleN_RNA,fill = Feature),stat = "identity")+
          scale_fill_manual(values = c("UnionN"=alpha("black",alpha = 0.5),"Background"=alpha("#666666",alpha = 0.3),"RNA_unique"=alpha("black",alpha = 1)))+
          #scale_fill_manual(values = c("UnionN"="#138F6A","Background"=alpha("#666666",alpha = 0.3)))+
          xlab("")+
          ylab(paste0("RNA Outliers: ",cutoff,"\nTotal Number: ",sampleN_RNA))+
          guides(fill=FALSE)+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(10,5,-10,5),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            axis.line.y = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            #axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            axis.text.y = element_text(face="bold",  color="black", size=18, angle = 0,hjust=1,vjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0.0","","0.5","","1.0"),expand = c(0,0),limits = c(0,1))
        
        
        #different features contribution
        RNA_Features <- Outlier_plot[-grep("UnionN",Outlier_plot$Feature),]
        RNA_Features$Gene <- factor(RNA_Features$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        RNA_Features$Gene_symbol <- factor(RNA_Features$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        p_Features_RNA <- ggplot(RNA_Features,aes(x=Gene_symbol,y=OutlierN/sampleN_RNA,fill=Feature))+
          geom_bar(stat = "identity",position = "fill")+
          xlab("")+
          ylab(paste0("Features: ",cutoff,"\nTotal Number: ",sampleN_RNA))+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(5,5,5,5),units="pt"),
            legend.position="bottom",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            #axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            #axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0.01),limits = c(0,1))+
          scale_fill_jco(alpha = 0.8)
        
        library(ggpubr)
        #ggarrange(p_Union_RNA, p_Features_RNA,
        #          ncol = 1, nrow = 2,heights = c(5,5),align = c("v"))
      }
      
      #summarize outliers (DNA)
      {
        #input features
        DNA_CNV_cutoff <- DNA_CNV_targetgenes[grep(cutoff,colnames(DNA_CNV_targetgenes))]
        DNA_CNV_cutoff[,3] <- gsub(".wgs", "", DNA_CNV_cutoff[,3])
        DNA_WPS_cutoff <- DNA_WPS_targetgenes[grep(cutoff,colnames(DNA_WPS_targetgenes))]
        DNA_WPS_cutoff[,3] <- gsub(".wgs", "", DNA_WPS_cutoff[,3])
        DNA_MeDIP_cutoff <- DNA_MeDIP_targetgenes[grep(cutoff,colnames(DNA_MeDIP_targetgenes))]
        DNA_MeDIP_cutoff[,3] <- gsub(".me", "", DNA_MeDIP_cutoff[,3])
        
        features <- cbind(DNA_CNV_cutoff,DNA_WPS_cutoff,DNA_MeDIP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
        features$UnionN <- NA 
        features$UnionN_sample <- NA
        features$UnionN_early <- NA 
        features$UnionN_early_sample <- NA 
        
        #Union number of input features
        i=1
        while(i<=nrow(features)){
          x <- features[i,grep(paste0(cutoff,"_sample"),colnames(features))]
          print(x)
          y <- lapply(x,as.character)
          z <- paste(y,sep = "\n",collapse = "\n")
          u <- unique(unlist(strsplit(z,"\n")))
          u <- u[u!=""]
          u <- u[u!="NA"]
          features[i,which(colnames(features)=="UnionN")] <- length(u)
          features[i,which(colnames(features)=="UnionN_sample")] <- paste(u,collapse = "\n")
          
          #remove late stage samples
          u <- gsub("STAD-PKU-11","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-12","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-13","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-16","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-17","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-18","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-19","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-25","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-27","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-28","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-29","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-34","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-35","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-36","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-37","",u,fixed = TRUE)
          u <- gsub("STAD-PKU-38","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-10","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-27","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-29","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-30","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-6","",u,fixed = TRUE)
          u <- gsub("CRC-PKU-9","",u,fixed = TRUE)
          u <- unique(unlist(strsplit(u,"\n")))
          u <- u[u!=""]
          u <- u[u!="NA"]
          features[i,which(colnames(features)=="UnionN_early")] <- length(u)
          features[i,which(colnames(features)=="UnionN_early_sample")] <- paste(u,collapse = "\n")
          
          i=i+1
        }
        
        #melt
        library(reshape)
        Outlier_tmp <- features[,-grep("_sample",colnames(features))]
        Outlier <- Outlier_tmp[,-grep("_ratio",colnames(Outlier_tmp))]
        colnames(Outlier) <- c("DNA: Copy-number outliers","DNA: Window protection score outliers","DNA: Methylation outliers","UnionN","UnionN_early")
        Outlier[is.na(Outlier)] <- 0
        Outlier$Gene <- rownames(Outlier)
        Outlier_plot <- melt(Outlier)
        Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
        colnames(Outlier_plot) <- c("Gene","Feature","OutlierN","Gene_symbol")
        
        write.csv(Outlier_plot,paste0("./",genes,"_",group,"_Outlier_plot_all_DNA_with_early_stage.csv"),quote = FALSE)
        
        #plot
        #Union
        DNA_Union <- Outlier_plot[grep("UnionN",Outlier_plot$Feature),]
        DNA_Union$Gene <- factor(DNA_Union$Gene,levels=DNA_Union[order(DNA_Union$OutlierN,decreasing = FALSE),]$Gene)
        DNA_Union$Gene_symbol <- factor(DNA_Union$Gene_symbol,levels=DNA_Union[order(DNA_Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        p_Union_DNA <- ggplot(DNA_Union,aes(x=Gene_symbol,y=OutlierN/sampleN_DNA,fill = Feature))+
          scale_fill_manual(values = c("#EEC900"))+
          geom_bar(stat = "identity")+
          xlab("")+
          ylab(paste0("DNA Outliers: ",cutoff,"\nTotal Number: ",sampleN_DNA))+
          guides(fill=FALSE)+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(10,5,-10,5),units="pt"),
            legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.ticks.x = element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            axis.text.x = element_blank(),
            #axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            axis.text.y = element_text(face="bold",  color="black", size=18, angle = 0,hjust=1,vjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0),limits = c(0,1))
        
        #reversed barplot
        {
          DNA_Union$Gene <- factor(DNA_Union$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
          DNA_Union$Gene_symbol <- factor(DNA_Union$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
          DNA_Union$background <- NA
          background <- Union
          background$background <- Union$OutlierN
          background$Feature <- "Background"
          background$OutlierN <- NA
          DNA_Union <- rbind(DNA_Union,background)
          p_Union_DNA_rev <- ggplot(DNA_Union)+
            geom_bar(aes(x=Gene_symbol,y=background/sampleN_DNA,fill = Feature),stat = "identity")+
            geom_bar(aes(x=Gene_symbol,y=OutlierN/sampleN_DNA,fill = Feature),stat = "identity")+
            scale_fill_manual(values = c("UnionN"=alpha("black",alpha = 0.5),"Background"=alpha("#666666",alpha = 0.3)))+
            #scale_fill_manual(values = c("UnionN"="#EEC900","Background"=alpha("#666666",alpha = 0.3)))+
            xlab("")+
            ylab(paste0("DNA Outlier: ",cutoff,"\nTotal Number: ",sampleN_DNA))+
            guides(fill=FALSE)+
            theme_bw()+
            theme(#legend.position="right",
              plot.margin = unit(x=c(5,5,10,5),units="pt"),
              legend.position="right",
              panel.grid=element_blank(),
              panel.border=element_blank(),
              legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
              legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
              plot.title = element_text(hjust = 0.5,size=24,face="bold"),
              axis.ticks.x = element_blank(),
              #axis.text.x = element_blank(),
              axis.line.y = element_line(color = "black"),
              axis.text.x = element_text(color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
              axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
              axis.title.x = element_text(face="bold", color="black", size=20),
              axis.title.y = element_text(face="bold",color="black", size=20))+
            scale_y_reverse(breaks = c(0,0.25,0.50,0.75,1),labels = c("0.0","","0.5","","1.0"),expand = c(0,0),limits = c(1,0))
          }
        
        #different features contribution
        DNA_Features <- Outlier_plot[-grep("UnionN",Outlier_plot$Feature),]
        DNA_Features$Gene <- factor(DNA_Features$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        DNA_Features$Gene_symbol <- factor(DNA_Features$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        p_Features_DNA <- ggplot(DNA_Features,aes(x=Gene_symbol,y=OutlierN/sampleN_DNA,fill=Feature))+
          geom_bar(stat = "identity",position = "fill")+
          xlab("")+
          ylab(paste0("Features: ",cutoff,"\nTotal Number: ",sampleN_DNA))+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(5,5,5,5),units="pt"),
            legend.position="bottom",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            #axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            #axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0.01),limits = c(0,1))+
          scale_fill_aaas(alpha = 0.8)
        
        library(ggpubr)
        #ggarrange(p_Union_DNA, p_Features_DNA,
        #          ncol = 1, nrow = 2,heights = 5,align = c("v"))
      }
      
      #summarize outliers (RNA+DNA)
      {
        #input features
        Alt_cutoff <- Alt_targetgenes[grep(cutoff,colnames(Alt_targetgenes))]
        Alt_cutoff[,3] <- gsub(".pico", "", Alt_cutoff[,3])
        
        Expression_cutoff <- Expression_targetgenes[grep(cutoff,colnames(Expression_targetgenes))]
        Expression_cutoff[,3] <- gsub(".pico", "", Expression_cutoff[,3])
        
        Splicing_cutoff <- Splicing_targetgenes[grep(cutoff,colnames(Splicing_targetgenes))]
        Splicing_cutoff[,3] <- gsub(".pico", "", Splicing_cutoff[,3])
        
        #APA_cutoff <- APA_targetgenes[grep(cutoff,colnames(APA_targetgenes))]
        #APA_cutoff[,3] <- gsub(".pico", "", APA_cutoff[,3])
        
        chimeric_cutoff <- chimeric_targetgenes[grep(cutoff,colnames(chimeric_targetgenes))]
        chimeric_cutoff[,3] <- gsub(".pico", "", chimeric_cutoff[,3])
        
        Editing_cutoff <- Editing_targetgenes[grep(cutoff,colnames(Editing_targetgenes))]
        Editing_cutoff[,3] <- gsub(".pico", "", Editing_cutoff[,3])
        
        ASE_cutoff <- ASE_targetgenes[grep(cutoff,colnames(ASE_targetgenes))]
        ASE_cutoff[,3] <- gsub(".pico", "", ASE_cutoff[,3])
        
        SNP_cutoff <- SNP_targetgenes[grep(cutoff,colnames(SNP_targetgenes))]
        SNP_cutoff[,3] <- gsub(".pico", "", SNP_cutoff[,3])
        
        DNA_CNV_cutoff <- DNA_CNV_targetgenes[grep(cutoff,colnames(DNA_CNV_targetgenes))]
        DNA_CNV_cutoff[,3] <- gsub(".wgs", "", DNA_CNV_cutoff[,3])
        
        DNA_WPS_cutoff <- DNA_WPS_targetgenes[grep(cutoff,colnames(DNA_WPS_targetgenes))]
        DNA_WPS_cutoff[,3] <- gsub(".wgs", "", DNA_WPS_cutoff[,3])
        
        DNA_MeDIP_cutoff <- DNA_MeDIP_targetgenes[grep(cutoff,colnames(DNA_MeDIP_targetgenes))]
        DNA_MeDIP_cutoff[,3] <- gsub(".me", "", DNA_MeDIP_cutoff[,3])
        
        features <- cbind(DNA_CNV_cutoff,DNA_WPS_cutoff,DNA_MeDIP_cutoff,Alt_cutoff,Expression_cutoff,Splicing_cutoff,chimeric_cutoff,Editing_cutoff,ASE_cutoff,SNP_cutoff) ## this sequential is important! It determines the label seqeuntial mannually set at melt step!
        
        #melt
        library(reshape)
        Outlier <- features[,grep("_ratio",colnames(features))]
        colnames(Outlier) <- c("DNA: Copy-number outlier","DNA: Window protection score outliers","DNA: Methylation outliers","Alternative promoter outliers","Expression outliers","Splice outliers","Chimeric RNA","RNA editing outliers","ASE outliers","RNA SNP")
        #colnames(Outlier) <- c("DNA_CNV","DNA_nucleosome","DNA_MeDIP","Altpromoter","Expression","Splicing","APA","chimeric","Editing","ASE","SNP")
        Outlier[is.na(Outlier)] <- 0
        Outlier$Gene <- rownames(Outlier)
        Outlier_plot <- melt(Outlier)
        Outlier_plot$Gene_symbol <- as.character(lapply(strsplit(as.character(Outlier_plot$Gene)," |",fixed = TRUE),function(x) x[1]))
        colnames(Outlier_plot) <- c("Gene","Feature","OutlierRatio","Gene_symbol")
        
        
        #plot
        #different features contribution
        Features <- Outlier_plot
        Features$Gene <- factor(Features$Gene,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene)
        Features$Gene_symbol <- factor(Features$Gene_symbol,levels=Union[order(Union$OutlierN,decreasing = FALSE),]$Gene_symbol)
        Features_color <- c("DNA: Copy-number outlier"="#FFE600",
                            "DNA: Window protection score outliers"="#60AFFE",
                            "DNA: Methylation outliers"="#CDC9C9",
                            "Alternative promoter outliers"="#E9C2A6",
                            "Expression outliers"="#A5435C",
                            "Splice outliers"="#C5E3BF",
                            #"Alternative polyadenyltion outliers"="#003F87",
                            "Chimeric RNA"="#FF3D0D",
                            "RNA editing outliers"="#324F17",
                            "ASE outliers"="#800080",
                            "RNA SNP"="#333333")
        p_Features <- ggplot(Features,aes(x=Gene_symbol,y=OutlierRatio,fill=Feature))+
          geom_bar(stat = "identity",position = "fill")+
          xlab("")+
          ylab("")+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(-10,5,-25,5),units="pt"),
            legend.position="bottom",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            #axis.text.x = element_text(face="bold", color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            #axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0.01),limits = c(0,1))+
          scale_fill_manual(values = Features_color)
        
        library(ggpubr)
        #ggarrange(p_Union_RNA,p_Features, p_Union_DNA_rev,
        #          ncol = 1, nrow = 3,heights = c(5,2.5,7),align = c("v"), legend = "right", common.legend = TRUE)
        
        p_Features <- ggplot(Features,aes(x=Gene_symbol,y=OutlierRatio,fill=Feature))+
          geom_bar(stat = "identity",position = "fill")+
          xlab("")+
          ylab("")+
          theme_bw()+
          theme(#legend.position="right",
            plot.margin = unit(x=c(-10,5,0,5),units="pt"),
            legend.position="bottom",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
            legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
            plot.title = element_text(hjust = 0.5,size=24,face="bold"),
            axis.ticks = element_blank(),
            #axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(color="black", size=20, angle = 90,hjust = 1,vjust = 0.5),
            #axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
            axis.title.x = element_text(face="bold", color="black", size=20),
            axis.title.y = element_text(face="bold",color="black", size=20))+
          scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0.01),limits = c(0,1))+
          scale_fill_manual(values = Features_color)
        #ggarrange(p_Union_RNA,p_Features,
        #          ncol = 1, nrow = 2,heights = c(5,5),align = c("v"), legend = "right", common.legend = TRUE)
      }
    }
    write.csv(Outlier_plot,paste0("./",genes,"_",group,"_Outlier_plot_2.csv"),quote = FALSE)
    write.csv(DNA_Union,paste0("./",genes,"_",group,"_DNA_Union.csv"),quote = FALSE)
    write.csv(RNA_Union,paste0("./",genes,"_",group,"_RNA_Union.csv"),quote = FALSE)
  }
  
  #Outliers summary
  {
    Outlier_summary <- read.csv("Early_stage_All_genes_outliers_summary.csv",header = TRUE,row.names = NULL,stringsAsFactors = FALSE)
    total_sample_number <- unique(Outlier_summary$Total_sample_number)
    #total_sample_number <- 31
    Outlier_summary$ensembl <- as.character(lapply(strsplit(as.character(Outlier_summary$Gene),"| ",fixed = TRUE),function(x) x[2]))
    
    Outlier_summary <- Outlier_summary %>% mutate(Gene_symbol = ifelse(Gene_symbol %in% "", ensembl, Gene_symbol))
    
    Outlier_summary
    #top30 RNA and DNA altered genes
    {
    top30_RNA <- head(Outlier_summary[order(Outlier_summary$RNA.detected.ratio,decreasing = TRUE),],30)
    top30_DNA <- head(Outlier_summary[order(Outlier_summary$DNA.detected.ratio,decreasing = TRUE),],30)
    
    top30 <- head(Outlier_summary[order(Outlier_summary$OutlierN_all,decreasing = TRUE),],30)
    
    #Top30 RNA+DNA
    {
      top30$Gene_symbol <- factor(top30$Gene_symbol,levels = as.factor(head(Outlier_summary[order(Outlier_summary$OutlierN_all,decreasing = TRUE),]$Gene_symbol,30)))
      Top_altered_RNA <- ggplot(top30,aes(x=top30$Gene_symbol,y=1,fill = top30$RNA.detected.ratio))+geom_bar(stat = "identity")+
        coord_flip()+
        scale_fill_gradient2(low="light blue",high="dark red",mid = "red",midpoint = 0.9,limits = c(0,1))+
        theme_bw()+
        theme(
          plot.margin = unit(x=c(15,0,10,0),units="pt"),
          legend.position="none",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face="bold", color="black", size=18),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(color="black",family = "Arial", size=10),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      Top_altered_DNA <- ggplot(top30,aes(x=top30$Gene_symbol,y=1,fill = top30$DNA.detected.ratio))+geom_bar(stat = "identity")+
        coord_flip()+theme_bw()+
        scale_fill_gradient2(low="light blue",high="dark red",mid = "red",midpoint = 0.9,limits = c(0,1))+
        theme(
          plot.margin = unit(x=c(15,0,10,0),units="pt"),
          legend.position="left",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face="bold", color="black", size=18),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(color="black",family = "Arial", size=10),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color="black", size=10, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      p <- ggarrange(Top_altered_DNA,Top_altered_RNA,align = "h",widths = c(7,1))
      ggsave(p,filename = "Top_altered_genes_20220918.pdf",width = 2.25,height = 4.01,device = "pdf")
    }
    
    #Top30 RNA
    {
    top30_RNA$Gene_symbol <- factor(top30_RNA$Gene_symbol,levels = as.factor(head(Outlier_summary[order(Outlier_summary$RNA.detected.ratio,decreasing = TRUE),]$Gene_symbol,30)))
    Top_altered_RNA <- ggplot(top30_RNA,aes(x=top30_RNA$Gene_symbol,y=1,fill = top30_RNA$RNA.detected.ratio))+geom_bar(stat = "identity")+
      coord_flip()+
      scale_fill_gradient2(low="light blue",high="dark red",mid = "red",midpoint = 0.9,limits = c(0,1))+
      theme_bw()+
      theme(
        plot.margin = unit(x=c(15,0,10,0),units="pt"),
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold", color="black", size=18),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        legend.text= element_text(color="black",family = "Arial", size=10),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    Top_altered_DNA <- ggplot(top30_RNA,aes(x=top30_RNA$Gene_symbol,y=1,fill = top30_RNA$DNA.detected.ratio))+geom_bar(stat = "identity")+
      coord_flip()+theme_bw()+
      scale_fill_gradient2(low="light blue",high="dark red",mid = "red",midpoint = 0.9,limits = c(0,1))+
      theme(
        plot.margin = unit(x=c(15,0,10,0),units="pt"),
        legend.position="left",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold", color="black", size=18),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        legend.text= element_text(color="black",family = "Arial", size=10),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=10, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    p <- ggarrange(Top_altered_DNA,Top_altered_RNA,align = "h",widths = c(7,1))
    ggsave(p,filename = "Top_altered_RNA_20220918.pdf",width = 2.25,height = 4.01,device = "pdf")
    }
    
    #Top30 DNA
    {
      top30_DNA$Gene_symbol <- factor(top30_DNA$Gene_symbol,levels = as.factor(head(Outlier_summary[order(Outlier_summary$DNA.detected.ratio,decreasing = TRUE),]$Gene_symbol,30)))
      Top_altered_DNA <- ggplot(top30_DNA,aes(x=top30_DNA$Gene_symbol,y=1,fill = top30_DNA$DNA.detected.ratio))+geom_bar(stat = "identity")+
        coord_flip()+
        scale_fill_gradient2(low="light blue",high="dark red",mid = "red",midpoint = 0.9,limits = c(0,1))+
        theme_bw()+
        theme(
          plot.margin = unit(x=c(15,0,10,0),units="pt"),
          legend.position="left",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face="bold", color="black", size=18),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(color="black",family = "Arial", size=10),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color="black", size=10, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      Top_altered_RNA <- ggplot(top30_DNA,aes(x=top30_DNA$Gene_symbol,y=1,fill = top30_DNA$RNA.detected.ratio))+geom_bar(stat = "identity")+
        coord_flip()+theme_bw()+
        scale_fill_gradient2(low="light blue",high="dark red",mid = "red",midpoint = 0.9,limits = c(0,1))+
        theme(
          plot.margin = unit(x=c(15,0,10,0),units="pt"),
          legend.position="none",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face="bold", color="black", size=18),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(color="black",family = "Arial", size=10),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      p <- ggarrange(Top_altered_DNA,Top_altered_RNA,align = "h",widths = c(7.5,1))
      ggsave(p,filename = "Top_altered_DNA_20220918.pdf",width = 2.81,height = 4.01,device = "pdf")
    }
    
    }
    #summary
    {
    Outlier_summary_forplot <- as.data.frame(array(dim=c(10,4)))
    colnames(Outlier_summary_forplot) <- c("Mean","Sd","Gene_group","Molecular_type")
    
    Outlier_summary_forplot[1,"Mean"] <- mean(Outlier_summary$DNA.unique/total_sample_number)
    Outlier_summary_forplot[1,"Sd"] <- sd(Outlier_summary$DNA.unique/total_sample_number)
    Outlier_summary_forplot[1,"Gene_group"] <- "Background"
    Outlier_summary_forplot[1,"Molecular_type"] <- "DNA"
    
    Outlier_summary_forplot[2,"Mean"] <- mean(Outlier_summary[Outlier_summary$OutlierN_all>0.75*total_sample_number,]$DNA.unique/total_sample_number)
    Outlier_summary_forplot[2,"Sd"] <- sd(Outlier_summary[Outlier_summary$OutlierN_all>0.75*total_sample_number,]$DNA.unique/total_sample_number)
    Outlier_summary_forplot[2,"Gene_group"] <- "Severely altered genes"
    Outlier_summary_forplot[2,"Molecular_type"] <- "DNA"
    
    Outlier_summary_forplot[3,"Mean"] <- mean(Outlier_summary[Outlier_summary$OutlierN_all<=0.75*total_sample_number & Outlier_summary$OutlierN_all>=0.25*total_sample_number,]$DNA.unique/total_sample_number)
    Outlier_summary_forplot[3,"Sd"] <- sd(Outlier_summary[Outlier_summary$OutlierN_all<=0.75*total_sample_number & Outlier_summary$OutlierN_all>=0.25*total_sample_number,]$DNA.unique/total_sample_number)
    Outlier_summary_forplot[3,"Gene_group"] <- "Moderately altered genes"
    Outlier_summary_forplot[3,"Molecular_type"] <- "DNA"
    
    Outlier_summary_forplot[4,"Mean"] <- mean(Outlier_summary[Outlier_summary$OutlierN_all<=0.25*total_sample_number,]$DNA.unique/total_sample_number)
    Outlier_summary_forplot[4,"Sd"] <- sd(Outlier_summary[Outlier_summary$OutlierN_all<=0.25*total_sample_number,]$DNA.unique/total_sample_number)
    Outlier_summary_forplot[4,"Gene_group"] <- "Minimal altered genes"
    Outlier_summary_forplot[4,"Molecular_type"] <- "DNA"
    
    Outlier_summary_forplot[5,"Mean"] <- mean(as.numeric(Outlier_summary$RNA.unique/total_sample_number))
    Outlier_summary_forplot[5,"Sd"] <- sd(Outlier_summary$RNA.unique/total_sample_number)
    Outlier_summary_forplot[5,"Gene_group"] <- "Background"
    Outlier_summary_forplot[5,"Molecular_type"] <- "RNA"
    
    Outlier_summary_forplot[6,"Mean"] <- mean(Outlier_summary[Outlier_summary$OutlierN_all>0.75*total_sample_number,]$RNA.unique/total_sample_number)
    Outlier_summary_forplot[6,"Sd"] <- sd(Outlier_summary[Outlier_summary$OutlierN_all>0.75*total_sample_number,]$RNA.unique/total_sample_number)
    Outlier_summary_forplot[6,"Gene_group"] <- "Severely altered genes"
    Outlier_summary_forplot[6,"Molecular_type"] <- "RNA"
    
    Outlier_summary_forplot[7,"Mean"] <- mean(Outlier_summary[Outlier_summary$OutlierN_all<=0.75*total_sample_number & Outlier_summary$OutlierN_all>=0.25*total_sample_number,]$RNA.unique/total_sample_number)
    Outlier_summary_forplot[7,"Sd"] <- sd(Outlier_summary[Outlier_summary$OutlierN_all<=0.75*total_sample_number & Outlier_summary$OutlierN_all>=0.25*total_sample_number,]$RNA.unique/total_sample_number)
    Outlier_summary_forplot[7,"Gene_group"] <- "Moderately altered genes"
    Outlier_summary_forplot[7,"Molecular_type"] <- "RNA"
    
    Outlier_summary_forplot[8,"Mean"] <- mean(Outlier_summary[Outlier_summary$OutlierN_all<=0.25*total_sample_number,]$RNA.unique/total_sample_number)
    Outlier_summary_forplot[8,"Sd"] <- sd(Outlier_summary[Outlier_summary$OutlierN_all<=0.25*total_sample_number,]$RNA.unique/total_sample_number)
    Outlier_summary_forplot[8,"Gene_group"] <- "Minimal altered genes"
    Outlier_summary_forplot[8,"Molecular_type"] <- "RNA"
    
    Outlier_summary_forplot[9,"Mean"] <- mean(Outlier_summary[Outlier_summary$OutlierN_all>0.9*total_sample_number,]$DNA.unique/total_sample_number)
    Outlier_summary_forplot[9,"Sd"] <- sd(Outlier_summary[Outlier_summary$OutlierN_all>0.9*total_sample_number,]$DNA.unique/total_sample_number)
    Outlier_summary_forplot[9,"Gene_group"] <- "Most severely altered genes"
    Outlier_summary_forplot[9,"Molecular_type"] <- "DNA"
    
    Outlier_summary_forplot[10,"Mean"] <- mean(Outlier_summary[Outlier_summary$OutlierN_all>0.9*total_sample_number,]$RNA.unique/total_sample_number)
    Outlier_summary_forplot[10,"Sd"] <- sd(Outlier_summary[Outlier_summary$OutlierN_all>0.9*total_sample_number,]$RNA.unique/total_sample_number)
    Outlier_summary_forplot[10,"Gene_group"] <- "Most severely altered genes"
    Outlier_summary_forplot[10,"Molecular_type"] <- "RNA"
    
    Outlier_summary_forplot$`Mean-Sd` <- Outlier_summary_forplot$Mean - Outlier_summary_forplot$Sd
    Outlier_summary_forplot$`Mean+Sd` <- Outlier_summary_forplot$Mean + Outlier_summary_forplot$Sd
    
    Outlier_summary_forplot[Outlier_summary_forplot$`Mean-Sd`<0,"Mean-Sd"] <- 0
    }
    
    Outlier_summary_forplot$Gene_group <- factor(Outlier_summary_forplot$Gene_group,levels=c("Background","Minimal altered genes","Moderately altered genes","Severely altered genes","Most severely altered genes"))
    
    Outlier_summary_forplot_all_stage <- Outlier_summary_forplot
    Outlier_summary_forplot_all_stage$Molecular_type <- paste0(Outlier_summary_forplot_all_stage$Molecular_type,"_all_stage")
    Outlier_summary_forplot_early_stage <- Outlier_summary_forplot
    Outlier_summary_forplot_early_stage$Molecular_type <- paste0(Outlier_summary_forplot_early_stage$Molecular_type,"_early_stage")
    Outlier_summary_forplot <- rbind(Outlier_summary_forplot_early_stage,Outlier_summary_forplot_all_stage)
    
    Outlier_summary_forplot$Molecular_type <- factor(Outlier_summary_forplot$Molecular_type,levels=c("DNA_early_stage","","","",""))
    ggplot(Outlier_summary_forplot,aes(x=Gene_group,y=Mean,fill=Molecular_type))+
      geom_bar(stat = "identity",colour = "transparent",position = "dodge")+
      #geom_text(aes(x=variable,y=Mean+1,label=round(Mean,digits = 0)),size = 4,angle = 0,position = "dodge")+
      geom_errorbar(aes(x=Gene_group, ymin=`Mean-Sd`, ymax=`Mean+Sd`), width=0.9, colour="black", alpha=0.5, size=0.5,position = "dodge")+
      xlab("")+
      labs(title = "")+
      scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels = c("0","","0.2","","0.4","","0.6",""),expand = c(0,0),limits = c(0,0.7))+
      #scale_fill_aaas()+
      scale_fill_manual(values = c("DNA_all_stage"="light grey","DNA_early_stage"=alpha("light grey",0.7),"RNA_all_stage"="black","RNA_early_stage"=alpha("black",0.7)))+
      #geom_vline(xintercept=1.5,linetype = 2)+
      #geom_vline(xintercept=3.5,linetype = 2)+
      xlab("")+
      ylab("Proportion of \ndetected patients")+
      theme_bw()+
      theme(
        plot.margin = unit(x=c(15,5,-10,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold", color="black", size=18),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0,size=24,face="bold",vjust=0),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_text(color="black", size=12, angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(color="black", size=12, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(color="black", size=16))
  }
  
  #Severely altered genes function
  {
    ensembel_genes <- as.character(lapply(strsplit(as.character(Outlier_summary[Outlier_summary$OutlierN_all>0.75*53,]$Gene),"| ",fixed = TRUE),function(x) x[2]))
    
    output_KEGG <- "./Severely_altered_genes_function_KEGG.txt"
    output_GO_BP <- "./Severely_altered_genes_function_BP.txt"
    output_GO_MF <- "./Severely_altered_genes_function_MF.txt"
    output_GO_CC <- "./Severely_altered_genes_function_CC.txt"
    KEGG_GO(ensembel_genes,output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC)
    
    severely_altered_genes_KEGG <- read.table("./Severely_altered_genes_function_KEGG.txt",sep = "\t",check.names = FALSE,row.names = 1)
    severely_altered_genes_KEGG[severely_altered_genes_KEGG$p.adjust < 0.05,]$Description
    selected_pathways<-c("Colorectal cancer","Transcriptional misregulation in cancer","PD-L1 expression and PD-1 checkpoint pathway in cancer","B cell receptor signaling pathway","T cell receptor signaling pathway","C-type lectin receptor signaling pathway","Natural killer cell mediated cytotoxicity","NOD-like receptor signaling pathway")
    severely_altered_genes_KEGG_selected <- severely_altered_genes_KEGG[severely_altered_genes_KEGG$Description%in%selected_pathways,]
    
    selected_pathways <- gsub("PD-L1 expression and PD-1 checkpoint pathway in cancer","PD-L1 expression and PD-1 checkpoint pathway",selected_pathways)
    severely_altered_genes_KEGG_selected$Description <- gsub("PD-L1 expression and PD-1 checkpoint pathway in cancer","PD-L1 expression and PD-1 checkpoint pathway",severely_altered_genes_KEGG_selected$Description)
    
    severely_altered_genes_KEGG_selected$Description <- factor(severely_altered_genes_KEGG_selected$Description,levels = rev(selected_pathways))
    p <- ggplot(severely_altered_genes_KEGG_selected,aes(x=Description,y=Count,fill=p.adjust))+
      geom_bar(stat = "identity",width = 0.8)+coord_flip()+
      theme_bw()+
      #scale_x_discrete(limits = c(0,8))+
      scale_y_continuous(breaks = c(0,10,20,30,40),labels = c("0","10","20","30","40"),limits = c(0,40),expand = c(0,0))+
      theme(#legend.position="right",
        plot.margin = unit(x=c(5,5,5,5),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.title = element_text(color="black",family = "Arial", size=12),
        legend.text= element_text(color="black",family = "Arial", size=12),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        #axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_line(size = 0.5,color = "black"), 
        axis.line.x = element_line(size = 0.5,color = "black"), 
        axis.text.y = element_text(color="black", size=12, angle = 0,hjust = 1,vjust = 0.5),
        axis.text.x = element_text(color="black", size=12, angle = 0,hjust = 0,vjust = 0),
        #axis.text.y = element_text(face="bold",  color="black", size=18, angle = 90,hjust=0.5),
        axis.title.x = element_text(face="bold", color="black", size=0),
        axis.title.y = element_text(face="bold",color="black", size=0))
    ggsave(p,filename = "./severely_altered_genes_function_KEGG.pdf",width = 6.52,height = 2.99,device = "pdf")
  }
  
  #Severely altered RNA function
  {
    ensembel_genes <- as.character(lapply(strsplit(as.character(Outlier_summary[Outlier_summary$OutlierN_RNA>0.75*53,]$Gene),"| ",fixed = TRUE),function(x) x[2]))
    
    output_KEGG <- "./Severely_altered_RNA_function_KEGG.txt"
    output_GO_BP <- "./Severely_altered_RNA_function_BP.txt"
    output_GO_MF <- "./Severely_altered_RNA_function_MF.txt"
    output_GO_CC <- "./Severely_altered_RNA_function_CC.txt"
    KEGG_GO(ensembel_genes,output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC)
  }
  
  #Severely altered DNA function
  {
    ensembel_genes <- as.character(lapply(strsplit(as.character(Outlier_summary[Outlier_summary$OutlierN_DNA>0.75*53,]$Gene),"| ",fixed = TRUE),function(x) x[2]))
    
    output_KEGG <- "./Severely_altered_DNA_function_KEGG.txt"
    output_GO_BP <- "./Severely_altered_DNA_function_BP.txt"
    output_GO_MF <- "./Severely_altered_DNA_function_MF.txt"
    output_GO_CC <- "./Severely_altered_DNA_function_CC.txt"
    KEGG_GO(ensembel_genes,output_KEGG,output_GO_BP,output_GO_MF,output_GO_CC)
  }
}
}

#Figure 5
{
  #tissue-paired plasma correlation
  {
    #KEGG pathway genes
    PATH_ID_NAME <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/09.Published_data/pathway/PATH_ID_NAME_modified.csv",header = TRUE,row.names = 1)
    #multiomics plasma
    counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Plasma_TPM.txt",sep="\t",header = T, row.names = 1)
    #multiomics tissue
    counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Tissue_TPM.txt",sep="\t",header = T, row.names = 1)
    
    counts <- log2(as.matrix(counts)+1)
    
    # set progress
    pb <- progress_bar$new(
      format = "  Processing [:bar] :percent eta: :eta",
      total = length(unique(PATH_ID_NAME$DESCRPTION)), clear = FALSE, width= 60) 
    #get all genes in each pathway in KEGG
    n=1
    pathway_count={}
    while(n <= length(unique(PATH_ID_NAME$DESCRPTION))){
      pathway_name <- as.character(unique(PATH_ID_NAME$DESCRPTION)[n])
      pathway <- subset(PATH_ID_NAME,DESCRPTION==pathway_name)
      
      candidate <- pathway$ensembl_gene_id #for colnames with ensemble id
      #candidate <- pathway$hgnc_symbol #for colnames with gene symbol
      
      candidate <- unique(candidate)
      
      #summary one pathway total count in each sample and make 1 row dataframe
      j=1
      pathway_gene_count={}
      while(j<=length(candidate)){
        target <- candidate[j]
        if(length(grep(target,rownames(counts)))==0) {
          #print(paste0("No ",target," in this dataset."))
          j=j+1
        } else {
          #temp <- counts[which(rownames(counts)==target),]
          temp <- counts[grep(target,rownames(counts),fixed=TRUE),]
          pathway_gene_count <- rbind(pathway_gene_count,temp)
          j=j+1
        }
      }
      if(is.null(pathway_gene_count)){
        n=n+1
        next;
      } else {
        one_pathway <- as.data.frame(t(colMeans(pathway_gene_count))) #colSums
        rownames(one_pathway) <- pathway_name
        pathway_count <- rbind(pathway_count,one_pathway)
        n=n+1
      }
      pb$tick()
      Sys.sleep(1 / 100)
      
    }
    write.table(pathway_count,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Tissue_pathway_mean.txt",sep = "\t", quote = FALSE)
    
    
    Plasma <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/Plasma_pathway_mean.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
    Tissue <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/Tissue_pathway_mean.txt",sep = "\t", row.names = 1, header = TRUE,check.names = FALSE)
    
    gene_ids <- rownames(Plasma)
    
    sample_ids <- c("CRC-PKU-27","CRC-PKU-28","CRC-PKU-29","CRC-PKU-30","CRC-PKU-32","CRC-PKU-34",
                    "CRC-PKU-35","CRC-PKU-36","CRC-PKU-37","CRC-PKU-38","CRC-PKU-39","CRC-PKU-40","CRC-PKU-41")
    
    Tissue_T_forcor <- Tissue[gene_ids,paste0(sample_ids,"-T")]
    colnames(Tissue_T_forcor) <- gsub("-T","",colnames(Tissue_T_forcor))
    Tissue_N_forcor <- Tissue[gene_ids,paste0(sample_ids,"-N")]
    colnames(Tissue_N_forcor) <- gsub("-N","",colnames(Tissue_N_forcor))
    Plasma_forcor <- Plasma[gene_ids,paste0(sample_ids,"-pico")]
    colnames(Plasma_forcor) <- gsub("-pico","",colnames(Plasma_forcor))
    
    forcor1 <- Plasma_forcor
    forcor2 <- Tissue_T_forcor
    cor_methods <- "spearman"
    output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/tissue_paired_analysis/Spearman_correlation/Plasma_vs_Tissue_pathway_mean/"
    dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/tissue_paired_analysis/Spearman_correlation/Plasma_vs_Tissue_pathway_mean/",recursive = TRUE)
    paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,1)
    
    forcor1 <- Plasma_forcor
    forcor2 <- Tissue_N_forcor
    cor_methods <- "spearman"
    output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/tissue_paired_analysis/Spearman_correlation/Plasma_vs_TumorAdjacentNormal_pathway_mean/"
    dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/tissue_paired_analysis/Spearman_correlation/Plasma_vs_TumorAdjacentNormal_pathway_mean/",recursive = TRUE)
    paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,1)
    
    forcor1 <- (Plasma_forcor+1)/(rowSums(Plasma[,-grep("CRC|STAD|NC-PKU-mix.-pico",colnames(Plasma))])+1)
    forcor2 <- (Tissue_T_forcor+1)/(Tissue_N_forcor+1)
    cor_methods <- "spearman"
    output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/tissue_paired_analysis/Spearman_correlation/Plasma_vs_Tumor_TvsN_pathway_mean/"
    dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/tissue_paired_analysis/Spearman_correlation/Plasma_vs_Tumor_TvsN_pathway_mean/",recursive = TRUE)
    paired_correlation(forcor1,forcor2,sample_ids,gene_ids,cor_methods,output_dir,1)
  }
  
  #LM22 radar plot:immune cell fraction of tissue and plasma
  {
    #plasma
    {
    {
      composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/res_cibersort_multiomics_plasma_PBMC_tissue_LM22_abs_20220608_perm1000.txt",sep="\t",header=TRUE,row.names=1)
      composition <- composition[,-which(colnames(composition)=="RMSE")]
      composition <- composition[,-which(colnames(composition)=="Correlation")]
      composition <- composition[,-which(colnames(composition)=="P.value")]
      composition <- composition[,-which(colnames(composition)=="Absolute.score..sig.score.")]
      
      composition$B.cells <- composition$B.cells.naive+composition$B.cells.memory
      composition <- composition[,-which(colnames(composition)=="B.cells.naive")]
      composition <- composition[,-which(colnames(composition)=="B.cells.memory")]
      composition$T.cells.CD4 <- composition$T.cells.CD4.naive+composition$T.cells.CD4.memory.resting+composition$T.cells.CD4.memory.activated+composition$T.cells.follicular.helper+composition$T.cells.regulatory..Tregs.
      composition <- composition[,-which(colnames(composition)=="T.cells.CD4.naive")]
      composition <- composition[,-which(colnames(composition)=="T.cells.CD4.memory.resting")]
      composition <- composition[,-which(colnames(composition)=="T.cells.CD4.memory.activated")]
      composition <- composition[,-which(colnames(composition)=="T.cells.follicular.helper")]
      composition <- composition[,-which(colnames(composition)=="T.cells.regulatory..Tregs.")]
      #composition$NK.cells <- composition$NK.cells.resting+composition$NK.cells.activated
      #composition <- composition[,-which(colnames(composition)=="NK.cells.resting")]
      #composition <- composition[,-which(colnames(composition)=="NK.cells.activated")]
      composition$Macrophages <- composition$Macrophages.M0+composition$Macrophages.M1+composition$Macrophages.M2
      composition <- composition[,-which(colnames(composition)=="Macrophages.M0")]
      composition <- composition[,-which(colnames(composition)=="Macrophages.M1")]
      composition <- composition[,-which(colnames(composition)=="Macrophages.M2")]
      composition$Dendritic.cells <- composition$Dendritic.cells.resting+composition$Dendritic.cells.activated
      composition <- composition[,-which(colnames(composition)=="Dendritic.cells.resting")]
      composition <- composition[,-which(colnames(composition)=="Dendritic.cells.activated")]
      composition$Mast.cells <- composition$Mast.cells.resting+composition$Mast.cells.activated
      composition <- composition[,-which(colnames(composition)=="Mast.cells.resting")]
      composition <- composition[,-which(colnames(composition)=="Mast.cells.activated")]
    }
    
    #composition <- as.data.frame(t(apply(composition, 1, function(x) x / sum(as.numeric(x)) * 10^2)))
    composition <- na.omit(composition)
    composition<-composition[grep("pico",rownames(composition)),]
    composition<-composition[-grep("mix..pico",rownames(composition)),]
    composition<-composition[-grep("CRC-PKU-5-pico",rownames(composition)),]
    composition<-composition[-grep("NC-PKU-mix17-pico",rownames(composition)),]
    composition<-composition[-grep("STAD-PKU-4-pico",rownames(composition)),]
    composition<-composition[-grep("STAD",rownames(composition)),]
    #paired tissue and plasma
    composition<-composition[grep("NC|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
    
    composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[1]))
    composition$group<-gsub("NC","HD",composition$group)
    
    for_wilcox <- composition[,-which(colnames(composition)=="group")]
    positive_ids <- rownames(composition[composition$group=="CRC",])
    negative_ids <- rownames(composition[composition$group=="HD",])
    x=1
    test_result <- data.frame(matrix(nrow=ncol(for_wilcox),ncol = 3))
    colnames(test_result) <- c("Cell type","pvalue","log2(fold change)")
    while(x<=ncol(for_wilcox)){
      cell_type <- colnames(for_wilcox)[x]
      message(cell_type)
      FC <- log2(mean(for_wilcox[positive_ids,cell_type]) / mean(for_wilcox[negative_ids,cell_type]))
      result <- wilcox.test(for_wilcox[positive_ids,cell_type],for_wilcox[negative_ids,cell_type])
      
      test_result[x,"Cell type"] <- cell_type
      test_result[x,"pvalue"] <- result$p.value
      test_result[x,"log2(fold change)"] <- FC
      
      x=x+1
    }
    write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],"./Figure 5/LM22/LM22_Plasma_test.csv")
    
    radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=median)
    radar <- as.tibble(radar)
    o <- rev(order(radar[2,-1]))+1
    #x <- c(8,3,10,5,9,6,11,13,7,12,4,2)
    radar <- radar[,c(1,o)]
    
    colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
    colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
    colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
    colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
    colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
    #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
    p <- ggradar(radar[c(2,1),],grid.min = 0,grid.mid = 0.25, grid.max = 0.5,
            #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
            #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
            grid.line.width = 2,
            base.size = 58,
            values.radar = c("", "0.25", "0.5"),
            plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
            font.radar = "Arial",
            group.point.size = 7,
            group.line.width = 4,
            legend.position = "right",
            background.circle.transparency = 0.0,
            group.colours = c(CRC="#8B5F65",HD="#87CEEB",Tumor="#FCB514",Normal="blue")) #STAD="red",HCC="red",LUAD="red",ESCA="red",HD="#87CEEB",GIC="#FCB514",Tumor="#FCB514",Normal="blue"
    ggsave(plot=p,filename = "./Figure 5/LM22/Plasma_radar.pdf",device = "pdf",width = 15,height = 12)
    }
    
    #tissue
    {
      {
        composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/res_cibersort_multiomics_plasma_PBMC_tissue_LM22_abs_20220608_perm1000.txt",sep="\t",header=TRUE,row.names=1)
        composition <- composition[,-which(colnames(composition)=="RMSE")]
        composition <- composition[,-which(colnames(composition)=="Correlation")]
        composition <- composition[,-which(colnames(composition)=="P.value")]
        composition <- composition[,-which(colnames(composition)=="Absolute.score..sig.score.")]
        
        composition$B.cells <- composition$B.cells.naive+composition$B.cells.memory
        composition <- composition[,-which(colnames(composition)=="B.cells.naive")]
        composition <- composition[,-which(colnames(composition)=="B.cells.memory")]
        composition$T.cells.CD4 <- composition$T.cells.CD4.naive+composition$T.cells.CD4.memory.resting+composition$T.cells.CD4.memory.activated+composition$T.cells.follicular.helper+composition$T.cells.regulatory..Tregs.
        composition <- composition[,-which(colnames(composition)=="T.cells.CD4.naive")]
        composition <- composition[,-which(colnames(composition)=="T.cells.CD4.memory.resting")]
        composition <- composition[,-which(colnames(composition)=="T.cells.CD4.memory.activated")]
        composition <- composition[,-which(colnames(composition)=="T.cells.follicular.helper")]
        composition <- composition[,-which(colnames(composition)=="T.cells.regulatory..Tregs.")]
        #composition$NK.cells <- composition$NK.cells.resting+composition$NK.cells.activated
        #composition <- composition[,-which(colnames(composition)=="NK.cells.resting")]
        #composition <- composition[,-which(colnames(composition)=="NK.cells.activated")]
        composition$Macrophages <- composition$Macrophages.M0+composition$Macrophages.M1+composition$Macrophages.M2
        composition <- composition[,-which(colnames(composition)=="Macrophages.M0")]
        composition <- composition[,-which(colnames(composition)=="Macrophages.M1")]
        composition <- composition[,-which(colnames(composition)=="Macrophages.M2")]
        composition$Dendritic.cells <- composition$Dendritic.cells.resting+composition$Dendritic.cells.activated
        composition <- composition[,-which(colnames(composition)=="Dendritic.cells.resting")]
        composition <- composition[,-which(colnames(composition)=="Dendritic.cells.activated")]
        composition$Mast.cells <- composition$Mast.cells.resting+composition$Mast.cells.activated
        composition <- composition[,-which(colnames(composition)=="Mast.cells.resting")]
        composition <- composition[,-which(colnames(composition)=="Mast.cells.activated")]
      }
      
      #composition <- as.data.frame(t(apply(composition, 1, function(x) x / sum(as.numeric(x)) * 10^2)))
      composition <- na.omit(composition)
      composition<-composition[-grep("PBMC",rownames(composition)),]
      composition<-composition[-grep("pico",rownames(composition)),]
      
      #paired tissue and plasma
      composition<-composition[grep("NC|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
      
      composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[length(x)]))
      composition$group<-gsub("T","Tumor",composition$group)
      composition$group<-gsub("N","Normal",composition$group)
      
      for_wilcox <- composition[,-which(colnames(composition)=="group")]
      positive_ids <- rownames(composition[composition$group=="Tumor",])
      negative_ids <- rownames(composition[composition$group=="Normal",])
      x=1
      test_result <- data.frame(matrix(nrow=ncol(for_wilcox),ncol = 3))
      colnames(test_result) <- c("Cell type","pvalue","log2(fold change)")
      while(x<=ncol(for_wilcox)){
        cell_type <- colnames(for_wilcox)[x]
        message(cell_type)
        FC <- log2(mean(for_wilcox[positive_ids,cell_type]) / mean(for_wilcox[negative_ids,cell_type]))
        result <- wilcox.test(for_wilcox[positive_ids,cell_type],for_wilcox[negative_ids,cell_type])
        
        test_result[x,"Cell type"] <- cell_type
        test_result[x,"pvalue"] <- result$p.value
        test_result[x,"log2(fold change)"] <- FC
        
        x=x+1
      }
      write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],"./Figure 5/LM22/LM22_Tissue_test.csv")
      
      radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=median)
      radar <- as.tibble(radar)
      #o <- rev(order(radar[2,-1]))+1
      #o <- c(8,3,10,5,9,6,11,13,7,12,4,2)
      radar <- radar[,c(1,o)]
      
      colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
      colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
      colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
      colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
      colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
      #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
      p <- ggradar(radar[c(2,1),],grid.min = 0,grid.mid = 0.35, grid.max = 0.7,
              #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
              #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
              grid.line.width = 2,
              base.size = 58,
              values.radar = c("", "0.35", "0.7"),
              plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
              font.radar = "Arial",
              group.point.size = 7,
              group.line.width = 4,
              legend.position = "right",
              background.circle.transparency = 0.0,
              group.colours = c(CRC="#8B5F65",HD="#87CEEB",Tumor="#FCB514",Normal="blue")) #STAD="red",HCC="red",LUAD="red",ESCA="red",HD="#87CEEB",GIC="#FCB514",Tumor="#FCB514",Normal="blue"
      ggsave(plot=p,filename = "./Figure 5/LM22/Tissue_radar.pdf",device = "pdf",width = 15,height = 12)
    }
  }
  
  #EPIC radar plot: immune + CAFs fraction of tissue and plasma
  {
    #plasma
    {
      composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Plasma_EPIC_cellFraction_TRef.txt",sep="\t",header=TRUE,row.names=1)
      composition<-composition[-grep("CRC-PKU-5-pico",rownames(composition)),]
      composition<-composition[-grep("NC-PKU-mix17-pico",rownames(composition)),]
      composition<-composition[-grep("STAD-PKU-4-pico",rownames(composition)),]
      composition<-composition[-grep("mix..pico",rownames(composition)),]
      composition<-composition[grep("NC|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
      
      composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[1]))
      composition$group<-gsub("NC","HD",composition$group)
      
      for_wilcox <- composition[,-which(colnames(composition)=="group")]
      positive_ids <- rownames(composition[composition$group=="CRC",])
      negative_ids <- rownames(composition[composition$group=="HD",])
      x=1
      test_result <- data.frame(matrix(nrow=ncol(for_wilcox),ncol = 3))
      colnames(test_result) <- c("Cell type","pvalue","log2(fold change)")
      while(x<=ncol(for_wilcox)){
        cell_type <- colnames(for_wilcox)[x]
        message(cell_type)
        FC <- log2(mean(for_wilcox[positive_ids,cell_type]) / mean(for_wilcox[negative_ids,cell_type]))
        result <- wilcox.test(for_wilcox[positive_ids,cell_type],for_wilcox[negative_ids,cell_type])
        
        test_result[x,"Cell type"] <- cell_type
        test_result[x,"pvalue"] <- result$p.value
        test_result[x,"log2(fold change)"] <- FC
        
        x=x+1
      }
      write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],"./Figure 5/EPIC_Plasma_test.csv")
      
      radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=mean)
      radar <- as.tibble(radar)
      o <- rev(order(radar[2,-1]))+1
      #o <- c(8,3,10,5,9,6,11,13,7,12,4,2)
      radar <- radar[,c(1,o)]
      
      colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
      #colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
      #colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
      #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
      ggradar(radar[c(2,1),-c(2)],grid.min = 0,grid.mid = 0.02, grid.max = 0.04,
              #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
              #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
              grid.line.width = 1,
              base.size = 58,
              values.radar = c("", "2%", "4%"),
              plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
              font.radar = "Arial",
              group.point.size = 3,
              group.line.width = 2,
              legend.position = "right",
              background.circle.transparency = 0.0,
              group.colours = c(CRC="#8B5F65",HD="#87CEEB",Tumor="#FCB514",Normal="blue")) #STAD="red",HCC="red",LUAD="red",ESCA="red",HD="#87CEEB",GIC="#FCB514",Tumor="#FCB514",Normal="blue"
    }
    
    #tissue
    {
      composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Tissue_EPIC_cellFraction_TRef.txt",sep="\t",header=TRUE,row.names=1)
      composition<-composition[-grep("2384058",rownames(composition)),]
      composition<-composition[-grep("2399129",rownames(composition)),]
      #composition<-composition[-grep("STAD-PKU-4-pico",rownames(composition)),]
      #composition<-composition[-grep("mix..pico",rownames(composition)),]
      composition<-composition[grep("NC|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
      
      composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[length(x)]))
      composition$group<-gsub("T","Tumor",composition$group)
      composition$group<-gsub("N","Normal",composition$group)
      
      for_wilcox <- composition[,-which(colnames(composition)=="group")]
      positive_ids <- rownames(composition[composition$group=="Tumor",])
      negative_ids <- rownames(composition[composition$group=="Normal",])
      x=1
      test_result <- data.frame(matrix(nrow=ncol(for_wilcox),ncol = 3))
      colnames(test_result) <- c("Cell type","pvalue","log2(fold change)")
      while(x<=ncol(for_wilcox)){
        cell_type <- colnames(for_wilcox)[x]
        message(cell_type)
        FC <- log2(mean(for_wilcox[positive_ids,cell_type]) / mean(for_wilcox[negative_ids,cell_type]))
        result <- wilcox.test(for_wilcox[positive_ids,cell_type],for_wilcox[negative_ids,cell_type])
        
        test_result[x,"Cell type"] <- cell_type
        test_result[x,"pvalue"] <- result$p.value
        test_result[x,"log2(fold change)"] <- FC
        
        x=x+1
      }
      write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],"./Figure 5/EPIC_Tissue_test.csv")
      
      radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=mean)
      radar <- as.tibble(radar)
      #o <- rev(order(radar[2,-1]))+1
      #x <- c(2,6,4,7,5,8,3)
      radar <- radar[,c(1,o)]
      
      colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
      #colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
      #colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
      #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
      ggradar(radar[c(2,1),-c(2)],grid.min = 0,grid.mid = 0.075, grid.max = 0.15,
              #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
              #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
              grid.line.width = 1,
              base.size = 58,
              values.radar = c("", "7.5%", "15%"),
              plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
              font.radar = "Arial",
              group.point.size = 3,
              group.line.width = 2,
              legend.position = "right",
              background.circle.transparency = 0.0,
              group.colours = c(CRC="#8B5F65",HD="#87CEEB",Tumor="#FCB514",Normal="blue")) #STAD="red",HCC="red",LUAD="red",ESCA="red",HD="#87CEEB",GIC="#FCB514",Tumor="#FCB514",Normal="blue"
    }
  }
  
  #Colorectal scRNA seq plot: immune and CAFs fraction of tissue and plasma
  {
    signature_gene_number <- 10
    #Plasma
    {
    #composition <-read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/20220913_Multiomics_CRC_scRNAseq_CD4_CD8_CAFs_detailed_top",signature_gene_number,".txt"),sep="\t",header=TRUE,row.names=1)
    composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/20220907_Multiomics_abs_CRC_scRNAseq_cell_type.txt",sep="\t",header=TRUE,row.names=1)
    composition <- composition[,-which(colnames(composition)=="RMSE")]
    composition <- composition[,-which(colnames(composition)=="Correlation")]
    composition <- composition[,-which(colnames(composition)=="P.value")]
    composition <- composition[,-which(colnames(composition)=="Absolute.score..sig.score.")]
    #composition <- composition[,-which(colnames(composition)=="NK.cells")]
    
    composition<-composition[grep("pico",rownames(composition)),]
    composition<-composition[-grep("mix..pico",rownames(composition)),]
    composition<-composition[-grep("CRC-PKU-5-pico",rownames(composition)),]
    composition<-composition[-grep("NC-PKU-mix17-pico",rownames(composition)),]
    composition<-composition[-grep("STAD-PKU-4-pico",rownames(composition)),]
    composition<-composition[-grep("STAD",rownames(composition)),]
    #paired tissue and plasma
    composition<-composition[grep("NC|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
    
    composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[1]))
    
    composition$group<-gsub("NC","HD",composition$group)
    
    for_wilcox <- composition[,-which(colnames(composition)=="group")]
    positive_ids <- rownames(composition[composition$group=="CRC",])
    negative_ids <- rownames(composition[composition$group=="HD",])
    x=1
    test_result <- data.frame(matrix(nrow=ncol(for_wilcox),ncol = 3))
    colnames(test_result) <- c("Cell type","pvalue","log2(fold change)")
    while(x<=ncol(for_wilcox)){
    cell_type <- colnames(for_wilcox)[x]
    message(cell_type)
    FC <- log2(mean(for_wilcox[positive_ids,cell_type]) / mean(for_wilcox[negative_ids,cell_type]))
    result <- wilcox.test(for_wilcox[positive_ids,cell_type],for_wilcox[negative_ids,cell_type])
    
    test_result[x,"Cell type"] <- cell_type
    test_result[x,"pvalue"] <- result$p.value
    test_result[x,"log2(fold change)"] <- FC
    
    x=x+1
    }
    #write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_detailed_top",signature_gene_number,"_Plasma_test.csv"))
    write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],"./Figure 5/CRC_scRNAseq_radar/Cell_type_Plasma_test.csv")
    
    radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=median)
    radar <- as.tibble(radar)
    o <- rev(order(radar[2,-1]))+1
    #x <- c(8,3,10,5,9,6,11,13,7,12,4,2)
    radar <- radar[,c(1,o)]
    
    colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
    #colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
    #colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
    #colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
    #colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
    #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
    p <- ggradar(radar[c(2,1),],grid.min = 0,grid.mid = 50, grid.max = 100,
            #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
            #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
            grid.line.width = 1,
            base.size = 58,
            #values.radar = c("", "0.25", "0.5"),
            plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
            font.radar = "Arial",
            group.point.size = 3,
            group.line.width = 2,
            legend.position = "right",
            background.circle.transparency = 0.0,
            group.colours = c(CRC="#8B5F65",HD="#87CEEB",Tumor="#FCB514",Normal="blue")) #STAD="red",HCC="red",LUAD="red",ESCA="red",HD="#87CEEB",GIC="#FCB514",Tumor="#FCB514",Normal="blue"
    #ggsave(plot=p,filename = paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_detailed_top",signature_gene_number,"_Plasma_radar.pdf"),device = "pdf",width = 15,height = 12)
    ggsave(plot=p,filename = "./Figure 5/CRC_scRNAseq_radar/Cell_type_Plasma_radar.pdf",device = "pdf",width = 15,height = 12)
    }
    
    #Tissue
    {
      #composition <-read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/20220913_Multiomics_CRC_scRNAseq_CD4_CD8_CAFs_detailed_top",signature_gene_number,".txt"),sep="\t",header=TRUE,row.names=1)
      composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/20220907_Multiomics_abs_CRC_scRNAseq_cell_type.txt",sep="\t",header=TRUE,row.names=1)
      composition <- composition[,-which(colnames(composition)=="RMSE")]
      composition <- composition[,-which(colnames(composition)=="Correlation")]
      composition <- composition[,-which(colnames(composition)=="P.value")]
      composition <- composition[,-which(colnames(composition)=="Absolute.score..sig.score.")]
      
      composition<-composition[-grep("pico",rownames(composition)),]
      composition<-composition[-grep("PBMC",rownames(composition)),]
      
      #paired tissue and plasma
      composition<-composition[grep("NC|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
      
      composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[length(x)]))
      composition$group<-gsub("T","Tumor",composition$group)
      composition$group<-gsub("N","Normal",composition$group)
      
      for_wilcox <- composition[,-which(colnames(composition)=="group")]
      positive_ids <- rownames(composition[composition$group=="Tumor",])
      negative_ids <- rownames(composition[composition$group=="Normal",])
      x=1
      test_result <- data.frame(matrix(nrow=ncol(for_wilcox),ncol = 3))
      colnames(test_result) <- c("Cell type","pvalue","log2(fold change)")
      while(x<=ncol(for_wilcox)){
        cell_type <- colnames(for_wilcox)[x]
        message(cell_type)
        FC <- log2(mean(for_wilcox[positive_ids,cell_type]) / mean(for_wilcox[negative_ids,cell_type]))
        result <- wilcox.test(for_wilcox[positive_ids,cell_type],for_wilcox[negative_ids,cell_type])
        
        test_result[x,"Cell type"] <- cell_type
        test_result[x,"pvalue"] <- result$p.value
        test_result[x,"log2(fold change)"] <- FC
        
        x=x+1
      }
      #write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_top",signature_gene_number,"_detailed_tissue_test.csv"))
      write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],"./Figure 5/CRC_scRNAseq_radar/Cell_type_tissue_test.csv")
      
      radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=mean)
      radar <- as.tibble(radar)
      #x <- rev(order(radar[2,-1]))+1
      #x <- c(8,3,10,5,9,6,11,13,7,12,4,2)
      radar <- radar[,c(1,o)]
      
      colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
      #colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
      #colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
      #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
      p <- ggradar(radar[c(2,1),],grid.min = 0,grid.mid = 50, grid.max = 100,
              #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
              #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
              grid.line.width = 1,
              base.size = 58,
              #values.radar = c("", "0.25", "0.5"),
              plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
              font.radar = "Arial",
              group.point.size = 3,
              group.line.width = 2,
              legend.position = "right",
              background.circle.transparency = 0.0,
              group.colours = c(CRC="#8B5F65",HD="#87CEEB",Tumor="#FCB514",Normal="blue")) #STAD="red",HCC="red",LUAD="red",ESCA="red",HD="#87CEEB",GIC="#FCB514",Tumor="#FCB514",Normal="blue"
      #ggsave(plot=p,filename = paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_detailed_top",signature_gene_number,"_tissue_radar.pdf"),device = "pdf",width = 15,height = 12)
      ggsave(plot=p,filename = paste0("./Figure 5/CRC_scRNAseq_radar/Cell_type_tissue_radar.pdf"),device = "pdf",width = 15,height = 12)
    }  
  }
  
  #2021-Oncogene
  {
    signature_genes <- read.table("./Figure 5/2021_Oncogene_Single_cell_type_signatures.txt",sep = "\t",check.names = FALSE,header =  TRUE,row.names = 1)
    values <- rownames(signature_genes)
    translated <- AnnotationDbi::select(EnsDb.Hsapiens.v86, key=values,columns=c("SYMBOL","GENEID"),keytype="SYMBOL")
    translated <- translated[grep("ENSG",translated$GENEID),]
    signature_genes$SYMBOL <- rownames(signature_genes)
    signature_genes <- left_join(translated,signature_genes,by=c("SYMBOL"="SYMBOL"))
    signature_genes <- signature_genes[,-which(colnames(signature_genes)=="SYMBOL")]
    write.table(signature_genes,"./Figure 5/2021_Oncogene_Single_cell_type_signatures_gene_ensembl.txt",sep = "\t",quote = FALSE,row.names = FALSE)
  }
  
  #2021-Oncogene: bulk radar plot
  {
    #Plasma
    {
      #composition <-read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/20220913_Multiomics_CRC_scRNAseq_CD4_CD8_CAFs_detailed_top",signature_gene_number,".txt"),sep="\t",header=TRUE,row.names=1)
      composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/20220913_Multiomics_Oncogene_bulk_cell_type_abs.txt",sep="\t",header=TRUE,row.names=1)
      composition <- composition[,-which(colnames(composition)=="RMSE")]
      composition <- composition[,-which(colnames(composition)=="Correlation")]
      composition <- composition[,-which(colnames(composition)=="P.value")]
      composition <- composition[,-which(colnames(composition)=="Absolute.score..sig.score.")]
      #composition <- composition[,-which(colnames(composition)=="NK.cells")]
      
      composition<-composition[grep("pico",rownames(composition)),]
      composition<-composition[-grep("mix..pico",rownames(composition)),]
      composition<-composition[-grep("CRC-PKU-5-pico",rownames(composition)),]
      composition<-composition[-grep("NC-PKU-mix17-pico",rownames(composition)),]
      composition<-composition[-grep("STAD-PKU-4-pico",rownames(composition)),]
      composition<-composition[-grep("STAD",rownames(composition)),]
      #paired tissue and plasma
      composition<-composition[grep("NC|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
      
      composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[1]))
      
      composition$group<-gsub("NC","HD",composition$group)
      
      for_wilcox <- composition[,-which(colnames(composition)=="group")]
      positive_ids <- rownames(composition[composition$group=="CRC",])
      negative_ids <- rownames(composition[composition$group=="HD",])
      x=1
      test_result <- data.frame(matrix(nrow=ncol(for_wilcox),ncol = 3))
      colnames(test_result) <- c("Cell type","pvalue","log2(fold change)")
      while(x<=ncol(for_wilcox)){
        cell_type <- colnames(for_wilcox)[x]
        message(cell_type)
        FC <- log2(mean(for_wilcox[positive_ids,cell_type]) / mean(for_wilcox[negative_ids,cell_type]))
        result <- wilcox.test(for_wilcox[positive_ids,cell_type],for_wilcox[negative_ids,cell_type])
        
        test_result[x,"Cell type"] <- cell_type
        test_result[x,"pvalue"] <- result$p.value
        test_result[x,"log2(fold change)"] <- FC
        
        x=x+1
      }
      #write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_detailed_top",signature_gene_number,"_Plasma_test.csv"))
      write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],"./Figure 5/2021_oncogene/Bulk_cell_type_Plasma_test.csv")
      
      radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=median)
      radar <- as.tibble(radar)
      o <- rev(order(radar[2,-1]))+1
      #x <- c(8,3,10,5,9,6,11,13,7,12,4,2)
      radar <- radar[,c(1,o)]
      
      colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
      #colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
      #colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
      #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
      p <- ggradar(radar[c(2,1),],grid.min = 0,grid.mid = 0.6, grid.max = 1.2,
                   #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
                   #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
                   grid.line.width = 1,
                   base.size = 58,
                   #values.radar = c("", "0.25", "0.5"),
                   plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
                   font.radar = "Arial",
                   group.point.size = 3,
                   group.line.width = 2,
                   legend.position = "right",
                   background.circle.transparency = 0.0,
                   group.colours = c(CRC="#8B5F65",HD="#87CEEB",Tumor="#FCB514",Normal="blue")) #STAD="red",HCC="red",LUAD="red",ESCA="red",HD="#87CEEB",GIC="#FCB514",Tumor="#FCB514",Normal="blue"
      #ggsave(plot=p,filename = paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_detailed_top",signature_gene_number,"_Plasma_radar.pdf"),device = "pdf",width = 15,height = 12)
      ggsave(plot=p,filename = "./Figure 5/2021_oncogene/cell_type_Plasma_radar.pdf",device = "pdf",width = 15,height = 12)
    }
    
    #Tissue
    {
      #composition <-read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/20220913_Multiomics_CRC_scRNAseq_CD4_CD8_CAFs_detailed_top",signature_gene_number,".txt"),sep="\t",header=TRUE,row.names=1)
      composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/20220913_Multiomics_Oncogene_bulk_cell_type_abs.txt",sep="\t",header=TRUE,row.names=1)
      composition <- composition[,-which(colnames(composition)=="RMSE")]
      composition <- composition[,-which(colnames(composition)=="Correlation")]
      composition <- composition[,-which(colnames(composition)=="P.value")]
      composition <- composition[,-which(colnames(composition)=="Absolute.score..sig.score.")]
      
      composition<-composition[-grep("pico",rownames(composition)),]
      composition<-composition[-grep("PBMC",rownames(composition)),]
      
      #paired tissue and plasma
      composition<-composition[grep("NC|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
      
      composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[length(x)]))
      composition$group<-gsub("T","Tumor",composition$group)
      composition$group<-gsub("N","Normal",composition$group)
      
      for_wilcox <- composition[,-which(colnames(composition)=="group")]
      positive_ids <- rownames(composition[composition$group=="Tumor",])
      negative_ids <- rownames(composition[composition$group=="Normal",])
      x=1
      test_result <- data.frame(matrix(nrow=ncol(for_wilcox),ncol = 3))
      colnames(test_result) <- c("Cell type","pvalue","log2(fold change)")
      while(x<=ncol(for_wilcox)){
        cell_type <- colnames(for_wilcox)[x]
        message(cell_type)
        FC <- log2(mean(for_wilcox[positive_ids,cell_type]) / mean(for_wilcox[negative_ids,cell_type]))
        result <- wilcox.test(for_wilcox[positive_ids,cell_type],for_wilcox[negative_ids,cell_type])
        
        test_result[x,"Cell type"] <- cell_type
        test_result[x,"pvalue"] <- result$p.value
        test_result[x,"log2(fold change)"] <- FC
        
        x=x+1
      }
      #write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_top",signature_gene_number,"_detailed_tissue_test.csv"))
      write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],"./Figure 5/2021_oncogene/Bulk_cell_type_tissue_test.csv")
      
      radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=mean)
      radar <- as.tibble(radar)
      #x <- rev(order(radar[2,-1]))+1
      #x <- c(8,3,10,5,9,6,11,13,7,12,4,2)
      radar <- radar[,c(1,o)]
      
      colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
      #colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
      #colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
      #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
      p <- ggradar(radar[c(2,1),],grid.min = 0,grid.mid = 0.6, grid.max = 1.2,
                   #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
                   #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
                   grid.line.width = 1,
                   base.size = 58,
                   #values.radar = c("", "0.25", "0.5"),
                   plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
                   font.radar = "Arial",
                   group.point.size = 3,
                   group.line.width = 2,
                   legend.position = "right",
                   background.circle.transparency = 0.0,
                   group.colours = c(CRC="#8B5F65",HD="#87CEEB",Tumor="#FCB514",Normal="blue")) #STAD="red",HCC="red",LUAD="red",ESCA="red",HD="#87CEEB",GIC="#FCB514",Tumor="#FCB514",Normal="blue"
      #ggsave(plot=p,filename = paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_detailed_top",signature_gene_number,"_tissue_radar.pdf"),device = "pdf",width = 15,height = 12)
      ggsave(plot=p,filename = paste0("./Figure 5/2021_oncogene/Bulk_cell_type_tissue_radar.pdf"),device = "pdf",width = 15,height = 12)
    }  
  }
  
  #2021-Oncogene: scRNA radar plot
  {
    #Plasma
    {
      #composition <-read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/20220913_Multiomics_CRC_scRNAseq_CD4_CD8_CAFs_detailed_top",signature_gene_number,".txt"),sep="\t",header=TRUE,row.names=1)
      composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/20220913_Multiomics_Oncogene_single_cell_type_abs.txt",sep="\t",header=TRUE,row.names=1)
      composition <- composition[,-which(colnames(composition)=="RMSE")]
      composition <- composition[,-which(colnames(composition)=="Correlation")]
      composition <- composition[,-which(colnames(composition)=="P.value")]
      composition <- composition[,-which(colnames(composition)=="Absolute.score..sig.score.")]
      #composition <- composition[,-which(colnames(composition)=="NK.cells")]
      
      composition<-composition[grep("pico",rownames(composition)),]
      composition<-composition[-grep("mix..pico",rownames(composition)),]
      composition<-composition[-grep("CRC-PKU-5-pico",rownames(composition)),]
      composition<-composition[-grep("NC-PKU-mix17-pico",rownames(composition)),]
      composition<-composition[-grep("STAD-PKU-4-pico",rownames(composition)),]
      composition<-composition[-grep("STAD",rownames(composition)),]
      #paired tissue and plasma
      composition<-composition[grep("NC|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
      
      composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[1]))
      
      composition$group<-gsub("NC","HD",composition$group)
      
      for_wilcox <- composition[,-which(colnames(composition)=="group")]
      positive_ids <- rownames(composition[composition$group=="CRC",])
      negative_ids <- rownames(composition[composition$group=="HD",])
      x=1
      test_result <- data.frame(matrix(nrow=ncol(for_wilcox),ncol = 3))
      colnames(test_result) <- c("Cell type","pvalue","log2(fold change)")
      while(x<=ncol(for_wilcox)){
        cell_type <- colnames(for_wilcox)[x]
        message(cell_type)
        FC <- log2(mean(for_wilcox[positive_ids,cell_type]) / mean(for_wilcox[negative_ids,cell_type]))
        result <- wilcox.test(for_wilcox[positive_ids,cell_type],for_wilcox[negative_ids,cell_type])
        
        test_result[x,"Cell type"] <- cell_type
        test_result[x,"pvalue"] <- result$p.value
        test_result[x,"log2(fold change)"] <- FC
        
        x=x+1
      }
      #write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_detailed_top",signature_gene_number,"_Plasma_test.csv"))
      write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],"./Figure 5/2021_oncogene/Single_cell_type_Plasma_test.csv")
      
      radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=mean)
      radar <- as.tibble(radar)
      o <- rev(order(radar[2,-1]))+1
      #x <- c(8,3,10,5,9,6,11,13,7,12,4,2)
      radar <- radar[,c(1,o)]
      
      colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
      #colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
      #colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
      #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
      p <- ggradar(radar[c(2,1),-c(9,14)],grid.min = 0,grid.mid = 0.35, grid.max = 0.7,
                   #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
                   #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
                   grid.line.width = 2,
                   base.size = 58,
                   values.radar = c("", "0.35", "0.7"),
                   plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
                   font.radar = "Arial",
                   group.point.size = 7,
                   group.line.width = 4,
                   legend.position = "right",
                   background.circle.transparency = 0.0,
                   group.colours = c(CRC="#8B5F65",HD="#87CEEB",Tumor="#FCB514",Normal="blue")) #STAD="red",HCC="red",LUAD="red",ESCA="red",HD="#87CEEB",GIC="#FCB514",Tumor="#FCB514",Normal="blue"
      #ggsave(plot=p,filename = paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_detailed_top",signature_gene_number,"_Plasma_radar.pdf"),device = "pdf",width = 15,height = 12)
      ggsave(plot=p,filename = "./Figure 5/2021_oncogene/Single_type_Plasma_radar.pdf",device = "pdf",width = 15,height = 12)
    }
    
    #Tissue
    {
      #composition <-read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/20220913_Multiomics_CRC_scRNAseq_CD4_CD8_CAFs_detailed_top",signature_gene_number,".txt"),sep="\t",header=TRUE,row.names=1)
      composition <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/20220913_Multiomics_Oncogene_single_cell_type_abs.txt",sep="\t",header=TRUE,row.names=1)
      composition <- composition[,-which(colnames(composition)=="RMSE")]
      composition <- composition[,-which(colnames(composition)=="Correlation")]
      composition <- composition[,-which(colnames(composition)=="P.value")]
      composition <- composition[,-which(colnames(composition)=="Absolute.score..sig.score.")]
      
      composition<-composition[-grep("pico",rownames(composition)),]
      composition<-composition[-grep("PBMC",rownames(composition)),]
      
      #paired tissue and plasma
      composition<-composition[grep("NC|CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(composition)),]
      
      composition$group <- as.character(lapply(strsplit(as.character(rownames(composition)),"-",fixed = TRUE),function(x) x[length(x)]))
      composition$group<-gsub("T","Tumor",composition$group)
      composition$group<-gsub("N","Normal",composition$group)
      
      for_wilcox <- composition[,-which(colnames(composition)=="group")]
      positive_ids <- rownames(composition[composition$group=="Tumor",])
      negative_ids <- rownames(composition[composition$group=="Normal",])
      x=1
      test_result <- data.frame(matrix(nrow=ncol(for_wilcox),ncol = 3))
      colnames(test_result) <- c("Cell type","pvalue","log2(fold change)")
      while(x<=ncol(for_wilcox)){
        cell_type <- colnames(for_wilcox)[x]
        message(cell_type)
        FC <- log2(mean(for_wilcox[positive_ids,cell_type]) / mean(for_wilcox[negative_ids,cell_type]))
        result <- wilcox.test(for_wilcox[positive_ids,cell_type],for_wilcox[negative_ids,cell_type])
        
        test_result[x,"Cell type"] <- cell_type
        test_result[x,"pvalue"] <- result$p.value
        test_result[x,"log2(fold change)"] <- FC
        
        x=x+1
      }
      #write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_top",signature_gene_number,"_detailed_tissue_test.csv"))
      write.csv(test_result[order(test_result$pvalue,decreasing = FALSE),],"./Figure 5/2021_oncogene/Single_cell_type_tissue_test.csv")
      
      radar <- aggregate(composition[,-which(colnames(composition)=="group")],list(Group=composition$group), FUN=mean)
      radar <- as.tibble(radar)
      #x <- rev(order(radar[2,-1]))+1
      #x <- c(8,3,10,5,9,6,11,13,7,12,4,2)
      radar <- radar[,c(1,o)]
      
      colnames(radar) <- gsub("."," ",fixed=TRUE,colnames(radar))
      #colnames(radar) <- gsub("T cells CD8","CD8 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells CD4","CD4 T cells",colnames(radar))
      #colnames(radar) <- gsub("T cells gamma delta","γδ T cells",colnames(radar))
      #colnames(radar) <- gsub("Dendritic cells","Dendritic\ncells",colnames(radar))
      #radar <- radar[,-which(colnames(radar)=="Megakaryocytes")]
      p <- ggradar(radar[c(2,1),-c(9,14)],grid.min = 0,grid.mid = 0.2, grid.max = 0.4,
                   #axis.labels = c("B cells","Plasma cells","CD8 T cells","CD4 T cells","γδ T cells","NK cells",
                   #                "Monocytes","Macrophages","Dendritic cells","Mast cells","Eosinophils","Neutrophils"),
                   grid.line.width = 2,
                   base.size = 58,
                   values.radar = c("", "0.2", "0.4"),
                   plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
                   font.radar = "Arial",
                   group.point.size = 7,
                   group.line.width = 4,
                   legend.position = "right",
                   background.circle.transparency = 0.0,
                   group.colours = c(CRC="#8B5F65",HD="#87CEEB",Tumor="#FCB514",Normal="blue")) #STAD="red",HCC="red",LUAD="red",ESCA="red",HD="#87CEEB",GIC="#FCB514",Tumor="#FCB514",Normal="blue"
      #ggsave(plot=p,filename = paste0("./Figure 5/CRC_scRNAseq_radar/CAFs_CD4_CD8_detailed_top",signature_gene_number,"_tissue_radar.pdf"),device = "pdf",width = 15,height = 12)
      ggsave(plot=p,filename = paste0("./Figure 5/2021_oncogene/Single_cell_type_tissue_radar.pdf"),device = "pdf",width = 15,height = 12)
    }  
  }
  
  #scatter plot
  {
    scatter1 <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Plasma_EPIC_cellFraction_TRef.txt",sep="\t",header=TRUE,row.names=1)
    scatter2 <-read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Tissue_EPIC_cellFraction_TRef.txt",sep="\t",header=TRUE,row.names=1)
    
    scatter1<-scatter1[grep("CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(scatter1)),]
    scatter2<-scatter2[grep("CRC-PKU-27|CRC-PKU-28|CRC-PKU-29|CRC-PKU-30|CRC-PKU-32|CRC-PKU-34|CRC-PKU-35|CRC-PKU-36|CRC-PKU-37|CRC-PKU-38|CRC-PKU-39|CRC-PKU-40|CRC-PKU-41",rownames(scatter2)),]
    scatter2<-scatter2[grep("T",rownames(scatter2)),]
    
    colnames(scatter1) <- paste0(colnames(scatter1),"_plasma")
    colnames(scatter2) <- paste0(colnames(scatter2),"_tumor")
    
    rownames(scatter1) <- gsub("-pico","",rownames(scatter1))
    rownames(scatter2) <- gsub("-T","",rownames(scatter2))
    
    scatter1$ID <- rownames(scatter1)
    scatter2$ID <- rownames(scatter2)
    
    scatter <- left_join(scatter1,scatter2,by=c("ID"="ID"))
    
    #scatter$CAFs_plasma <- log10(scatter$CAFs_plasma)
    ggscatter(scatter, x = "CAFs_plasma", y = "CAFs_tumor", 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.coeff.args = list(method = "spearman", label.sep = "\n",size = 8))+
      theme(#legend.position="right",
        legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=20),
        axis.text.y = element_text(face="bold",  color="black", size=20),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))
  }
}

#Figure 6
{
  #clinical correlation
  {
    
    #read in immune fractions
    {
      Immune_fraction <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/res_cibersort_multiomics_plasma_PBMC_tissue_LM22_abs_20220608_perm1000.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
      Immune_fraction2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/res_cibersort_GSE174302_LM22_abs_20220614_perm1000.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
      Immune_fraction <- rbind(Immune_fraction,Immune_fraction2)
      
      Immune_fraction <- Immune_fraction[,-grep("P-value|Correlation|RMSE|Absolute score",colnames(Immune_fraction))]
      {
        Immune_fraction$`B cells` <- Immune_fraction$`B cells naive`+Immune_fraction$`B cells memory`
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="B cells naive")]
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="B cells memory")]
        Immune_fraction$`T cells CD4` <- Immune_fraction$`T cells CD4 naive`+Immune_fraction$`T cells CD4 memory resting`+Immune_fraction$`T cells CD4 memory activated`+Immune_fraction$`T cells follicular helper`+Immune_fraction$`T cells regulatory (Tregs)`
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="T cells CD4 naive")]
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="T cells CD4 memory resting")]
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="T cells CD4 memory activated")]
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="T cells follicular helper")]
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="T cells regulatory (Tregs)")]
        #Immune_fraction$`NK cells` <- Immune_fraction$`NK cells resting`+Immune_fraction$`NK cells activated`
        #Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="NK cells resting")]
        #Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="NK cells activated")]
        Immune_fraction$Macrophages <- Immune_fraction$`Macrophages M0`+Immune_fraction$`Macrophages M1`+Immune_fraction$`Macrophages M2`
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Macrophages M0")]
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Macrophages M1")]
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Macrophages M2")]
        Immune_fraction$`Dendritic cells` <- Immune_fraction$`Dendritic cells resting`+Immune_fraction$`Dendritic cells activated`
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Dendritic cells resting")]
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Dendritic cells activated")]
        Immune_fraction$`Mast cells` <- Immune_fraction$`Mast cells resting`+Immune_fraction$`Mast cells activated`
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Mast cells resting")]
        Immune_fraction <- Immune_fraction[,-which(colnames(Immune_fraction)=="Mast cells activated")]
      }
      #Immune_fraction <- Immune_fraction[-grep("CRC.PKU.15|STAD.PKU.23|STAD.PKU.36",rownames(Immune_fraction)),]
      #Immune_fraction <- Immune_fraction[grep("STAD",rownames(Immune_fraction)),]
      Immune_fraction <- Immune_fraction[-which(rowSums(Immune_fraction)==0),] # This step will remove samples that have no immune cells fraction
      Immune_fraction <- Immune_fraction/rowSums(Immune_fraction)
      Immune_fraction$ID <- rownames(Immune_fraction)
      Immune_fraction$ID <- gsub(".","-",fixed = TRUE,Immune_fraction$ID)
      rownames(Immune_fraction) <- gsub(".","-",fixed = TRUE,rownames(Immune_fraction))
    }
    
    rank1 <- Immune_fraction[Immune_fraction$`T cells gamma delta`>0.0,]
    rank2 <- Immune_fraction[Immune_fraction$`T cells gamma delta`<=0.0,]
    final_rank <- c(rank1[order(rank1$`T cells gamma delta`+rank1$`B cells`+rank1$`Plasma cells`,decreasing = TRUE),]$ID,
                    rank2[order(rank2$`B cells`+rank2$`Plasma cells`,decreasing = TRUE),]$ID)
    final_rank <- Immune_fraction[order(Immune_fraction$`T cells gamma delta`,-Immune_fraction$`NK cells resting`,decreasing = TRUE),]$ID
    
    Immune_fraction_forplot <- reshape2::melt(Immune_fraction,id.vars = "ID")
    #Immune_fraction_forplot <- Immune_fraction_forplot[-grep("CRC",Immune_fraction_forplot$ID),]
    Immune_fraction_forplot <- Immune_fraction_forplot[-grep("NC|PBMC|-T|-N|HCC|LUAD|ESCA|Tumor|Normal|mix|CRC-PKU-5|NC-PKU-mix17|STAD-PKU-4",Immune_fraction_forplot$ID),]
    #Immune_fraction_forplot <- Immune_fraction_forplot[-grep("CRC-PKU-15|STAD-PKU-23|STAD-PKU-36|STAD-PKU-24|STAD-PKU-35",Immune_fraction_forplot$ID),]
    
    Immune_fraction_forplot$ID <- factor(Immune_fraction_forplot$ID,levels = final_rank)
    #Immune_fraction_forplot$ID <- factor(Immune_fraction_forplot$ID,levels = Immune_fraction[order((Immune_fraction$`T cells gamma delta`-Immune_fraction$`CD8 T cells`),decreasing = FALSE),]$ID)
    
    
    #read in clinical data
    {
      clinical <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Sample_clinical_information_x.csv",header = TRUE,check.names = FALSE,row.names = 1)
      clinical$ID <- paste0(rownames(clinical),"-pico")
      clinical <- clinical[-grep("NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",clinical$ID),]
      
      clinical2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Sample_clinical_information_GSE174302.csv",header = TRUE,check.names = FALSE,row.names = 1)
      rownames(clinical2) <- gsub("CRC-2410911","CRC-2410966",rownames(clinical2))
      clinical2$ID <- rownames(clinical2)
      clinical <- rbind(clinical,clinical2)
      
      
      #clinical <- clinical[-grep("NC|PBMC|Tumor|Normal|mix",clinical$ID),]
      rownames(clinical) <- clinical$ID
      clinical <- clinical[as.character(unique(Immune_fraction_forplot$ID)),]
      clinical$ID <- factor(clinical$ID,levels=levels(Immune_fraction_forplot$ID))
      #clinical$ID <- factor(clinical$ID,levels=c(clinical[order(clinical$`Tumor size`,decreasing = TRUE),]$ID))
      #clinical$ID <- factor(clinical$ID,levels=c(clinical[order(clinical$T,decreasing = TRUE),]$ID))
      #Immune_fraction_forplot$ID <- factor(Immune_fraction_forplot$ID,levels = levels(clinical$ID))
      colnames(clinical) <- paste0(colnames(clinical),"_clinical")
    }
    
    #read in TIDE scores
    { 
      TIDE <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Plasma_TIDE.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
      TIDE$ID <- rownames(TIDE)
      TIDE <- TIDE[-grep("NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",TIDE$ID),]
      
      TIDE2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/GSE174302_TIDE.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
      TIDE2$ID <- rownames(TIDE2)
      
      TIDE <- rbind(TIDE,TIDE2)
      
      TIDE$group <- as.character(lapply(strsplit(TIDE$ID,"-"), function(x) x[1]))
      TIDE$group <- gsub("NC","HD",TIDE$group)
      
      
      #my_comparisons <- list(c("STAD","HD"),c("CRC","HD"))
      #Tumor_comparison <- 
      #ggplot(TIDE,aes(x=group,y=Dysfunction,fill = group))+
      #  geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      #  geom_point(size = 1, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
      #  scale_fill_manual(values=c("STAD"="red","CRC"="#FCB514","HD"="blue")) +
      #  theme_bw()+
      #  ylab("Dysfunction score")+
      #  xlab("")+
      #  theme(#legend.position="right",
      #    panel.grid=element_blank(),
      #    panel.border=element_blank(),
      #    axis.line = element_line(size=1, colour = "black"),
      #    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      #    legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
      #    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      #    axis.text.x = element_blank(),
      #    #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
      #    axis.text.y = element_text(face="bold",  color="black", size=24),
      #    axis.title.x = element_text(face="bold", color="black", size=24),
      #    axis.title.y = element_text(face="bold",color="black", size=24))+
      #  stat_compare_means(comparisons = my_comparisons,
      #                     method = "wilcox.test",
      #                     size = 13,
      #                     vjust = 0.6,
      #                     method.args = list(alternative = "two.sided",paired = TRUE),
      #                     label = "p.signif")
      #ggsave(plot = Tumor_comparison, filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Tumor_dysfunction_comparison.pdf",device = "pdf",width = 3.8,height = 5)
      
      TIDE <- TIDE[as.character(unique(Immune_fraction_forplot$ID)),]
      TIDE$group <- factor(TIDE$group,levels = c("CRC","STAD","HD"))
      TIDE$ID <- factor(TIDE$ID,levels=levels(Immune_fraction_forplot$ID))
      
      colnames(TIDE) <- paste0(colnames(TIDE),"_TIDE")
    }
    
    #read in EPIC scores
    {
      EPIC <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Plasma_EPIC_cellFraction_TRef.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
      EPIC$ID <- rownames(EPIC)
      EPIC <- EPIC[-grep("NC-PKU-mix.-pico|CRC-PKU-mix1-pico|CRC-PKU-5-pico|NC-PKU-mix17-pico|STAD-PKU-4-pico",EPIC$ID),]
      
      EPIC2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/GSE174302_EPIC_cellFraction_TRef.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
      EPIC2$ID <- rownames(EPIC2)
      EPIC <- rbind(EPIC,EPIC2)
      
      EPIC <- EPIC[as.character(unique(Immune_fraction_forplot$ID)),]
      EPIC$ID <- factor(EPIC$ID,levels = levels(Immune_fraction_forplot$ID))
      colnames(EPIC) <- paste0(colnames(EPIC),"_EPIC")
    }
    
    #read in gene expression signatures
    {
      T_cell_receptor_signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/T cell receptor signaling pathway_matrix/Multiomics_20211113/Multiomics_20211113_T cell receptor signaling pathway_log2sum.txt",sep = "\t")
      T_cell_receptor_signature2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/T cell receptor signaling pathway_matrix/GSE174302/GSE174302_T cell receptor signaling pathway_log2sum.txt",sep = "\t")
      T_cell_receptor_signature <- rbind(T_cell_receptor_signature,T_cell_receptor_signature2)
      
      B_cell_receptor_signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/B cell receptor signaling pathway_matrix/Multiomics_20211113/Multiomics_20211113_B cell receptor signaling pathway_log2sum.txt",sep = "\t")
      B_cell_receptor_signature2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/B cell receptor signaling pathway_matrix/GSE174302/GSE174302_B cell receptor signaling pathway_log2sum.txt",sep = "\t")
      B_cell_receptor_signature <- rbind(B_cell_receptor_signature,B_cell_receptor_signature2)
      
      Naive_T_cell_signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP_matrix/Multiomics_20211113/Multiomics_20211113_GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP_log2sum.txt",sep = "\t")
      Naive_T_cell_signature2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP_matrix/GSE174302/GSE174302_GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_UP_log2sum.txt",sep = "\t")
      Naive_T_cell_signature <- rbind(Naive_T_cell_signature, Naive_T_cell_signature2)
      
      FOXO_regulated_signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP_matrix/Multiomics_20211113/Multiomics_20211113_GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP_log2sum.txt",sep = "\t")
      FOXO_regulated_signature2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP_matrix/GSE174302/GSE174302_GSE21678_WT_VS_FOXO1_FOXO3_KO_TREG_UP_log2sum.txt",sep = "\t")
      FOXO_regulated_signature <- rbind(FOXO_regulated_signature,FOXO_regulated_signature2)
      
      Effector_T_cell_receptor_signature <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP_matrix/Multiomics_20211113/Multiomics_20211113_GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP_log2sum.txt",sep = "\t")
      Effector_T_cell_receptor_signature2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/01.RP_mRNA_wilcox/GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP_matrix/GSE174302/GSE174302_GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP_log2sum.txt",sep = "\t")
      Effector_T_cell_receptor_signature <- rbind(Effector_T_cell_receptor_signature,Effector_T_cell_receptor_signature2)
      
      Gene_signature <- data.frame(row.names = as.character(unique(Immune_fraction_forplot$ID)),
                                   "T_cell_receptor_signature"=T_cell_receptor_signature[as.character(unique(Immune_fraction_forplot$ID)),],
                                   "B_cell_receptor_signature"=B_cell_receptor_signature[as.character(unique(Immune_fraction_forplot$ID)),],
                                   "Naive_T_cell_signature"=Naive_T_cell_signature[unique(as.character(Immune_fraction_forplot$ID)),],
                                   "FOXO_regulated_signature"=FOXO_regulated_signature[unique(as.character(Immune_fraction_forplot$ID)),],
                                   "Effector_T_cell_receptor_signature"=Effector_T_cell_receptor_signature[unique(as.character(Immune_fraction_forplot$ID)),]
      )
      
      colnames(Gene_signature) <- paste0(colnames(Gene_signature),"_Gene_signature")
    }
    
    #read in CRC scRNAseq deconvolution
    {
      #read in immune fractions
      {
        CRCscRNAseq_fraction <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/20220907_Multiomics_abs_CRC_scRNAseq_cell_type.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
        CRCscRNAseq_fraction2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/output/03.cibersort/20220912_GSE174302_abs_CRC_scRNAseq_cell_type.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
        CRCscRNAseq_fraction <- rbind(CRCscRNAseq_fraction,CRCscRNAseq_fraction2)
        
        CRCscRNAseq_fraction <- CRCscRNAseq_fraction[,-grep("P-value|Correlation|RMSE|Absolute score",colnames(CRCscRNAseq_fraction))]
        
        #CRCscRNAseq_fraction <- CRCscRNAseq_fraction[-grep("CRC.PKU.15|STAD.PKU.23|STAD.PKU.36",rownames(CRCscRNAseq_fraction)),]
        #CRCscRNAseq_fraction <- CRCscRNAseq_fraction[grep("STAD",rownames(CRCscRNAseq_fraction)),]
        #CRCscRNAseq_fraction <- CRCscRNAseq_fraction[-which(rowSums(CRCscRNAseq_fraction)==0),] # This step will remove samples that have no immune cells fraction
        #CRCscRNAseq_fraction <- CRCscRNAseq_fraction/rowSums(CRCscRNAseq_fraction)
        CRCscRNAseq_fraction$ID <- rownames(CRCscRNAseq_fraction)
        CRCscRNAseq_fraction$ID <- gsub(".","-",fixed = TRUE,CRCscRNAseq_fraction$ID)
        rownames(CRCscRNAseq_fraction) <- gsub(".","-",fixed = TRUE,rownames(CRCscRNAseq_fraction))
      }
      
      CRCscRNAseq_fraction <- CRCscRNAseq_fraction[as.character(unique(Immune_fraction_forplot$ID)),]
      #CRCscRNAseq_fraction$group <- factor(CRCscRNAseq_fraction$group,levels = c("CRC","STAD","HD"))
      CRCscRNAseq_fraction$ID <- factor(CRCscRNAseq_fraction$ID,levels=levels(Immune_fraction_forplot$ID))
      colnames(CRCscRNAseq_fraction) <- paste0(colnames(CRCscRNAseq_fraction),"_scRNAseq")
    }
    
    #read in 2021-oncogene bulk cell type
    {
      #read in immune fractions
      {
        bulk_cell_type_fraction <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/20220913_Multiomics_Oncogene_bulk_cell_type_abs.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
        bulk_cell_type_fraction2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/20220913_GSE174302_Oncogene_bulk_cell_type_abs.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
        bulk_cell_type_fraction <- rbind(bulk_cell_type_fraction,bulk_cell_type_fraction2)
        
        bulk_cell_type_fraction <- bulk_cell_type_fraction[,-grep("P-value|Correlation|RMSE|Absolute score",colnames(bulk_cell_type_fraction))]
        
        #bulk_cell_type_fraction <- bulk_cell_type_fraction[-grep("CRC.PKU.15|STAD.PKU.23|STAD.PKU.36",rownames(bulk_cell_type_fraction)),]
        #bulk_cell_type_fraction <- bulk_cell_type_fraction[grep("STAD",rownames(bulk_cell_type_fraction)),]
        #bulk_cell_type_fraction <- bulk_cell_type_fraction[-which(rowSums(bulk_cell_type_fraction)==0),] # This step will remove samples that have no immune cells fraction
        #bulk_cell_type_fraction <- bulk_cell_type_fraction/rowSums(bulk_cell_type_fraction)
        bulk_cell_type_fraction$ID <- rownames(bulk_cell_type_fraction)
        bulk_cell_type_fraction$ID <- gsub(".","-",fixed = TRUE,bulk_cell_type_fraction$ID)
        rownames(bulk_cell_type_fraction) <- gsub(".","-",fixed = TRUE,rownames(bulk_cell_type_fraction))
      }
      
      bulk_cell_type_fraction <- bulk_cell_type_fraction[as.character(unique(Immune_fraction_forplot$ID)),]
      #bulk_cell_type_fraction$group <- factor(bulk_cell_type_fraction$group,levels = c("CRC","STAD","HD"))
      bulk_cell_type_fraction$ID <- factor(bulk_cell_type_fraction$ID,levels=levels(Immune_fraction_forplot$ID))
      colnames(bulk_cell_type_fraction) <- paste0(colnames(bulk_cell_type_fraction),"_oncogene_bulk")
    }
    
    #read in 2021-oncogene single cell type
    {
      #read in immune fractions
      {
        single_cell_type_fraction <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/20220913_Multiomics_Oncogene_single_cell_type_abs.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
        single_cell_type_fraction2 <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 5/20220913_GSE174302_Oncogene_single_cell_type_abs.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
        single_cell_type_fraction <- rbind(single_cell_type_fraction,single_cell_type_fraction2)
        
        single_cell_type_fraction <- single_cell_type_fraction[,-grep("P-value|Correlation|RMSE|Absolute score",colnames(single_cell_type_fraction))]
        
        #single_cell_type_fraction <- single_cell_type_fraction[-grep("CRC.PKU.15|STAD.PKU.23|STAD.PKU.36",rownames(single_cell_type_fraction)),]
        #single_cell_type_fraction <- single_cell_type_fraction[grep("STAD",rownames(single_cell_type_fraction)),]
        #single_cell_type_fraction <- single_cell_type_fraction[-which(rowSums(single_cell_type_fraction)==0),] # This step will remove samples that have no immune cells fraction
        #single_cell_type_fraction <- single_cell_type_fraction/rowSums(single_cell_type_fraction)
        single_cell_type_fraction$ID <- rownames(single_cell_type_fraction)
        single_cell_type_fraction$ID <- gsub(".","-",fixed = TRUE,single_cell_type_fraction$ID)
        rownames(single_cell_type_fraction) <- gsub(".","-",fixed = TRUE,rownames(single_cell_type_fraction))
      }
      
      single_cell_type_fraction <- single_cell_type_fraction[as.character(unique(Immune_fraction_forplot$ID)),]
      #single_cell_type_fraction$group <- factor(single_cell_type_fraction$group,levels = c("CRC","STAD","HD"))
      single_cell_type_fraction$ID <- factor(single_cell_type_fraction$ID,levels=levels(Immune_fraction_forplot$ID))
      colnames(single_cell_type_fraction) <- paste0(colnames(single_cell_type_fraction),"_oncogene_single")
    }
    
    #generate mean signature
    #Paired correlation between tissue and plasma (pathway level, mean of pathway genes expression)
    {
      #KEGG pathway genes
      PATH_ID_NAME <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/09.Published_data/pathway/PATH_ID_NAME_modified.csv",header = TRUE,row.names = 1)
      #multiomics plasma
      counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Tissue_TPM.txt",sep="\t",header = T, row.names = 1)
      #multiomics tissue
      counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/18.differential analysis for each alterations/Expression/matrix/Tissue_TPM.txt",sep="\t",header = T, row.names = 1)
      #multiomics GSE174302
      counts <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/05.RP_mRNA/reference/Expression_matrix_TPM/GSE174302_intron-spanning-available-samples_TPM.txt",sep="\t",header = T, row.names = 1)
      
      counts <- log2(as.matrix(counts)+1)
      
      # set progress
      pb <- progress_bar$new(
        format = "  Processing [:bar] :percent eta: :eta",
        total = length(unique(PATH_ID_NAME$DESCRPTION)), clear = FALSE, width= 60) 
      #get all genes in each pathway in KEGG
      n=1
      pathway_count={}
      while(n <= length(unique(PATH_ID_NAME$DESCRPTION))){
        pathway_name <- as.character(unique(PATH_ID_NAME$DESCRPTION)[n])
        pathway <- subset(PATH_ID_NAME,DESCRPTION==pathway_name)
        
        candidate <- pathway$ensembl_gene_id #for colnames with ensemble id
        #candidate <- pathway$hgnc_symbol #for colnames with gene symbol
        
        candidate <- unique(candidate)
        
        #summary one pathway total count in each sample and make 1 row dataframe
        j=1
        pathway_gene_count={}
        while(j<=length(candidate)){
          target <- candidate[j]
          if(length(grep(target,rownames(counts)))==0) {
            #print(paste0("No ",target," in this dataset."))
            j=j+1
          } else {
            #temp <- counts[which(rownames(counts)==target),]
            temp <- counts[grep(target,rownames(counts),fixed=TRUE),]
            pathway_gene_count <- rbind(pathway_gene_count,temp)
            j=j+1
          }
        }
        if(is.null(pathway_gene_count)){
          n=n+1
          next;
        } else {
          one_pathway <- as.data.frame(t(colMeans(pathway_gene_count))) #colSums
          rownames(one_pathway) <- pathway_name
          pathway_count <- rbind(pathway_count,one_pathway)
          n=n+1
        }
        pb$tick()
        Sys.sleep(1 / 100)
        
      }
      write.table(pathway_count,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/GSE174302_pathway_mean.txt",sep = "\t", quote = FALSE)
    }
    #mean sigantures
    {
      Gene_signature <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/Plasma_pathway_mean.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
      Gene_signature$ID <- rownames(Gene_signature)
      Gene_signature2 <- read.table("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/22.tissue_paired_analysis/GSE174302_pathway_mean.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t') #读入数据
      Gene_signature2$ID <- rownames(Gene_signature2)
      Gene_signature <- left_join(Gene_signature,Gene_signature2,by=c("ID"="ID"))
      rownames(Gene_signature) <- Gene_signature$ID
      Gene_signature <- Gene_signature[,-which(colnames(Gene_signature)=="ID")]
      Gene_signature <- as.data.frame(t(Gene_signature))
      rownames(Gene_signature) <- gsub(".","-",fixed = TRUE,rownames(Gene_signature))
      Gene_signature$ID <- rownames(Gene_signature)
      Gene_signature <- Gene_signature[as.character(unique(Immune_fraction_forplot$ID)),]
      Gene_signature$ID <- factor(Gene_signature$ID,levels = levels(Immune_fraction_forplot$ID))
      Gene_signature <- Gene_signature[,-which(colnames(Gene_signature)=="ID")]
    }
    
    #Correlation and comparison
    {
      forCorrlation <- cbind(clinical,Immune_fraction[unique(as.character(Immune_fraction_forplot$ID)),],TIDE,EPIC,Gene_signature,CRCscRNAseq_fraction,single_cell_type_fraction,bulk_cell_type_fraction)
      
      forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="Stage_clinical")]
      forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="ID_clinical")]
      forCorrlation <- forCorrlation[,-grep("ID_EPIC|ID_TIDE|ID_oncogene_bulk|ID_oncogene_single|ID_scRNAseq|No benefits|MSI Score|Patient ID|MSI/MMR|group|Diagnosis|CTL.flag|Responder",colnames(forCorrlation))]
      #CRC
      #forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="M")]
      #forCorrlation <- forCorrlation[,-grep("No benefits|MSI Score|Type|Patient ID|MSI/MMR|group|Diagnosis|CTL.flag|Responder|Multiple primary|MUC2|MUC5AC|EBER|EGFR|CDNA2",colnames(forCorrlation))]
      #STAD
      #forCorrlation <- forCorrlation[,-grep("No benefits|MSI Score|Type|Patient ID|MSI/MMR|group|Diagnosis|CTL.flag|Responder|Right half|EGFR|MLH1|MSH2|MSH6|PMS2",colnames(forCorrlation))]
      #forCorrlation <- forCorrlation[-grep("CRC",rownames(forCorrlation)),]
      
      {
        forCorrlation$Stage_simplified_clinical <- factor(forCorrlation$Stage_simplified_clinical,levels = c("early","late","x"))
        #forCorrlation$`Multiple primary cancer(Yes/No)` <- factor(forCorrlation$`Multiple primary cancer(Yes/No)`,levels = c())
        #forCorrlation$`Right half / left half` <- factor(forCorrlation$`Right half / left half`,levels = c())
        #forCorrlation$`Tumor size` <- factor(forCorrlation$,levels = c())
        forCorrlation$`T_clinical` <- factor(forCorrlation$`T_clinical`,levels = c("Tis","1","1a","1b","2","3","4","4a","4b","x"))
        #forCorrlation$`N` <- factor(forCorrlation$,levels = c())
        #forCorrlation$`M` <- factor(forCorrlation$,levels = c())
        forCorrlation$`Vascular tumor thrombus_clinical` <- factor(forCorrlation$`Vascular tumor thrombus_clinical`,levels = c("No","Yes","x"))
        forCorrlation$`Neurological invasion_clinical` <- factor(forCorrlation$`Neurological invasion_clinical`,levels = c("No","Yes","x"))
        forCorrlation$`Tumor deposition/cancer nodules_clinical` <- factor(forCorrlation$`Tumor deposition/cancer nodules_clinical`,levels = c("No","Yes","x"))
        #forCorrlation$`CEA(HE)` <- factor(forCorrlation$,levels = c())
        #forCorrlation$`CA199(HE)` <- factor(forCorrlation$,levels = c())
        #forCorrlation$`EGFR` <- factor(forCorrlation$,levels = c())
        #forCorrlation$`HER2` <- factor(forCorrlation$,levels = c())
        #forCorrlation$`CDx2` <- factor(forCorrlation$,levels = c())
        #forCorrlation$P53 <- factor(forCorrlation$,levels = c())
        #forCorrlation$`TOP IIa` <- factor(forCorrlation$,levels = c())
        #forCorrlation$Ki67 <- factor(forCorrlation$,levels = c())
        #forCorrlation$`MMR(p/d)` <- factor(forCorrlation$,levels = c())
        #forCorrlation$MLH1 <- factor(forCorrlation$,levels = c())
        #forCorrlation$PMS2 <- factor(forCorrlation$,levels = c())
        #forCorrlation$MSH2 <- factor(forCorrlation$,levels = c())
        #forCorrlation$MSH6 <- factor(forCorrlation$,levels = c())
        #forCorrlation$MUC2 <- factor(forCorrlation$,levels = c())
        #forCorrlation$MUC5AC <- factor(forCorrlation$,levels = c())
        #forCorrlation$EBER <- factor(forCorrlation$,levels = c())
      }
      
      correlation_method <- "spearman"
      output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 6/Expression_correlation/Multiomics+GSE174302_final/"
      enmuerate_correlation(forCorrlation,correlation_method,output_dir,ncol(clinical))
      
      #Visualization
      #method1
      {
        correlation_method <- "spearman"
        output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation/GSE174302_STAD_"
        #enmuerate_correlation(forCorrlation,correlation_method,output_dir)
        forCorrlation <- forCorrlation[-grep("pico",rownames(forCorrlation)),]
        enmuerate_correlation(forCorrlation[grep("STAD",rownames(forCorrlation)),],correlation_method,output_dir)
        
        result <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation/CRC_STAD_Result_final.csv",header = TRUE, row.names = 1,check.names = FALSE)
        result$column <- as.character(lapply(strsplit(as.character(result$DataNames)," vs ",fixed = TRUE), function(x) x[1]))
        result$row <- as.character(lapply(strsplit(as.character(result$DataNames)," vs ",fixed = TRUE), function(x) x[2]))
        R <- data.frame(matrix(nrow=length(colnames(forCorrlation)),ncol = length(colnames(forCorrlation))))
        rownames(R) = colnames(forCorrlation)
        colnames(R) = colnames(forCorrlation)
        pvalue <- data.frame(matrix(nrow=length(colnames(forCorrlation)),ncol = length(colnames(forCorrlation))))
        rownames(pvalue) = colnames(forCorrlation)
        colnames(pvalue) = colnames(forCorrlation)
        
        i=1
        while(i<=nrow(result)){
          pvalue[result[i,"column"],result[i,"row"]] <- as.numeric(as.character(result[i,"Pvalue"]))
          pvalue[result[i,"row"],result[i,"column"]] <- as.numeric(as.character(result[i,"Pvalue"]))
          R[result[i,"column"],result[i,"row"]] <- as.numeric(as.character(result[i,"R"]))
          R[result[i,"row"],result[i,"column"]] <- as.numeric(as.character(result[i,"R"]))
          i=i+1
        }
        
        i=1
        while(i<=nrow(R)){
          pvalue[i,i] <- 0
          R[i,i] <- 1
          i=i+1
        }
        
        R[is.na(R)] <- 0
        pvalue[is.na(pvalue)] <- 1
        
        r = rbind(c('Plasma cells','Type','Mast cells','EBER'),
                  c('TIDE','Type','TAM M2','EBER'),
                  c('Bcells_EPIC','Type','otherCells_EPIC','EBER'),
                  c('T_cell_receptor_signature','Type','Effector_T_cell_receptor_signature','EBER'),
                  c('Plasma cells','Type','Effector_T_cell_receptor_signature','Age'),
                  c('Plasma cells','Stage_raw','Effector_T_cell_receptor_signature','Multiple primary cancer'),
                  c('Plasma cells','Hgb','Effector_T_cell_receptor_signature','PT'),
                  c('Plasma cells','CEA','Effector_T_cell_receptor_signature','CA199(HE)'),
                  c('Plasma cells','EGFR','Effector_T_cell_receptor_signature','EBER')
        )
        corrplot(corr = as.matrix(R),tl.col="black",order="original",tl.pos = "ld",tl.cex=0.7,tl.srt = 45, type="lower", col = col2(200), #mar = c(1, 1, 1, 1),
                 p.mat = as.matrix(pvalue), sig.level = 0.05,insig = "blank",pch.cex = 2) %>% corrRect(namesMat = r)
      }
      
      #method2
      {
        library(corrplot)
        library(RColorBrewer)
        library("Hmisc")
        library(magrittr)
        col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                       "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                       "#4393C3", "#2166AC", "#053061")))
        res2 <- rcorr(as.matrix(sapply(forCorrlation, as.numeric))  ,type="spearman")
        
        res2$r[is.na(res2$r)] <- 1
        res2$P[is.na(res2$P)] <- 0
        
        r = rbind(c('Plasma cells','Type','Mast cells','EBER'),
                  c('TIDE','Type','TAM M2','EBER'),
                  c('Bcells','Type','otherCells','EBER'),
                  c('T_cell_receptor_signature','Type','Effector_T_cell_receptor_signature','EBER'),
                  c('Plasma cells','Type','Effector_T_cell_receptor_signature','Age'),
                  c('Plasma cells','Stage_raw','Effector_T_cell_receptor_signature','Multiple primary cancer'),
                  c('Plasma cells','Hgb','Effector_T_cell_receptor_signature','PT'),
                  c('Plasma cells','CEA','Effector_T_cell_receptor_signature','CA199(HE)'),
                  c('Plasma cells','EGFR','Effector_T_cell_receptor_signature','EBER')
        )
        
        r = rbind(c('B cells','Gender','Neutrophils','EBER'),
                  c('TIDE','Gender','TAM M2','EBER'),
                  c('Bcells','Gender','otherCells','EBER'),
                  c('T_cell_receptor_signature','Gender','Effector_T_cell_receptor_signature','EBER'),
                  c('B cells','Gender','Effector_T_cell_receptor_signature','Age'),
                  c('B cells','Stage_raw','Effector_T_cell_receptor_signature','Multiple primary cancer'),
                  c('B cells','Hgb','Effector_T_cell_receptor_signature','PT'),
                  c('B cells','CEA','Effector_T_cell_receptor_signature','CA199(HE)'),
                  c('B cells','HER2','Effector_T_cell_receptor_signature','EBER')
        )
        
        r = rbind(c('B cells','Gender','Neutrophils','MSH6'),
                  c('TIDE','Gender','TAM M2','MSH6'),
                  c('Bcells','Gender','otherCells','MSH6'),
                  c('T_cell_receptor_signature','Gender','Effector_T_cell_receptor_signature','MSH6'),
                  c('B cells','Gender','Effector_T_cell_receptor_signature','Age'),
                  c('B cells','Stage_raw','Effector_T_cell_receptor_signature','Tumor deposition/cancer nodules'),
                  c('B cells','Hgb','Effector_T_cell_receptor_signature','PT'),
                  c('B cells','CEA','Effector_T_cell_receptor_signature','CA199(HE)'),
                  c('B cells','EGFR','Effector_T_cell_receptor_signature','MSH6')
        )
        
        corrplot(corr = res2$r,tl.col="black",order="original",tl.pos = "ld",tl.cex=0.7,tl.srt = 45, type="lower", col = col2(200), #mar = c(1, 1, 1, 1),
                 p.mat = res2$P, sig.level = 0.05,insig = "blank",pch.cex = 2) %>% corrRect(namesMat = r)
      }
      
    }
    
  }
  
  #plot cibersort barplot
  {
    barplot_cibersort <- read.csv("Figure 6/Expression_correlation/barplot_LM22.csv")
    barplot_cibersort$Cell.fraction <- factor(barplot_cibersort$Cell.fraction,levels = barplot_cibersort[order(barplot_cibersort$Rank,decreasing = TRUE),]$Cell.fraction)
    #barplot_cibersort <- barplot_cibersort[-grep("Li et al., 2021: ",barplot_cibersort$Cell.fraction),]
    p1 <- ggplot(barplot_cibersort,aes(x = Cell.fraction, y=R))+
      geom_bar(stat = "identity",aes(fill = Trend,color = Sig.))+
      scale_fill_manual(values=c("Negative"="light blue","Positive"="pink"))+
      scale_color_manual(values=c("Non-significant"="white","Significant"="Red"))+
      scale_y_continuous(breaks = c(0,0.05,0.10,0.15,0.20,0.25,0.30),labels = c("0.00","0.05","0.10","0.15","0.20","0.25","0.30"),expand = c(0,0),limits = c(0,0.3))+
      theme_bw()+
      coord_flip()+
      theme(
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=16),
        legend.text = element_text(color="black",family = "Arial", size=16),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=20, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    barplot_cibersort <- read.csv("Figure 6/Expression_correlation/barplot_non-immune.csv")
    barplot_cibersort$Cell.fraction <- factor(barplot_cibersort$Cell.fraction,levels = barplot_cibersort[order(barplot_cibersort$Rank,decreasing = TRUE),]$Cell.fraction)
    #barplot_cibersort <- barplot_cibersort[-grep("Li et al., 2021: ",barplot_cibersort$Cell.fraction),]
    p2 <- ggplot(barplot_cibersort,aes(x = Cell.fraction, y=R))+
      geom_bar(stat = "identity",aes(fill = Trend,color = Sig.))+
      scale_fill_manual(values=c("Negative"="light blue","Positive"="pink"))+
      scale_color_manual(values=c("Non-significant"="white","Significant"="Red"))+
      scale_y_continuous(breaks = c(0,0.05,0.10,0.15,0.20,0.25,0.30),labels = c("0.00","0.05","0.10","0.15","0.20","0.25","0.30"),expand = c(0,0),limits = c(0,0.3))+
      theme_bw()+
      coord_flip()+
      theme(
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=16),
        legend.text = element_text(color="black",family = "Arial", size=16),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(color="black", size=16, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(color="black", size=20, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_text(color="black", size=20, angle = 0,hjust=0.5,vjust = 0.5),
        axis.title.y = element_blank())
    p <- ggarrange(p1,p2,ncol=1,heights = c(2,1),align = "v",legend = FALSE)
    ggsave(p,filename = "./Figure 6/Barplot_non-immune.pdf",device = "pdf",width = 8.56,height = 11)
  }
  
  #complex heatmap
  {
    Immune_fraction_forplot$variable <- factor(Immune_fraction_forplot$variable, levels = rev(c("T cells gamma delta","NK cells resting","NK cells activated",
                                                                                                "T cells CD4","T cells CD8","B cells","Plasma cells",
                                                                                                "Monocytes","Macrophages","Dendritic cells",
                                                                                                "Mast cells","Eosinophils","Neutrophils")))
    Immune_fraction_forplot_selected <- Immune_fraction_forplot[grep("CRC|STAD",Immune_fraction_forplot$ID),]
    clinical_selected <- clinical[grep("CRC|STAD",clinical$ID),]
    TIDE_selected <- TIDE[grep("CRC|STAD",TIDE$ID),]
    EPIC_selected <- EPIC[grep("CRC|STAD",EPIC$ID_EPIC),]
    
    simpsons <- c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF",
                  "#197EC0FF","#F05C3BFF","#46732EFF","#71D0F5FF","#370335FF","#075149FF",
                  "#C80813FF","#91331FFF","#1A9993FF","#FD8CC1FF")
    
    simpsons <- c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF",
                  "#197EC0FF",alpha("purple",alpha=0.4),alpha("brown",alpha=0.8),"#71D0F5FF","#EE2C2C","#370335FF","dark green",
                  "#C80813FF","#91331FFF","#1A9993FF","#FD8CC1FF")
    bar <- 
      ggplot(Immune_fraction_forplot_selected,aes(x=ID,y=value,fill=variable))+
      geom_bar(stat = "identity",position = "stack")+
      xlab("")+
      ylab("")+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,10,-20,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        #axis.text.x = element_text( color="black", size=16, angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(color="black", size=16, angle = 90,hjust=0.5, family = "Arial"),
        axis.title.x = element_text(face="bold", color="black", size=20, family = "Arial"),
        axis.title.y = element_text(face="bold",color="black", size=20, family = "Arial"))+
      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),labels = c("0","0.25","0.50","0.75","1"),expand = c(0,0),limits = c(0,1))+
      scale_fill_manual(values = simpsons)
    #scale_fill_rickandmorty(alpha = 1)
    
    
    #clinical_selected$`Tumor size` <- gsub("x","NA",clinical_selected$`Tumor size`)
    clinical_TumorSize <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= as.numeric(as.character(clinical_selected$`Tumor size_clinical`))))+
      geom_tile()+
      scale_fill_gradient2(low = "white",
                           mid = "white",
                           high = "red",
                           midpoint = 5,na.value = "grey")+
      #scale_fill_manual(values=c("-"="white","early"="sky blue","late"="red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    clinical_PLT <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= as.numeric(as.character(clinical_selected$PLT_clinical))))+#`Tumor size`))+
      geom_tile()+
      scale_fill_gradient(low="white",high = "orange")+
      #scale_fill_manual(values=c("-"="white","early"="sky blue","late"="red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    clinical_P53 <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= P53_clinical))+#`Tumor size`))+
      geom_tile()+
      #scale_fill_gradient(low="white",high = "black")+
      scale_fill_manual(values=c("NA"="grey","-"="white","+"="white","++"="dark grey","+++"="black"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    clinical_MMR <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= `MSI/MMR_clinical`))+#`Tumor size`))+
      geom_tile()+
      #scale_fill_gradient(low="white",high = "black")+
      scale_fill_manual(values=c("NA"="grey","x"="grey","Normal"="white","Unstable; lack"="#F05C3BFF","unstable; lack"="#F05C3BFF"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    clinical_HER2 <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= HER2_clinical))+#`Tumor size`))+
      geom_tile()+
      #scale_fill_gradient(low="white",high = "black")+
      scale_fill_manual(values=c("NA"="grey","-"="white","+"="white","++"="red","+++"="red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    clinical_Right <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= `Right half / left half_clinical`))+#`Tumor size`))+
      geom_tile()+
      #scale_fill_gradient(low="white",high = "black")+
      scale_fill_manual(values=c("NA"="grey","x"="grey","L"="white","R"="dark red","L+R"="dark red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    clinical_Stage <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= Stage_clinical))+#`Tumor size`))+
      geom_tile()+
      #scale_fill_gradient(low="white",high = "black")+
      scale_fill_manual(values=c("NA"="grey","x"="grey","1"="#90E0EF","1A"="#90E0EF","1B"="#90E0EF","2A"="#00B4D8","2B"="#00B4D8","2C"="#00B4D8","3A"="#03045E","3B"="#03045E","3C"="#03045E","4"="black"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    clinical_T <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= `T_clinical`))+#`Tumor size`))+
      geom_tile()+
      #scale_fill_gradient(low="white",high = "black")+
      scale_fill_manual(values=c("NA"="grey","Tis"="#CAF0F8","1"="#90E0EF","1a"="#90E0EF","1b"="#00B4D8","2"="#00B4D8","3"="#03045E","4a"="black","4b"="black"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    clinical_N <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= N_clinical))+
      geom_tile()+
      #scale_fill_gradient(low="white",high = "black")+
      scale_fill_manual(values=c("x"="grey","0"="#CAF0F8","1"="#90E0EF","1a"="#90E0EF","1b"="#90E0EF","2"="#00B4D8","2a"="#00B4D8","2b"="#00B4D8","3a"="#03045E","3b"="#03045E"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    clinical_M <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= as.factor(clinical_selected$M_clinical)))+
      geom_tile()+
      #scale_fill_gradient(low="white",high = "black")+
      scale_fill_manual(values=c("0"="white","1"="black"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    clinical_Cancer <-
      ggplot(clinical_selected,aes(x=ID_clinical,y=0,fill= Type_clinical))+
      geom_tile()+
      #scale_fill_gradient(low="white",high = "black")+
      scale_fill_manual(values=c("Colorectal Cancer"="#FCB514","Stomach Cancer"="red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    TIDE_Responder <-
      ggplot(TIDE_selected,aes(x=ID_TIDE,y=0,fill= Responder_TIDE))+
      geom_tile()+
      #scale_fill_gradient(low="white",high = "blue")+
      scale_fill_manual(values=c("True"="dark blue","False"="white"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    TIDE_IFNG <-
      ggplot(TIDE_selected,aes(x=ID_TIDE,y=0,fill= IFNG_TIDE))+
      geom_tile()+
      scale_fill_gradient(low="white",high = "sky blue")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    TIDE_Dysfunction <-
      ggplot(TIDE_selected,aes(x=ID_TIDE,y=0,fill= Dysfunction_TIDE))+
      geom_tile()+
      scale_fill_gradient2(low="white",mid="white",midpoint = -0.25,high = "dark green")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      #title("MDSC+CAF")+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    TIDE_Exclusion <-
      ggplot(TIDE_selected,aes(x=ID_TIDE,y=0,fill= Exclusion_TIDE))+
      geom_tile()+
      scale_fill_gradient2(low="white",mid="white",midpoint = 0,high = "orange")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      #title("MDSC+CAF")+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    TIDE_CAF <-
      ggplot(TIDE_selected,aes(x=ID_TIDE,y=0,fill= CAF_TIDE))+
      geom_tile()+
      scale_fill_gradient2(low="white",mid="white",midpoint = 0,high = "dark orange")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      #title("MDSC+CAF")+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    TIDE_TAM <-
      ggplot(TIDE_selected,aes(x=ID_TIDE,y=0,fill= `TAM M2_TIDE`))+
      geom_tile()+
      scale_fill_gradient2(low="white",mid="white",midpoint = 0,high = "dark green")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      #title("MDSC+CAF")+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    TIDE_TIDE <-
      ggplot(TIDE_selected,aes(x=ID_TIDE,y=0,fill= TIDE_TIDE))+
      geom_tile()+
      scale_fill_gradient2(low="white",mid="white",midpoint = 0,high = "orange")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      #title("MDSC+CAF")+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    EPIC_CAFsEndothelial <-
      ggplot(EPIC_selected,aes(x=ID_EPIC,y=0,fill= EPIC_selected$CAFs_EPIC+EPIC_selected$Endothelial_EPIC))+
      geom_tile()+
      scale_fill_gradient(low="white",high = "blue")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    EPIC_CAFs <-
      ggplot(EPIC_selected,aes(x=ID_EPIC,y=0,fill= CAFs_EPIC))+
      geom_tile()+
      scale_fill_gradient2(low="white",mid="white",midpoint = 0,high = "dark orange")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    EPIC_Macrophages <-
      ggplot(EPIC_selected,aes(x=ID_EPIC,y=0,fill= Macrophages_EPIC))+
      geom_tile()+
      scale_fill_gradient(low="white",high = "blue")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    EPIC_CAFs <-
      ggplot(EPIC_selected,aes(x=ID_EPIC,y=0,fill= CAFs_EPIC))+
      geom_tile()+
      scale_fill_gradient2(low="white",mid="white",midpoint = 0,high = "dark orange")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    EPIC_Bcells <-
      ggplot(EPIC_selected,aes(x=ID_EPIC,y=0,fill= Bcells_EPIC))+
      geom_tile()+
      scale_fill_gradient(low="white",high = "blue")+
      #scale_fill_manual(values=c("True"="dark blue","False"="red"))+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(0,0,0,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        #legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x= element_blank(),
        #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    ggarrange(bar,clinical_TumorSize,clinical_PLT,clinical_P53,clinical_MMR,clinical_HER2,clinical_Right,clinical_T,clinical_N,clinical_M,TIDE_Responder,TIDE_IFNG,TIDE_Dysfunction,TIDE_Exclusion,EPIC_CAFsEndothelial,
              ncol = 1,nrow = 14,align = "v",heights = c(100,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5),legend = "none")
    
    ggarrange(bar,clinical_TumorSize,#TIDE_bar3,
              ncol = 1,nrow = 2,align = "v",heights = c(100,10),legend = "none")
    
    ggarrange(clinical_Cancer,bar,TIDE_TAM,EPIC_CAFs,clinical_Right,clinical_MMR,clinical_TumorSize,clinical_M,clinical_Stage,
              ncol = 1,nrow = 9,align = "v",heights = c(5,100,10,10,10,10,10,10,10),legend = "none")
    
  }
  
  #boxplot
  #comparison
  {
    forCorrlation$boxplot_group <- NA
    forCorrlation[forCorrlation$`T cells gamma delta`>0,]$boxplot_group <- "gamma_delta_postive"
    forCorrlation[forCorrlation$`T cells gamma delta`<=0,]$boxplot_group <- "gamma_delta_negative"
    
    forCorrlation$Stage_group <- NA
    forCorrlation[forCorrlation$Stage_raw=="1",]$Stage_group <- "Stage I"
    forCorrlation[forCorrlation$Stage_raw=="1A",]$Stage_group <- "Stage I"
    forCorrlation[forCorrlation$Stage_raw=="1B",]$Stage_group <- "Stage I"
    forCorrlation[forCorrlation$Stage_raw=="2A",]$Stage_group <- "Stage II"
    forCorrlation[forCorrlation$Stage_raw=="2B",]$Stage_group <- "Stage II"
    forCorrlation[forCorrlation$Stage_raw=="2C",]$Stage_group <- "Stage II"
    forCorrlation[forCorrlation$Stage_raw=="3A",]$Stage_group <- "Stage III/IV"
    forCorrlation[forCorrlation$Stage_raw=="3B",]$Stage_group <- "Stage III/IV"
    forCorrlation[forCorrlation$Stage_raw=="3C",]$Stage_group <- "Stage III/IV"
    forCorrlation[forCorrlation$Stage_raw=="4",]$Stage_group <- "Stage III/IV"
    
    
    my_comparisons <- list(c("gamma_delta_postive","gamma_delta_negative"))
    my_comparisons <- list(c("early","late"))
    #my_comparisons <- list(c("Stage I","Stage II"),c("Stage I","Stage III"),c("Stage I","Stage IV"))
    my_comparisons <- list(c("Stage II","Stage III/IV"),c("Stage I","Stage III/IV"))
    forCorrlation_boxplot <- forCorrlation[grep("Colorectal Cancer",forCorrlation$Type),]
    forCorrlation_boxplot <- forCorrlation
    Tumor_comparison <- 
      ggplot(forCorrlation_boxplot[forCorrlation_boxplot$Stage_simplified!="x",],aes(x=Stage_group, y = `TAM M2_TIDE`, fill = Stage_group))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
      scale_fill_manual(values=c("gamma_delta_postive"="dark green","gamma_delta_negative"="white",
                                 "early"="sky blue","late"="black",
                                 "Stage I"="#90E0EF","Stage II"="#00B4D8","Stage III"="#03045E","Stage III/IV"="#03045E","Stage IV"="black")) +
      theme_bw()+
      ylab("M2 TAM")+
      xlab("")+
      theme(#legend.position="right",
        #element_text(family = "Arial"),
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(color="black", size=24,family = "Arial"),
        axis.title.x = element_text(color="black", size=24,family = "Arial"),
        axis.title.y = element_text(color="black", size=24,family = "Arial"))+
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         size = 13,
                         vjust = 0.6,
                         method.args = list(alternative = "less",paired = FALSE),
                         label = "p.signif")
    Tumor_comparison
    ggsave(plot = Tumor_comparison, filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Stage_TAM_CRC+STAD_multiomics+GSE174302_comparison.pdf",device = "pdf",width = 4,height = 5)
  }
  
  #get CAF genes
  {
    #method1:from EPIC, elife-26476-supp2-v2, CAF expression / all cell expression > 0.99 genes
    library(biomaRt)
    genes <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/NK_LM22_genes.csv",header=TRUE,check.names = FALSE)
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    gene_names <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol","chromosome_name"),
                        filters = "hgnc_symbol",
                        values=as.character(genes$Gene.Symbol), mart= mart,useCache = FALSE)
    gene_names <- gene_names[-grep("CHR",gene_names$chromosome_name),]
    gene_names <- left_join(gene_names,genes,by=c("hgnc_symbol"="Gene.Symbol"))
    write.csv(gene_names,"/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/genes_ensembl_NK_cells_LM22.csv")
    
    #method2:from TIDE, and the reference 83 Calon, A. et al. Dependency of colorectal cancer on a TGF-beta-driven program in stromal cells for metastasis initiation. Cancer cell 22, 571-584, doi:10.1016/j.ccr.2012.08.013 (2012).
    library(affy)
    library(limma)
    library(biomaRt)
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    
    path_to_cel <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/GSE39396_RAW/"
    datafiles <- ReadAffy(celfile.path = path_to_cel)
    eset <- rma(datafiles)
    View(eset)
    write.exprs(eset,file="/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/GSE39396.txt")
    
    Expression_matrix <- data.frame(exprs(eset))
    ph<-pData(eset) 
    design <- model.matrix(~factor(p_disease)) 
    colnames(design) <- c("case","control") 
    
    fit <- lmFit(eset, design) 
    fit <- eBayes(fit) 
    options(digits=2) 
    
    genes<- topTable(fit, coef=2, n=40, adjust="BH") 
    
    gene_names <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol","chromosome_name"),
                        filters = "affy_ht_hg_u133_plus_pm",
                        values=rownames(Expression_matrix), mart= mart,useCache = FALSE)
    
  }
  
  #gene level correlation
  {
    Plasma_TPM <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/Plasma_PBMC_tissue_TPM_noduplicated_gene.txt",sep = "\t",check.names = FALSE)
    GSE173402_TPM <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/TPM/GSE174302_intron-spanning-available-samples_TPM.txt",sep = "\t",check.names = FALSE)
    colnames(GSE173402_TPM) <- gsub(".","-",fixed=TRUE,colnames(GSE173402_TPM))
    Plasma_TPM$ID <- rownames(Plasma_TPM)
    GSE173402_TPM$ID <- rownames(GSE173402_TPM)
    
    Plasma_TPM <- left_join(Plasma_TPM,GSE173402_TPM, by = c("ID"="ID"))
    Plasma_TPM[is.na(Plasma_TPM)] <- 0
    rownames(Plasma_TPM) <- Plasma_TPM$ID
    Plasma_TPM <- Plasma_TPM[,-which(colnames(Plasma_TPM)=="ID")]
    
    #LM22_genes <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/16.origin/signature/LM22_ENSG.txt",sep = "\t")
    #Marker_genes <- Plasma_TPM[as.character(LM22_genes$Gene.symbol),as.character(unique(Immune_fraction_forplot$ID))]
    
    Interested_genes <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Genes/genes_ensembl_NK_cells_LM22.csv")
    Marker_genes <- Plasma_TPM[unique(as.character(Interested_genes$ensembl_gene_id)),as.character(unique(Immune_fraction_forplot$ID))]
    Marker_genes <- na.omit(Marker_genes)
    
    forCorrlation <- cbind(clinical,t(Marker_genes))
    
    #forCorrlation <- forCorrlation[grep("pico",rownames(forCorrlation)),]
    #forCorrlation <- forCorrlation[which(forCorrlation$Type=="Colorectal Cancer"),]
    
    #forCorrlation <- forCorrlation[which(forCorrlation$Type=="Stomach Cancer"),]
    forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="Stage")]
    forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="ID")]
    forCorrlation <- forCorrlation[,-which(colnames(forCorrlation)=="N")]
    forCorrlation <- forCorrlation[,-grep("Type|Patient ID|Gender|Age|Stage_simplified|T|M|Diagnosis|Tumor size|Vascular tumor thrombus|Neurological invasion",colnames(forCorrlation))]
    forCorrlation <- forCorrlation[,-grep("Tumor deposition/cancer nodules|Multiple primary cancer|Hgb|PLT|ALT|AST|ALB|PT|CA724|CA242|CEA|CA199|EGFR|HER2|CDNA2|TOP IIa|Right half / left half|MMR|MLH1|PMS2|MSH2|MSH6|MUC2|MUC5AC|EBER",colnames(forCorrlation))]
    
    correlation_method <- "spearman"
    output_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Multiomics_NK_cell_LM22_genes/Multiomics_NK_cell_LM22"
    dir.create("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/20.Immune/Correlation_forGenes/Multiomics_NK_cell_LM22_genes/")
    enmuerate_correlation(forCorrlation,correlation_method,output_dir)
  }
  
  #TCGA#Survival analysis & stage analysis
  {
    setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/21.Survival analysis")
    library(survival)
    library(survminer)
    library(TCGAutils)
    library(dplyr)
    
    #read sample sheet
    sample_sheet1 <- read.table('COAD_TCGA/COAD_sample_sheet.2019-06-08.tsv',header=T,sep='\t',check.names = FALSE)
    sample_sheet2 <- read.table('READ_TCGA/READ_sample_sheet.2019-06-08.tsv',header=T,sep='\t',check.names = FALSE)
    sample_sheet3 <- read.table('STAD_TCGA/STAD_sample_sheet.2019-06-09.tsv',header=T,sep='\t',check.names = FALSE)
    sample_sheet <- rbind(sample_sheet1,sample_sheet2,sample_sheet3)
    sample_sheet <- sample_sheet[grep("htseq.counts.gz",sample_sheet$`File Name`),]
    sample_sheet$`File Name` <- as.character(lapply(strsplit(as.character(sample_sheet$`File Name`),".",fixed = TRUE), function(x) x[1]))
    
    head(sample_sheet)
    
    # read RNA file 
    rna1 <- read.table('COAD_TCGA/TCGA_COAD_count_matrix_TPM.txt',header=T,row.names=1,sep='\t',check.names = FALSE)
    rna2 <- read.table('READ_TCGA/TCGA_READ_count_matrix_TPM.txt',header=T,row.names=1,sep='\t',check.names = FALSE)
    rna3 <- read.table('STAD_TCGA/TCGA_STAD_count_matrix_TPM.txt',header=T,row.names=1,sep='\t',check.names = FALSE)
    rna <- cbind(rna1,rna2,rna3)
    rna_filename <- data.frame("file_name"=colnames(rna))
    rna_caseid <- left_join(rna_filename,sample_sheet,by=c("file_name"="File Name"))
    
    # and read the Clinical file, in this case i transposed it to keep the clinical feature title as column name
    clinical1 <- t(read.table('COAD_TCGA/clinical.tsv',header=T, row.names=1, sep='\t'))
    clinical2 <- t(read.table('READ_TCGA/clinical.tsv',header=T, row.names=1, sep='\t'))
    clinical3 <- t(read.table('STAD_TCGA/clinical.tsv',header=T, row.names=1, sep='\t'))
    clinical <- cbind(clinical1,clinical2,clinical3)
    clinical_caseid <- UUIDtoBarcode(colnames(clinical), from_type = "case_id")
    
    complete_id <- left_join(rna_caseid,clinical_caseid,by=c("Case ID"="submitter_id"))
    complete_id <- complete_id[complete_id$`Sample Type`=="Primary Tumor",]
    
    #substitute
    rna_forSurvival <- rna[,complete_id[grep("TCGA-STAD|TCGA-READ|TCGA-COAD",complete_id$`Project ID`),]$file_name]
    colnames(rna_forSurvival) <- complete_id[grep("TCGA-STAD|TCGA-READ|TCGA-COAD",complete_id$`Project ID`),]$`Sample ID`
    clinical_forSurvival <- clinical[,complete_id[grep("TCGA-STAD|TCGA-READ|TCGA-COAD",complete_id$`Project ID`),]$case_id]
    colnames(clinical_forSurvival) <- complete_id[grep("TCGA-STAD|TCGA-READ|TCGA-COAD",complete_id$`Project ID`),]$`Sample ID`
    
    clinical_forSurvival.t <- as.data.frame(t(clinical_forSurvival))
    No_status <- paste(names(grep("--|Not Reported",clinical_forSurvival.t$vital_status,value = TRUE)),collapse = "|")
    complete_id <- complete_id[-grep(No_status,complete_id$`Sample ID`),]
    #rm(rna)
    #rm(clinical)
    
    # first I remove genes whose expression is == 0 in more than 50% of the samples:
    rem <- function(x){
      x <- as.matrix(x)
      x <- t(apply(x,1,as.numeric))
      r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
      remove <- which(r > dim(x)[2]*0.5)
      return(remove)
    }
    remove <- rem(rna_forSurvival)
    rna_forSurvival <- rna_forSurvival[-remove,]
    
    
    #Now I need to identify normal and tumor samples. this is done using the TCGA barcode (https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode). The two digits at position 14-15 of the barcode will indicate teh sample type, from the link:
    #  "Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29."
    
    # see the values
    table(substr(colnames(rna_forSurvival),14,14))
    # 0      1 
    # 534   72
    
    clinical_forSurvival.t <- as.data.frame(t(clinical_forSurvival))
    #for alive indvivduals, time = days to last follow up; for dead indvivuals, time = days to death
    clinical_forSurvival.t$time <- NA
    clinical_forSurvival.t[clinical_forSurvival.t$vital_status=="Alive",]$time <- as.numeric(as.character(clinical_forSurvival.t[clinical_forSurvival.t$vital_status=="Alive",]$days_to_last_follow_up))
    clinical_forSurvival.t[clinical_forSurvival.t$vital_status=="Dead",]$time <- as.numeric(as.character(clinical_forSurvival.t[clinical_forSurvival.t$vital_status=="Dead",]$days_to_death))
    
    rna_forSurvival.t <- as.data.frame(t(rna_forSurvival))
    
    forSurvival <- cbind(clinical_forSurvival.t,rna_forSurvival.t)
    forSurvival$days_to_death <- as.numeric(as.character(forSurvival$days_to_death))
    forSurvival$vital_status <- as.character(forSurvival$vital_status)
    forSurvival$vital_status <- gsub("Alive","1",forSurvival$vital_status)
    forSurvival$vital_status <- gsub("Dead","2",forSurvival$vital_status)
    forSurvival$vital_status <- as.integer(forSurvival$vital_status)
    
    forSurvival$Customized_group <- NA
    
    #grep("TP53",colnames(forSurvival),value = TRUE)
    #forSurvival[forSurvival[["ENSG00000141510|TP53|2653"]] >= median(forSurvival[["ENSG00000141510|TP53|2653"]]),]$Customized_group <- paste0("ENSG00000141510|TP53|2653", "|High expressed group (50%)")
    #forSurvival[forSurvival[["ENSG00000141510|TP53|2653"]] < median(forSurvival[["ENSG00000141510|TP53|2653"]]),]$Customized_group <- paste0("ENSG00000141510|TP53|2653", "|Low expressed group (50%)")
    
    #CAFs
    forSurvival$Plasma_exclusion_score <-  0.2084667051*forSurvival$`ENSG00000108821|COL1A1|5914`+0.0699843508*forSurvival$`ENSG00000171345|KRT19|1390`+0.3286087154*forSurvival$`ENSG00000168542|COL3A1|5490`-0.0002842341*forSurvival$`ENSG00000162366|PDZK1IP1|838`+ 
      0.4808993020*forSurvival$`ENSG00000164692|COL1A2|5993`-0.0673166255*forSurvival$`ENSG00000172426|RSPH9|2586`+0.0889751549*forSurvival$`ENSG00000141448|GATA6|3624`+0.2194534345*forSurvival$`ENSG00000091986|CCDC80|12301`+
      0.6918664090*forSurvival$`ENSG00000174939|ASPHD1|1816`+0.0339712745*forSurvival$`ENSG00000153495|TEX29|835`+0.0429682274*forSurvival$`ENSG00000115380|EFEMP1|3024`-0.5139937159*forSurvival$`ENSG00000178531|CTXN1|1237`+ 
      0.0243741768*forSurvival$`ENSG00000186832|KRT16|1658`-0.0472012043*forSurvival$`ENSG00000155622|XAGE2|627`+0.1070337255*forSurvival$`ENSG00000006042|TMEM98|4218`+0.1790269326*forSurvival$`ENSG00000060718|COL11A1|7327`+  
      0.2421727418*forSurvival$`ENSG00000175745|NR2F1|3843`+0.5693757555*forSurvival$`ENSG00000205076|LGALS7|565`-0.1713988864*forSurvival$`ENSG00000183569|SERHL2|2096`+0.3891344341*forSurvival$`ENSG00000047936|ROS1|8451`
    
    forSurvival$Plasma_exclusion_score <-  forSurvival$`ENSG00000108821|COL1A1|5914`+forSurvival$`ENSG00000171345|KRT19|1390`+forSurvival$`ENSG00000168542|COL3A1|5490`+forSurvival$`ENSG00000162366|PDZK1IP1|838`+ 
      forSurvival$`ENSG00000164692|COL1A2|5993`+forSurvival$`ENSG00000172426|RSPH9|2586`+forSurvival$`ENSG00000141448|GATA6|3624`+ forSurvival$`ENSG00000091986|CCDC80|12301`+
      forSurvival$`ENSG00000174939|ASPHD1|1816`+forSurvival$`ENSG00000153495|TEX29|835`+forSurvival$`ENSG00000115380|EFEMP1|3024`+forSurvival$`ENSG00000178531|CTXN1|1237`+ 
      forSurvival$`ENSG00000186832|KRT16|1658`+forSurvival$`ENSG00000155622|XAGE2|627`+forSurvival$`ENSG00000006042|TMEM98|4218`+forSurvival$`ENSG00000060718|COL11A1|7327`+  
      forSurvival$`ENSG00000175745|NR2F1|3843`+forSurvival$`ENSG00000205076|LGALS7|565`+forSurvival$`ENSG00000183569|SERHL2|2096`+forSurvival$`ENSG00000047936|ROS1|8451`
    
    
    forSurvival$Plasma_exclusion_score <- forSurvival$`ENSG00000108821|COL1A1|5914`
    forSurvival$Plasma_exclusion_score <- forSurvival$`ENSG00000153563|CD8A|3048`
    forSurvival[forSurvival[["Plasma_exclusion_score"]] >= median(forSurvival[["Plasma_exclusion_score"]]),]$Customized_group <- "High" #paste0("Plasma_exclusion_score", "|High expressed group (50%)")
    forSurvival[forSurvival[["Plasma_exclusion_score"]] < median(forSurvival[["Plasma_exclusion_score"]]),]$Customized_group <- "Low" #paste0("Plasma_exclusion_score", "|Low expressed group (50%)")
    
    #NK cell score
    forSurvival$Plasma_NK_cell_resting_score <- #forSurvival$`ENSG00000111537|IFNG|1211` + forSurvival$`ENSG00000081985|IL12RB2|5443`
      (forSurvival$`ENSG00000148600|CDHR1|6877`+forSurvival$`ENSG00000162676|GFI1|4554`+
      forSurvival$`ENSG00000075234|TTC38|2566`+forSurvival$`ENSG00000129566|TEP1|10775`)/
      (forSurvival$`ENSG00000205810|KLRC3|1025`+
         forSurvival$`ENSG00000081985|IL12RB2|5443`+forSurvival$`ENSG00000105374|NKG7|816`+
         forSurvival$`ENSG00000180739|S1PR5|2338`+forSurvival$`ENSG00000115085|ZAP70|3273`)
    
    forSurvival$Plasma_NK_cell_resting_score <-
      
    #forSurvival$`ENSG00000134545|KLRC1|1600`
    
    forSurvival$Plasma_NK_cell_resting_score <- forSurvival$`ENSG00000172232|AZU1|1713`+forSurvival$`ENSG00000101425|BPI|1845`+forSurvival$`ENSG00000164047|CAMP|745`+forSurvival$`ENSG00000148600|CDHR1|6877`+
      forSurvival$`ENSG00000124469|CEACAM8|2297`+forSurvival$`ENSG00000164821|DEFA4|587`+forSurvival$`ENSG00000197561|ELANE|1028`+forSurvival$`ENSG00000257335|MGAM|9172`+
      forSurvival$`ENSG00000149516|MS4A3|1618`+forSurvival$`ENSG00000086288|NME8|2308`+forSurvival$`ENSG00000166289|PLEKHF1|1699`+forSurvival$`ENSG00000129566|TEP1|10775`+
      forSurvival$`ENSG00000075234|TTC38|2566`+forSurvival$`ENSG00000176293|ZNF135|3349`
    
    forSurvival$Plasma_NK_cell_activated_score <- forSurvival$`ENSG00000213809|KLRK1|2064`
    
    #gammadeltaT 
    forSurvival$Plasma_gammadeltaT_score <- -6.059254e-04*forSurvival$`ENSG00000089012|SIRPG|1732`-5.665177e-04*forSurvival$`ENSG00000183918|SH2D1A|2554`-1.280746e-06*forSurvival$`ENSG00000271503|CCL5|1352`-8.250298e-05*forSurvival$`ENSG00000211751|TRBC1|760`+
      1.162115e-04*forSurvival$`ENSG00000172116|CD8B|4794`-6.216530e-04*forSurvival$`ENSG00000125245|GPR18|1875`-2.667684e-05*forSurvival$`ENSG00000211829|TRDC|720`-2.455918e-03*forSurvival$`ENSG00000174946|GPR171|1810`-
      7.147134e-05*forSurvival$`ENSG00000167286|CD3D|701`-7.103221e-04*forSurvival$`ENSG00000206561|COLQ|6346`-3.057003e-03*forSurvival$`ENSG00000158050|DUSP2|1685`-1.694696e-04*forSurvival$`ENSG00000139187|KLRG1|1874`-
      2.824135e-05*forSurvival$`ENSG00000077984|CST7|891`-7.939934e-04*forSurvival$`ENSG00000089692|LAG3|2587`-4.949101e-05*forSurvival$`ENSG00000145649|GZMA|896`-1.090442e-04*forSurvival$`ENSG00000112303|VNN2|2014`-
      1.743048e-04*forSurvival$`ENSG00000115607|IL18RAP|2773`+7.649479e-05*forSurvival$`ENSG00000113088|GZMK|1506`
    
    #gammadeltaT2
    forSurvival$Plasma_gammadeltaT_score <- -0.053551401*forSurvival$`ENSG00000089012|SIRPG|1732`+0.170103962*forSurvival$`ENSG00000122224|LY9|5083`+0.140654566*forSurvival$`ENSG00000162676|GFI1|4554`-0.042064351*forSurvival$`ENSG00000183918|SH2D1A|2554`+0.002171451*forSurvival$`ENSG00000213658|LAT|2443`
    
    forSurvival$Plasma_gammadeltaT_score <- -forSurvival$`ENSG00000089012|SIRPG|1732`+forSurvival$`ENSG00000122224|LY9|5083`+forSurvival$`ENSG00000162676|GFI1|4554`-forSurvival$`ENSG00000183918|SH2D1A|2554`+forSurvival$`ENSG00000213658|LAT|2443`
    
    #forSurvival$Plasma_gammadeltaT_score <- forSurvival$`ENSG00000172116|CD8B|4794`
    
    
    Group_by <- "Plasma_NK_cell_resting_score"
    forSurvival[forSurvival[[Group_by]] >= median(forSurvival[[Group_by]]),]$Customized_group <- "High" #paste0("Plasma_gammadeltaT_score", "|High expressed group (50%)")
    forSurvival[forSurvival[[Group_by]] < median(forSurvival[[Group_by]]),]$Customized_group <- "Low" #paste0("Plasma_gammadeltaT_score", "|Low expressed group (50%)")
    forSurvival$Customized_group <- factor(forSurvival$Customized_group, levels = c("Low","High"))
    forSurvival_plot <- forSurvival
    #forSurvival_plot <- forSurvival[-which(forSurvival$project_id=="TCGA-STAD"),]
    
    Survival_data <- forSurvival_plot[,c("gender","time","vital_status","project_id","Customized_group")]
    
    Survival_data <- na.omit(Survival_data)
    #Survival_data[Survival_data$time >= 1095,"time"] <- 1095
    #Survival_data[Survival_data$time >= 1095,"vital_status"] <- 1
    fit1 <- surv_fit(Surv(time, vital_status) ~ Customized_group, data = Survival_data)
    plot1 <- ggsurvplot(fit1, data = forSurvival, pval = TRUE, pval.method = TRUE, conf.int = TRUE, 
                        risk.table = "absolute", 
                        palette = "aaas", 
                        ggtheme = theme_survminer(base_size = 16, base_family = "Arial",
                                                  font.main = c(20,"plain","black"), font.tickslab = c(16,"plain","black"),
                                                  font.x = c(24,"plain","black"),font.y = c(24,"plain","black"),
                                                  font.legend = c(16,"plain","black"), legend = "right"))
    plot1
    ggsave(plot = plot1 , filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 6/Survival_plasma_CAFs_STAD_COAD+READ_641.pdf",device = "pdf",width = 4.5,height = 5)
    
    ggscatter(forSurvival, x = "tumor_stage", y = "Plasma_exclusion_score", 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.coeff.args = list(method = "spearman",label.sep = "\n",size = 8))+
      #xlab("Exclusion - γδ T cells - B cells - plasma cells")+
      #ylab("Tumor size (cm)")+
      theme(#legend.position="right",
        legend.position="right",
        panel.grid=element_blank(),
        plot.margin = unit(c(20,20,20,20),"pt"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=20, family = "Arial"),
        axis.text.y = element_text(face="bold",  color="black", size=20, family = "Arial"),
        axis.title.x = element_text(face="bold", color="black", size=16, family = "Arial"),
        axis.title.y = element_text(face="bold",color="black", size=24, family = "Arial"))
    
    forSurvival$Stage_simplified <- NA
    forSurvival[which(forSurvival$tumor_stage=="not reported"),]$Stage_simplified <- "Not reported"
    forSurvival[which(forSurvival$tumor_stage=="--"),]$Stage_simplified <- "Not reported"
    
    forSurvival[which(forSurvival$tumor_stage=="stage i"),]$Stage_simplified <- "Stage I"
    forSurvival[which(forSurvival$tumor_stage=="stage ia"),]$Stage_simplified <- "Stage I"
    forSurvival[which(forSurvival$tumor_stage=="stage ib"),]$Stage_simplified <- "Stage I"
    
    forSurvival[which(forSurvival$tumor_stage=="stage ii"),]$Stage_simplified <- "Stage II"
    forSurvival[which(forSurvival$tumor_stage=="stage iia"),]$Stage_simplified <- "Stage II"
    forSurvival[which(forSurvival$tumor_stage=="stage iib"),]$Stage_simplified <- "Stage II"
    forSurvival[which(forSurvival$tumor_stage=="stage iic"),]$Stage_simplified <- "Stage II"
    
    forSurvival[which(forSurvival$tumor_stage=="stage iii"),]$Stage_simplified <- "Stage III"
    forSurvival[which(forSurvival$tumor_stage=="stage iiia"),]$Stage_simplified <- "Stage III"
    forSurvival[which(forSurvival$tumor_stage=="stage iiib"),]$Stage_simplified <- "Stage III"
    forSurvival[which(forSurvival$tumor_stage=="stage iiic"),]$Stage_simplified <- "Stage III"
    
    forSurvival[which(forSurvival$tumor_stage=="stage iv"),]$Stage_simplified <- "Stage IV"
    forSurvival[which(forSurvival$tumor_stage=="stage iva"),]$Stage_simplified <- "Stage IV"
    forSurvival[which(forSurvival$tumor_stage=="stage ivb"),]$Stage_simplified <- "Stage IV"
    
    forSurvival <- forSurvival[which(forSurvival$Stage_simplified!="Not reported"),]
  
    p<-ggplot(forSurvival[grep("TCGA-COAD|TCGA-READ",forSurvival$project_id),],aes(x=Stage_simplified,y=log10(Plasma_gammadeltaT_score),fill = Stage_simplified))+
      geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
      geom_point(size = 1, position = position_jitterdodge(dodge.width=1.1,jitter.width = 0.8),color = alpha("black",alpha = 0.2))+
      #scale_fill_manual(values=c("CRC"="#FCB514","STAD"="red","HD"="blue","early"="dark green","late"="black"))+
      scale_fill_manual(values=c("Stage I"="#90E0EF","Stage II"="#00B4D8","Stage III"="#03045E","Stage III/IV"="#03045E","Stage IV"="black"))+
      #scale_fill_brewer(palette="Blues") +
      #ylim(0,25)+
      theme_bw()+
      #ylab(colnames(plot)[1])+
      theme(#legend.position="right",
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=36,face="bold"),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(color="black", size=24,angle = 45,hjust = 1),
        axis.text.y = element_text(color="black", size=24),
        axis.title.x = element_text(color="black", size=0),
        axis.title.y = element_text(color="black", size=16))+
      stat_compare_means(comparisons = 
                           list(c("Stage I","Stage III"),c("Stage I","Stage II"),c("Stage I","Stage IV")),
                         method = "wilcox.test",
                         method.args = list(alternative = "greater",paired = TRUE),
                         label = "p.signif",
                         size = 10,
                         vjust = 0.5)
    p
    ggsave(plot=p,filename = "/Users/yuhuan/Desktop/Seafile/Lulab/2022/01.multiomics/Figure 6/Gamma_delat_T_inTCGA.pdf",width=4.40,height = 5.69)
    }
}
