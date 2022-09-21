#! /usr/bin/env Rscript
library(reshape2)
library(dplyr)
input_dir <- "/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/multiomics_paired/SNP_VEP"
output_dir <- "/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/multiomics_paired/SNP_VEP/Summary"

dir.create(output_dir,recursive = TRUE)

samples <- grep("VEP|Summary",dir(input_dir),value = TRUE,invert = TRUE)
summary_all <- read.table(paste0(input_dir,"/",samples[1]),fill = TRUE,sep = "\t",row.names=NULL)

Variant_class <- grep("[Variant classes]",summary_all[,1],fixed = TRUE)
Severe_consequence <- grep("[Consequences (most severe)]",summary_all[,1],fixed = TRUE)
All_consequence <- grep("[Consequences (all)]",summary_all[,1],fixed = TRUE)
Coding_consequence <- grep("[Coding consequences]",summary_all[,1],fixed = TRUE)
SIFT_summary <- grep("[SIFT summary]",summary_all[,1],fixed = TRUE)
PolyPhen_summary <- grep("[PolyPhen summary]",summary_all[,1],fixed = TRUE)
Variants_by_chromosome <- grep("[Variants by chromosome]",summary_all[,1],fixed = TRUE)

summary_Variant_class <- summary_all[(Variant_class+1):(Severe_consequence-1),1,drop=FALSE]
summary_Severe_consequence <- summary_all[(Severe_consequence+1):(All_consequence-1),1,drop=FALSE]
summary_All_consequence <- summary_all[(All_consequence+1):(Coding_consequence-1),1,drop=FALSE]
summary_Coding_consequence <- summary_all[(Coding_consequence+1):(SIFT_summary-1),1,drop=FALSE]
summary_SIFT_summary <- summary_all[(SIFT_summary+1):(PolyPhen_summary-1),1,drop=FALSE]
summary_PolyPhen_summary <- summary_all[(PolyPhen_summary+1):(Variants_by_chromosome-1),1,drop=FALSE]
colnames(summary_Variant_class) <- c("Variant")
colnames(summary_Severe_consequence) <- c("Variant")
colnames(summary_All_consequence) <- c("Variant")
colnames(summary_Coding_consequence) <- c("Variant")
colnames(summary_SIFT_summary) <- c("Variant")
colnames(summary_PolyPhen_summary) <- c("Variant")

for(sample in samples){
print(sample)
summary <- read.table(paste0(input_dir,"/",sample),fill = TRUE,sep = "\t",row.names=NULL)

Variant_class <- grep("[Variant classes]",summary[,1],fixed = TRUE)
Severe_consequence <- grep("[Consequences (most severe)]",summary[,1],fixed = TRUE)
All_consequence <- grep("[Consequences (all)]",summary[,1],fixed = TRUE)
Coding_consequence <- grep("[Coding consequences]",summary[,1],fixed = TRUE)
SIFT_summary <- grep("[SIFT summary]",summary[,1],fixed = TRUE)
PolyPhen_summary <- grep("[PolyPhen summary]",summary[,1],fixed = TRUE)
Variants_by_chromosome <- grep("[Variants by chromosome]",summary[,1],fixed = TRUE)

summary_Variant_class_sample <- summary[(Variant_class+1):(Severe_consequence-1),1:2,drop=FALSE]
summary_Severe_consequence_sample <- summary[(Severe_consequence+1):(All_consequence-1),1:2,drop=FALSE]
summary_All_consequence_sample <- summary[(All_consequence+1):(Coding_consequence-1),1:2,drop=FALSE]
summary_Coding_consequence_sample <- summary[(Coding_consequence+1):(SIFT_summary-1),1:2,drop=FALSE]
summary_SIFT_summary_sample <- summary[(SIFT_summary+1):(PolyPhen_summary-1),1:2,drop=FALSE]
summary_PolyPhen_summary_sample <- summary[(PolyPhen_summary+1):(Variants_by_chromosome-1),1:2,drop=FALSE]

sample <- gsub(".rmEDIT.SNP.filtered.txt","",fixed = TRUE,sample)
colnames(summary_Variant_class_sample) <- c("Variant",sample)
colnames(summary_Severe_consequence_sample) <- c("Variant",sample)
colnames(summary_All_consequence_sample) <- c("Variant",sample)
colnames(summary_Coding_consequence_sample) <- c("Variant",sample)
colnames(summary_SIFT_summary_sample) <- c("Variant",sample)
colnames(summary_PolyPhen_summary_sample) <- c("Variant",sample)

summary_Variant_class <- full_join(summary_Variant_class,summary_Variant_class_sample,by = c("Variant"="Variant"))
summary_Severe_consequence <- full_join(summary_Severe_consequence,summary_Severe_consequence_sample,by = c("Variant"="Variant"))
summary_All_consequence <- full_join(summary_All_consequence,summary_All_consequence_sample, by = c("Variant"="Variant"))
summary_Coding_consequence <- full_join(summary_Coding_consequence,summary_Coding_consequence_sample, by = c("Variant"="Variant"))
summary_SIFT_summary <- full_join(summary_SIFT_summary,summary_SIFT_summary_sample, by = c("Variant"="Variant"))
summary_PolyPhen_summary <- full_join(summary_PolyPhen_summary,summary_PolyPhen_summary_sample, by = c("Variant"="Variant"))
}

summary_Variant_class[is.na(summary_Variant_class)] <- 0
summary_Severe_consequence[is.na(summary_Severe_consequence)] <- 0
summary_All_consequence[is.na(summary_All_consequence)] <- 0
summary_Coding_consequence[is.na(summary_Coding_consequence)] <- 0
summary_SIFT_summary[is.na(summary_SIFT_summary)] <- 0
summary_PolyPhen_summary[is.na(summary_PolyPhen_summary)] <- 0

write.table(summary_Variant_class,paste0(output_dir,"/Variant_class.txt"),quote=FALSE,sep = "\t",row.names = FALSE)
write.table(summary_Severe_consequence,paste0(output_dir,"/Severe_consequence.txt"),quote=FALSE,sep = "\t",row.names = FALSE)
write.table(summary_All_consequence,paste0(output_dir,"/All_consequence.txt"),quote=FALSE,sep = "\t",row.names = FALSE)
write.table(summary_Coding_consequence,paste0(output_dir,"/Coding_consequence.txt"),quote=FALSE,sep = "\t",row.names = FALSE)
write.table(summary_SIFT_summary,paste0(output_dir,"/SIFT_summary.txt"),quote=FALSE,sep = "\t",row.names = FALSE)
write.table(summary_PolyPhen_summary,paste0(output_dir,"/PolyPhen_summary.txt"),quote=FALSE,sep = "\t",row.names = FALSE)
