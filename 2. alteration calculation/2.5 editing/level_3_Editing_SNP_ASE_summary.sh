#!/bin/bash
#SBATCH -J Summary_Editing_SNP_ASE
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --exclude=biot[03]
#SBATCH --ntasks=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
tmp="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/tmp"
mkdir -p ${tmp}
input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/bam"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}"
mkdir -p ${output}
mkdir -p ${output}/bam_RG
mkdir -p ${output}/bam_Intronspanning
mkdir -p ${output}/REDIportal
mkdir -p ${output}/ASE
mkdir -p ${output}/ASE/dbSNP
mkdir -p ${output}/ASE/COSMIC
mkdir -p ${output}/SNP
mkdir -p ${output}/SNP/tmp
mkdir -p ${output}/matrix
ref="/data/taoyuhuan/projects/exOmics_RNA/bin"
log="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/log"
mkdir -p ${log}
mkdir -p ${log}/ReadGroup
mkdir -p ${log}/REDIportal
mkdir -p ${log}/ASE
mkdir -p ${log}/ASE/dbSNP
mkdir -p ${log}/ASE/COSMIC
mkdir -p ${log}/SNP


#summarize matrix
Rscript /data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/SNP_AF.R -i ${output}/SNP -o ${output}/matrix

Rscript /data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/Editing_ratio.R -i ${output}/REDIportal -o ${output}/matrix/Editing_ratio.txt

Rscript /data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/ASE_AE.R -i ${output}/ASE/dbSNP -o ${output}/matrix/ASE_dbSNP_${dataset}.txt
Rscript /data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/ASE_AE.R -i ${output}/ASE/COSMIC -o ${output}/matrix/ASE_COSMIC_${dataset}.txt
