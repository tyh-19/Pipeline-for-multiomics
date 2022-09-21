#!/bin/bash
#SBATCH -J Alt.promoter_usage
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

sample=$1
Rscript Promoter_usage.R -n /data/taoyuhuan/projects/exOmics_RNA/level_3_Alt.promoter/multiomics_paired/AlternativePromoter_multiomics_paired_normalized.txt -t /data/taoyuhuan/projects/exOmics_RNA/level_3_Alt.promoter/multiomics_paired/TPM-by-promoter_multiomics_paired.txt -o /data/taoyuhuan/projects/exOmics_RNA/level_3_Alt.promoter/multiomics_paired/promoter_usage -s ${sample}
