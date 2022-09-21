#!/bin/bash
#SBATCH -J VEP
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#input_dir=$1
#output_dir=$2
#sample=$3

dataset=$1
input_dir="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/SNP_filtered"
output_dir="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/SNP_VEP"
mkdir -p ${output_dir}

#for vcf in $(ls ${input_dir} | grep -E "mix.-.-pico") #| grep ${sample})
for vcf in $(ls ${input_dir})
do
echo "Start at `date`"
vep -i ${input_dir}/${vcf} --tab --species human --cache --dir_cache /data/taoyuhuan/reference/VEP --everything --stats_text --stats_file ${output_dir}/${vcf%.*} -o ${output_dir}/${vcf%.*}.VEP.txt --fork 16 --pick --force_overwrite
done
