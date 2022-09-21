#!/bin/bash
#SBATCH -J rMATs
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
positive=$2
negative=$3
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/$1"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/$1/${positive}vs${negative}"
mkdir -p ${output}
mkdir -p ${output}/tmp
gtf="/data/taoyuhuan/projects/exOmics_RNA/bin/genome/gtf/gencode.v27.annotation.gtf"

echo "Start call alternative splicing for $1: ${positive} vs ${negative} at `date`"
rmats.py --b1 ${input}/${positive}.txt --b2 ${input}/${negative}.txt --gtf ${gtf} --tmp ${output}/tmp \
         --od ${output} -t paired --libType fr-firststrand --readLength 150 --nthread 16 --variable-read-length --allow-clipping

#SMARTer v1: forward/fr-secondstrand
#SMARTer v2: reverse/fr-firststrand
#multiomics_plasma: reverse/fr-firststrand
#multiomics_PBMC:reverse/fr-firststrand
#multiomics_tissue: forward/fr-firststrand
